/**
 * @file hdf5_r6_minimal.cpp
 * @brief Minimal R6 wrappers for HDF5 dataset operations
 *
 * PRODUCTION VERSION - DEBUG outputs removed for performance
 * FIXED: Accepts 3 arguments (filename, group, dataset) as expected by R code
 *
 * Live-pointer registry (live_ptrs):
 * ------------------------------------
 * Every C++ hdf5Dataset object created by rcpp_hdf5dataset_open() is
 * inserted into a process-wide unordered_map<void*, SEXP>.  The entry is
 * removed — and the object deleted — exactly once, either by an explicit
 * rcpp_hdf5dataset_close() call or by the XPtr GC finalizer, whichever
 * fires first.  Any subsequent attempt to delete the same address (stale
 * XPtr, recycled address, double-close) is silently ignored.
 *
 * Storing SEXP alongside the raw pointer allows rcpp_hdf5_close_at_paths()
 * to clear the XPtr before deleting the C++ object, so that any R6
 * objects whose output datasets are overwritten see is_valid() == FALSE
 * immediately rather than holding a dangling pointer.
 *
 * This eliminates the family of double-free / use-after-free crashes that
 * can occur when:
 *   - R's allocator recycles a freed hdf5Dataset address for a new object
 *     while a stale XPtr still holds the old (now-invalid) raw address.
 *   - An on.exit handler or R6 finalizer fires after the object was already
 *     deleted via an explicit close() call.
 *   - A previous result object (e.g. eigen vectors) is still open when the
 *     same output path is overwritten by a new call with overwrite=TRUE.
 */

#include <BigDataStatMeth.hpp>
#include <unordered_map>
#include <mutex>

// ---------------------------------------------------------------------------
// Live-pointer registry
// ---------------------------------------------------------------------------
namespace {
    /**
     * @brief Process-wide map of currently-alive hdf5Dataset raw pointers
     *        to their corresponding R external pointer SEXP.
     *
     * Access is intentionally unsynchronized: all R/Rcpp code is
     * single-threaded from R's perspective (GC, user code, and C entry
     * points all run on the same thread).  The static local is
     * initialised once on first use.
     *
     * The SEXP values are kept alive by the R6 objects that own them —
     * we store them without PROTECT because the GC cannot collect an SEXP
     * that is still referenced by a live R object.
     */
    std::unordered_map<void*, SEXP>& live_ptrs() {
        static std::unordered_map<void*, SEXP> m;
        return m;
    }

    /** Register a newly created object with its XPtr SEXP. */
    inline void register_ptr(void* p, SEXP sexp) {
        live_ptrs()[p] = sexp;
    }

    /**
     * @brief Attempt to take ownership of a pointer for deletion.
     *
     * Removes the entry from the registry.  If found, also clears the
     * external pointer so that R_ExternalPtrAddr() returns NULL for any
     * remaining R references.
     *
     * @return true  — pointer was live; caller must delete the C++ object.
     * @return false — pointer was not live; caller must NOT delete.
     */
    inline bool claim_ptr(void* p) {
        auto it = live_ptrs().find(p);
        if (it == live_ptrs().end()) return false;
        SEXP sexp = it->second;
        live_ptrs().erase(it);
        // Clear the XPtr so any remaining R references see NULL immediately
        if (sexp != R_NilValue && TYPEOF(sexp) == EXTPTRSXP)
            R_ClearExternalPtr(sexp);
        return true;
    }
} // anonymous namespace



//' Close all open HDF5Dataset objects and HDF5 handles
//'
//' @description
//' Closes all C++ \code{hdf5Dataset} objects tracked in the live-pointer
//' registry and then calls \code{BigDataStatMeth::closeAllHDF5Handles()}
//' to close any remaining HDF5 handles at the C library level (files,
//' datasets, groups, datatypes, attributes) that were not tracked by
//' the registry. Equivalent in effect to \code{rhdf5::h5closeAll()}.
//'
//' Called automatically from \code{.onUnload()} when the package is
//' unloaded. Can also be called manually for diagnostic purposes via
//' \code{BigDataStatMeth:::rcpp_hdf5_close_all_registry()}.
//'
//' @return \code{NULL} invisibly.
//'
//' @keywords internal
// [[Rcpp::export]]
SEXP rcpp_hdf5_close_all_registry() {
    
    try {
        std::vector<void*> all_ptrs;
        for (auto& kv : live_ptrs())
            all_ptrs.push_back(kv.first);
        for (void* vp : all_ptrs)
            if (claim_ptr(vp))
                delete static_cast<BigDataStatMeth::hdf5Dataset*>(vp);
            // closeAllHDF5Handles() deliberately NOT called here:
            // it is reserved for .onUnload() only, where library state
            // does not matter. Calling it here corrupts HDF5 pre-defined types.
    } catch (...) {}
    return R_NilValue;
}



//' Open HDF5 dataset and return external pointer (R6 wrapper)
//'
//' @param filename Path to HDF5 file
//' @param group Group path (e.g., "data" or "/data")
//' @param dataset Dataset name within the group (e.g., "matrix")
//'
//' @return External pointer to hdf5Dataset object
//'
//' @keywords internal
// [[Rcpp::export]]
SEXP rcpp_hdf5dataset_open(std::string filename, std::string group, std::string dataset)
{
    try {
        H5::Exception::dontPrint();

        auto* ds = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
        ds->openDataset();

        if (ds->getDatasetptr() == nullptr) {
            delete ds;
            std::string msg = "Failed to open dataset: " + group + "/" + dataset;
            Rf_error("%s", msg.c_str());
        }

        SEXP ptr = R_MakeExternalPtr(ds, R_NilValue, R_NilValue);

        // Register BEFORE handing to R — must happen before any GC can run.
        register_ptr(static_cast<void*>(ds), ptr);

        R_RegisterCFinalizerEx(ptr, [](SEXP p) {
            // GC finalizer: only deletes if the pointer is still in the
            // live registry.  If rcpp_hdf5dataset_close() already ran,
            // claim_ptr() returns false and we do nothing — no double-free.
            // claim_ptr() also calls R_ClearExternalPtr internally, but
            // since we are inside the finalizer the XPtr is already being
            // collected — that call is a harmless no-op.
            auto* obj = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(p));
            if (obj && claim_ptr(static_cast<void*>(obj))) {
                delete obj;
            }
        }, (Rboolean)TRUE);

        return ptr;

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error opening dataset: %s", e.what());
    }
}


//' Get dimensions of HDF5 dataset (R6 wrapper)
//'
//' @param ptr_sexp External pointer to hdf5Dataset
//'
//' @return Integer vector c(nrows, ncols)
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_hdf5dataset_dim(SEXP ptr_sexp)
{
    try {
        auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_sexp));

        if (ds == nullptr)
            Rf_error("Invalid external pointer");
        if (!ds->isOpen())
            Rf_error("Dataset is closed");

        return Rcpp::IntegerVector::create(
            static_cast<int>(ds->nrows_r()),
            static_cast<int>(ds->ncols_r())
        );

    } catch (std::exception& e) {
        Rf_error("Error getting dimensions: %s", e.what());
    }
}


//' Get dataset information (R6 wrapper)
//'
//' @param ptr_sexp External pointer to hdf5Dataset
//'
//' @return Named list with filename, group, dataset, datatype
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_info(SEXP ptr_sexp)
{
    try {
        auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_sexp));

        if (ds == nullptr)
            Rf_error("Invalid external pointer");
        if (!ds->isOpen())
            Rf_error("Dataset is closed");

        return Rcpp::List::create(
            Rcpp::Named("filename") = ds->getFullPath(),
            Rcpp::Named("group")    = ds->getGroup(),
            Rcpp::Named("dataset")  = ds->getDatasetName(),
            Rcpp::Named("datatype") = ds->getDatatype()
        );

    } catch (std::exception& e) {
        Rf_error("Error getting info: %s", e.what());
    }
}


//' Check if dataset is valid and open (R6 wrapper)
//'
//' @param ptr_sexp External pointer to hdf5Dataset
//'
//' @return Logical: TRUE if valid and open, FALSE otherwise
//'
//' @keywords internal
// [[Rcpp::export]]
bool rcpp_hdf5dataset_is_valid(SEXP ptr_sexp)
{
    try {
        auto* ptr = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_sexp));

        if (ptr == nullptr) return false;
        if (!ptr->isOpen()) return false;

        try {
            hid_t id = ptr->getDatasetptr()->getId();
            return H5Iis_valid(id) > 0;
        } catch (...) {
            return false;
        }

    } catch(...) {
        return false;
    }
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_read_dimnames
// ---------------------------------------------------------------------------

//' Read dimension names (rownames / colnames) from an HDF5 dataset
//'
//' @description
//' Reads the row and column names stored alongside an HDF5 dataset following
//' the BigDataStatMeth convention:
//' \itemize{
//'   \item rownames stored at \code{group/.<dataset>_dimnames/1}
//'   \item colnames stored at \code{group/.<dataset>_dimnames/2}
//' }
//' When a component has not been written an empty \code{character(0)} is
//' returned for it.  The function uses \code{BigDataStatMeth::hdf5Dims} in
//' read mode (\code{bWrite = false}) so no data on disk is modified.
//'
//' @param ptr_sexp External pointer (SEXP) to an open \code{hdf5Dataset}
//'   object managed by the R6 class.
//'
//' @return Named list with two \code{character} elements:
//' \describe{
//'   \item{rownames}{Row names, or \code{character(0)} if absent}
//'   \item{colnames}{Column names, or \code{character(0)} if absent}
//' }
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_read_dimnames(SEXP ptr_sexp)
{
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("rownames") = Rcpp::CharacterVector(0),
        Rcpp::Named("colnames") = Rcpp::CharacterVector(0));

    try {
        H5::Exception::dontPrint();

        auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(
                       R_ExternalPtrAddr(ptr_sexp));

        if (ds == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!ds->isOpen())
            throw std::runtime_error("Dataset is closed");

        result = ds->readDimnames();

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error reading dimnames: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error reading dimnames: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error reading dimnames: %s", e.what());
    }

    return result;
}

//' Write dimension names through the R6 dataset handle
//'
//' @description
//' Writes row and/or column names for an HDF5 dataset using the existing
//' open file handle managed by the R6 object.  Unlike
//' \code{bdWrite_hdf5_dimnames()}, this function operates through
//' \code{hdf5Dataset::writeDimnames()} so the long-lived R6 handle sees
//' the changes immediately - no metadata cache staleness.
//'
//' @param ptr_sexp External pointer (SEXP) to an open \code{hdf5Dataset}.
//' @param rownames Character vector of row names. Use \code{character(0)}
//'   to skip writing row names.
//' @param colnames Character vector of column names. Use \code{character(0)}
//'   to skip writing column names.
//'
//' @return \code{NULL} invisibly.
//'
//' @keywords internal
// [[Rcpp::export]]
SEXP rcpp_hdf5dataset_write_dimnames(SEXP ptr_sexp,
                                     Rcpp::StringVector rownames,
                                     Rcpp::StringVector colnames)
{
    try {
        H5::Exception::dontPrint();

        auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_sexp));

        if (ds == nullptr)
            Rf_error("Invalid external pointer");
        if (!ds->isOpen())
            Rf_error("Dataset is closed");

        ds->writeDimnames(rownames, colnames);

    } catch (std::exception& e) {
        Rf_error("Error writing dimnames: %s", e.what());
    }

    return R_NilValue;
}


//' Close and destroy an HDF5 dataset handle immediately.
//'
//' Uses the live-pointer registry to prevent double-free: if the pointer
//' is no longer in the registry (already closed by close() or GC), this
//' is a safe no-op.  Clears the external pointer so the GC finalizer
//' becomes a no-op too.
//'
//' @param ptr_sexp External pointer to hdf5Dataset
//' @keywords internal
// [[Rcpp::export]]
void rcpp_hdf5dataset_close(SEXP ptr_sexp) {
    // Guard: accept only a proper external pointer (not R NULL / NILSXP).
    if (ptr_sexp == R_NilValue || TYPEOF(ptr_sexp) != EXTPTRSXP) return;
    
    // Guard: si R_ClearExternalPtr ya limpió este SEXP, addr es NULL
    if (R_ExternalPtrAddr(ptr_sexp) == nullptr) return;

    auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(
        R_ExternalPtrAddr(ptr_sexp));

    // claim_ptr() removes from registry AND calls R_ClearExternalPtr internally.
    // We then delete the C++ object.
    if (ds && claim_ptr(static_cast<void*>(ds))) {
        delete ds;
    }
}


//' Close all live HDF5Matrix handles pointing to specific dataset paths.
//'
//' @description
//' Scans the live-pointer registry for any open \code{hdf5Dataset} objects
//' that match the given \code{filename} and any of the \code{paths}.
//' Each matching object is closed and its external pointer cleared, so
//' that any R6 \code{HDF5Matrix} objects holding those pointers will
//' return \code{FALSE} from \code{is_valid()} immediately.
//'
//' This is called automatically by R6 methods that use
//' \code{overwrite = TRUE} (e.g. \code{$eigen()}, \code{$svd()},
//' \code{$qr()}, \code{$chol()}, \code{$prcomp()}) to ensure that
//' previous result objects are safely invalidated before the HDF5 datasets
//' they reference are deleted and recreated.
//'
//' @param filename  Canonical filesystem path to the HDF5 file.
//' @param paths     Character vector of HDF5-internal paths
//'   (e.g. \code{c("EIGEN/sym/values", "EIGEN/sym/vectors")}).
//'
//' @return \code{NULL} invisibly.
//'
//' @keywords internal
// [[Rcpp::export]]
SEXP rcpp_hdf5_close_at_paths(std::string filename, Rcpp::StringVector paths)
{
    try {
        H5::Exception::dontPrint();

        // Build a fast lookup set of target paths
        std::unordered_set<std::string> path_set;
        for (int i = 0; i < paths.size(); ++i)
            path_set.insert(Rcpp::as<std::string>(paths[i]));

        // Collect matching raw pointers first (don't modify map while iterating)
        std::vector<void*> to_close;
        for (auto& kv : live_ptrs()) {
            auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(kv.first);
            try {
                // Match on canonical filename AND group/dataset path
                if (ds->getFullPath() != filename) continue;
                std::string ds_path = ds->getGroup() + "/" + ds->getDatasetName();
                if (path_set.count(ds_path) > 0)
                    to_close.push_back(kv.first);
            } catch (...) {
                // Defensive: skip any pointer that can't be inspected
            }
        }

        // claim_ptr clears the XPtr and removes from registry; then delete
        for (void* vp : to_close) {
            if (claim_ptr(vp))
                delete static_cast<BigDataStatMeth::hdf5Dataset*>(vp);
        }

    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5_close_at_paths: %s", e.what());
    }

    return R_NilValue;
}


//' Close all HDF5 handles for a specific file (R6 wrapper)
//'
//' @description
//' Closes all C++ objects tracked in the live-pointer registry that
//' belong to \code{filename}, then closes any remaining HDF5 handles
//' for that file at the HDF5 C library level.
//'
//' @param filename Absolute path to the HDF5 file (use
//'   \code{normalizePath()} in R before calling).
//'
//' @keywords internal
// [[Rcpp::export]]
void rcpp_hdf5_close_file_handles(std::string filename) {
     try {
         // 1. Close tracked C++ objects for this file
         std::vector<void*> to_close;
         for (auto& kv : live_ptrs()) {
             auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(kv.first);
             try {
                 if (ds && ds->getFullPath() == filename)
                     to_close.push_back(kv.first);
             } catch (...) {}
         }
         for (void* vp : to_close)
             if (claim_ptr(vp))
                 delete static_cast<BigDataStatMeth::hdf5Dataset*>(vp);
             
             // 2. Close any remaining HDF5 handles for this file at C library level
             BigDataStatMeth::closeHDF5HandlesForFile(filename);
             
     } catch (...) {}
}
