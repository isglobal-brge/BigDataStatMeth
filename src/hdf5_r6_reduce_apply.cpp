/**
 * @file hdf5_r6_reduce_apply.cpp
 * @brief R6 wrappers for group-reduce and apply-function operations on HDF5.
 *
 * ARCHITECTURAL RULE (roadmap section 1.4):
 *   Calls header functions directly, never Rcpp-exported bd* symbols.
 *
 * Functions provided:
 *   rcpp_hdf5dataset_reduce         -- reduce all datasets in a group with
 *                                     "+" or "-" (RcppReduce_dataset_hdf5)
 *   rcpp_hdf5dataset_apply_function -- apply an operation to a list of
 *                                     datasets (RcppApplyFunctionHdf5)
 *
 * Both operate on entire HDF5 groups (lists of datasets), not on a single
 * dataset reference, so they take strings (filename, group) rather than
 * external pointers.
 *
 * @note Uses Rf_error() in catch blocks to avoid macOS ARM64 double-free.
 */

#include <BigDataStatMeth.hpp>


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_reduce
// ---------------------------------------------------------------------------

/**
 * @brief Reduce all datasets in a group with a binary operation (R6 wrapper).
 *
 * Applies a reduction operation ("+" or "-") across all datasets found in
 * the given group, accumulating the result into a single output dataset.
 * Delegates to BigDataStatMeth::RcppReduce_dataset_hdf5().
 *
 * @param filename     Path to the HDF5 file.
 * @param group        Group containing the datasets to reduce.
 * @param out_group    Output group path (default "REDUCED").
 * @param out_dataset  Output dataset name (default "reduced").
 * @param func         Reduction function: "+" (default) or "-".
 * @param overwrite    Overwrite existing output dataset (default false).
 * @param remove_input Remove input datasets after reduction (default false).
 * @return Named list with elements "filename", "path", and "func".
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_reduce(std::string filename,
                                    std::string group,
                                    std::string out_group    = "REDUCED",
                                    std::string out_dataset  = "reduced",
                                    std::string func         = "+",
                                    bool        overwrite    = false,
                                    bool        remove_input = false)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "",
        Rcpp::Named("path")     = "",
        Rcpp::Named("func")     = func);

    try {
        H5::Exception::dontPrint();
        
        if (func != "+" && func != "-")
            throw std::runtime_error(
                    "rcpp_hdf5dataset_reduce: func must be \"+\" or \"-\"");
        
        // Get canonical filename
        std::unique_ptr<BigDataStatMeth::hdf5File> f(
                new BigDataStatMeth::hdf5File(filename, false));
        const std::string full_path = f->getFullPath();
        f.reset();
        
        // List all datasets in the group
        std::unique_ptr<BigDataStatMeth::hdf5File> fList(
                new BigDataStatMeth::hdf5File(full_path, false));
        fList->openFile("r");
        Rcpp::StringVector dsnames = fList->getDatasetNames(group, "", "");
        fList.reset();
        
        int ndatasets = dsnames.size();
        if (ndatasets == 0)
            throw std::runtime_error("rcpp_hdf5dataset_reduce: no datasets found in group");
        
        std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
        
        // Read first dataset to initialise accumulator
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsFirst(
                new BigDataStatMeth::hdf5Dataset(
                        full_path, group + "/" + std::string(dsnames[0]), false));
        dsFirst->openDataset();
        if (!dsFirst->getDatasetptr())
            throw std::runtime_error("rcpp_hdf5dataset_reduce: cannot open first dataset");
        
        const hsize_t nr = dsFirst->nrows();   // HDF5-native rows
        const hsize_t nc = dsFirst->ncols();   // HDF5-native cols
        const int comp   = static_cast<int>(dsFirst->getCompressionLevel());
        
        std::vector<double> accum(nr * nc);
        dsFirst->readDatasetBlock({0, 0}, {nr, nc}, stride, block, accum.data());
        dsFirst.reset();
        
        // Accumulate remaining datasets element-wise in HDF5-native space
        for (int i = 1; i < ndatasets; ++i) {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsCur(
                    new BigDataStatMeth::hdf5Dataset(
                            full_path, group + "/" + std::string(dsnames[i]), false));
            dsCur->openDataset();
            if (!dsCur->getDatasetptr())
                throw std::runtime_error("rcpp_hdf5dataset_reduce: cannot open dataset");
            if (dsCur->nrows() != nr || dsCur->ncols() != nc)
                throw std::runtime_error("rcpp_hdf5dataset_reduce: dimension mismatch");
            
            std::vector<double> cur(nr * nc);
            dsCur->readDatasetBlock({0, 0}, {nr, nc}, stride, block, cur.data());
            dsCur.reset();
            
            if (func == "+") {
                for (hsize_t k = 0; k < nr * nc; ++k) accum[k] += cur[k];
            } else {
                for (hsize_t k = 0; k < nr * nc; ++k) accum[k] -= cur[k];
            }
        }
        
        // Write result — same HDF5-native dimensions as inputs (nc, nr for R r×c)
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsOut(
                new BigDataStatMeth::hdf5Dataset(full_path, out_group, out_dataset, overwrite));
        dsOut->setCompressionLevel(comp);
        // createDataset(R_rows, R_cols): nr=HDF5-rows=R_cols, nc=HDF5-cols=R_rows
        dsOut->createDataset(nc, nr, "real");
        dsOut->writeDatasetBlock(accum, {0, 0}, {nr, nc}, stride, block);
        
        lst["filename"] = full_path;
        lst["path"]     = out_group + "/" + out_dataset;

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_reduce (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_reduce (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_reduce: %s", e.what());
    }

    return lst;
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_apply_function
// ---------------------------------------------------------------------------

/**
 * @brief Apply a statistical/algebraic function to a list of HDF5 datasets
 *        (R6 wrapper).
 *
 * Applies func to each dataset listed in datasets (found in group) and
 * writes results to out_group.
 * Delegates to BigDataStatMeth::RcppApplyFunctionHdf5().
 *
 * Valid func values: "QR", "CrossProd", "tCrossProd", "invChol",
 * "blockmult", "CrossProd_double", "tCrossProd_double", "solve",
 * "normalize", "sdmean", "descChol".
 *
 * @param filename     Path to the HDF5 file.
 * @param group        Group containing the input dataset(s).
 * @param datasets     Character vector of dataset names to process.
 * @param out_group    Output group path (default "APPLIED").
 * @param func         Function to apply (see description above).
 * @param b_group      Group for second datasets (e.g. for "blockmult"). "" = none.
 * @param b_datasets   Character vector of second dataset names. NULL = none.
 * @param overwrite    Overwrite existing output datasets (default false).
 * @param transp_a     Transpose input dataset A (default false).
 * @param transp_b     Transpose dataset B (default false).
 * @param full_matrix  Return full symmetric matrix where applicable (default false).
 * @param byrows       Apply operation by rows (default false).
 * @param threads      Integer or NULL; OpenMP threads.
 * @return Named list with elements "filename", "out_group", "func", "datasets".
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_apply_function(
    std::string          filename,
    std::string          group,
    Rcpp::StringVector   datasets,
    std::string          out_group  = "APPLIED",
    std::string          func       = "QR",
    std::string          b_group    = "",
    Rcpp::Nullable<Rcpp::StringVector> b_datasets = R_NilValue,
    bool                 overwrite  = false,
    bool                 transp_a   = false,
    bool                 transp_b   = false,
    bool                 full_matrix = false,
    bool                 byrows     = false,
    Rcpp::Nullable<int>  threads    = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename")  = "",
        Rcpp::Named("out_group") = out_group,
        Rcpp::Named("func")      = func,
        Rcpp::Named("datasets")  = datasets);

    try {
        H5::Exception::dontPrint();

        // Validate func
        Rcpp::CharacterVector valid_funcs = {
            "QR", "CrossProd", "tCrossProd",
            "invChol", "blockmult", "CrossProd_double", "tCrossProd_double",
            "solve", "normalize", "sdmean", "descChol"};
        bool found = false;
        for (int i = 0; i < valid_funcs.size(); ++i)
            if (Rcpp::as<std::string>(valid_funcs[i]) == func) { found = true; break; }
        if (!found)
            throw std::runtime_error(
                "rcpp_hdf5dataset_apply_function: unknown func '" + func + "'");

        // Build canonical filename
        std::unique_ptr<BigDataStatMeth::hdf5File> f(
            new BigDataStatMeth::hdf5File(filename, false));
        const std::string full_path = f->getFullPath();

        // Rcpp::Nullable has no operator= from CharacterVector -- construct inline
        Rcpp::Nullable<Rcpp::CharacterVector> n_bgroup(
            b_group.empty()
                ? R_NilValue
                : static_cast<SEXP>(Rcpp::CharacterVector::create(b_group)));

        // Delegate (Rule 1.4)
        // BigDataStatMeth::RcppApplyFunctionHdf5(
        //     full_path, group, datasets, out_group, func,
        //     n_bgroup, b_datasets,
        //     overwrite, transp_a, transp_b, full_matrix, byrows,
        //     threads);
        
        BigDataStatMeth::RcppApplyFunctionHdf5(
            full_path, group, datasets, out_group, func,
            n_bgroup, b_datasets,
            Rcpp::wrap(overwrite),
            Rcpp::wrap(transp_a),
            Rcpp::wrap(transp_b),
            Rcpp::wrap(full_matrix),
            Rcpp::wrap(byrows),
            threads);     
        

        lst["filename"] = full_path;
        
        // Rename outputs: RcppApplyFunctionHdf5 names them <dataset>,
        // convention expects <func>_<dataset>.
        // Uses renameElement() from hdf5Utilities.hpp via a single file handle.
        {
            std::unique_ptr<BigDataStatMeth::hdf5File> fRen(
                    new BigDataStatMeth::hdf5File(full_path, false));
            fRen->openFile("rw");
            for (int i = 0; i < datasets.size(); ++i) {
                std::string src = out_group + "/" + Rcpp::as<std::string>(datasets[i]);
                std::string dst = out_group + "/" + func + "_" +
                    Rcpp::as<std::string>(datasets[i]);
                try {
                    BigDataStatMeth::renameElement(fRen->getFileptr(), src, dst);
                } catch (...) {}  // dataset may not exist (e.g. sdmean creates mean./sd.)
            }
        }

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_apply_function (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_apply_function (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_apply_function: %s", e.what());
    }

    return lst;
}
