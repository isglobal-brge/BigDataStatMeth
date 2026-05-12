/**
 * @file hdf5_r6_multiply.cpp
 * @brief R6 wrappers for matrix algebra operations
 *
 * @details
 * Provides Rcpp wrappers for three block-wise matrix algebra operations
 * used by the HDF5Matrix R6/S3 interface:
 *
 *  - \c rcpp_hdf5dataset_multiply   : general matrix product  (A \%*\% B)
 *  - \c rcpp_hdf5dataset_crossprod  : cross product           (t(A) \%*\% B)
 *  - \c rcpp_hdf5dataset_tcrossprod : transposed cross product (A \%*\% t(B))
 *
 * All three functions follow the same pattern as the existing \c bdblockmult_hdf5,
 * \c bdCrossprod_hdf5 and \c bdtCrossprod_hdf5 wrappers:
 *  - inputs are external pointers (\c SEXP, accessed via R_ExternalPtrAddr)
 *  - ALL logic including validation is inside the try{} block
 *  - validation throws std::runtime_error (not Rcpp::stop) to ensure
 *    clean C++ stack unwinding before the catch{} handlers call Rcpp::stop()
 *  - C++ objects are created independently with \c unique_ptr (exception-safe)
 *  - the result is written to the same HDF5 file and its location returned as
 *    a \c List(filename, path)
 *
 * @note
 * CRITICAL: Rcpp::stop() must only be called inside catch{} handlers,
 * never in the try{} body. Rcpp::stop() throws Rcpp::exception which,
 * when propagated through END_RCPP, causes Rf_error() → longjmp that
 * bypasses C++ destructors (including Rcpp::RNGScope) → SIGABRT.
 * Use throw std::runtime_error() for all error conditions inside try{}.
 *
 * This file is part of the HDF5Matrix R6+S3 architecture introduced for the
 * JSS revision. The underlying C++ algorithms are unchanged.
 */

/**
 * IMPORTANT: Uses Rf_error() instead of Rcpp::stop() in catch blocks.
 * Reason: Rcpp::stop() with Rcpp 1.1.0 on macOS ARM64 causes memory 
 * corruption (double-free) when exceptions propagate through END_RCPP.
 * Rf_error() bypasses Rcpp exception handling → no corruption.
 */


#include <BigDataStatMeth.hpp>
#include <sstream>
#include <ctime>
#include <cstdlib>

// ---------------------------------------------------------------------------
// Internal helper: generate a unique temporary dataset name
// ---------------------------------------------------------------------------
// static std::string make_temp_name(const char* prefix)
// {
//     std::ostringstream ss;
//     ss << prefix << "_" << std::time(nullptr) << "_" << std::rand();
//     return ss.str();
// }


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_multiply
// ---------------------------------------------------------------------------

//' General matrix product for HDF5 datasets (R6 wrapper)
//'
//' @description
//' Computes \code{A \%*\% B} (or transposed variants) for two HDF5 datasets
//' referenced by external pointers, using the BigDataStatMeth block-wise
//' multiplication algorithm.
//'
//' @param ptr_a External pointer (SEXP) for matrix A
//' @param ptr_b External pointer (SEXP) for matrix B
//' @param transpose_a Logical; transpose A before multiplying
//' @param transpose_b Logical; transpose B before multiplying
//' @param paral Logical or NULL; enable OpenMP parallelisation
//' @param block_size Integer or NULL; block size (NULL = auto)
//' @param threads Integer or NULL; thread count when \code{paral = TRUE}
//' @param outgroup   Character or NULL. Output group in the HDF5 file.
//'   Default \code{"OUTPUT"}.
//' @param outdataset Character or NULL. Output dataset name.
//'   Default \code{"A_x_B"} where A and B are the input dataset names.
//'
//' @return Named list with \code{filename} (character) and \code{path}
//'   (character) locating the result dataset within the HDF5 file.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_multiply(SEXP ptr_a,
                                     SEXP ptr_b,
                                     bool transpose_a = false,
                                     bool transpose_b = false,
                                     Rcpp::Nullable<bool>        paral       = R_NilValue,
                                     Rcpp::Nullable<int>         block_size  = R_NilValue,
                                     Rcpp::Nullable<int>         threads     = R_NilValue,
                                     Rcpp::Nullable<int>         compression = R_NilValue,
                                     Rcpp::Nullable<std::string> outgroup    = R_NilValue,
                                     Rcpp::Nullable<std::string> outdataset  = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "", Rcpp::Named("path") = "");

    try {
        H5::Exception::dontPrint();

        auto* dsA_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(ptr_a));
        auto* dsB_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(ptr_b));

        if (dsA_raw == nullptr || dsB_raw == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!dsA_raw->isOpen() || !dsB_raw->isOpen())
            throw std::runtime_error("Dataset is closed");

        //.. 20260304 ..// const std::string filename = dsA_raw->getFileName();
        const std::string filename = dsA_raw->getFullPath();
        const std::string groupA   = dsA_raw->getGroup();
        const std::string nameA    = dsA_raw->getDatasetName();
        const std::string groupB   = dsB_raw->getGroup();
        const std::string nameB    = dsB_raw->getDatasetName();

        hsize_t inner_a = transpose_a ? dsA_raw->nrows_r() : dsA_raw->ncols_r();
        hsize_t inner_b = transpose_b ? dsB_raw->ncols_r() : dsB_raw->nrows_r();
        if (inner_a != inner_b)
            throw std::runtime_error("Non-conformable matrices");

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        const std::string out_group = outgroup.isNull()
            ? std::string("OUTPUT")
                : Rcpp::as<std::string>(outgroup);
        const std::string result_name = outdataset.isNull()
            ? (nameA + "_x_" + nameB)
            : Rcpp::as<std::string>(outdataset);
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, out_group, result_name, true));
        dsC->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());

        BigDataStatMeth::multiplication(dsA.get(), dsB.get(), dsC.get(),
                                        transpose_a, transpose_b,
                                        paral, block_size, threads);

        lst["filename"] = filename;
        lst["path"] = out_group + "/" + result_name;

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error: %s", e.what());
    }

    return lst;
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_crossprod
// ---------------------------------------------------------------------------

//' Cross product for HDF5 datasets (R6 wrapper)
//'
//' @description
//' Computes \code{t(A) \%*\% B} using the dedicated BigDataStatMeth
//' block-wise cross-product algorithm. When A and B refer to the same
//' dataset, the symmetric optimisation (\code{bisSymetric = TRUE}) is
//' applied automatically.
//'
//' @param ptr_a External pointer (SEXP) for matrix A
//' @param ptr_b External pointer (SEXP) for matrix B
//' @param paral Logical or NULL; enable OpenMP parallelisation
//' @param block_size Integer or NULL; block size (NULL = auto)
//' @param threads Integer or NULL; thread count when \code{paral = TRUE}
//' @param outgroup   Character or NULL. Output group in the HDF5 file.
//'   Default \code{"OUTPUT"}.
//' @param outdataset Character or NULL. Output dataset name.
//'   Default \code{"CrossProd_A"} (single matrix) or
//'   \code{"CrossProd_A_x_B"} (two matrices).
//'
//' @return Named list with \code{filename} and \code{path} of the result.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_crossprod(SEXP ptr_a,
                                      SEXP ptr_b,
                                      Rcpp::Nullable<bool> paral       = R_NilValue,
                                      Rcpp::Nullable<int>  block_size  = R_NilValue,
                                      Rcpp::Nullable<int>  threads     = R_NilValue,
                                      Rcpp::Nullable<int>  compression = R_NilValue,
                                      Rcpp::Nullable<std::string> outgroup    = R_NilValue,
                                      Rcpp::Nullable<std::string> outdataset  = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(Rcpp::Named("filename") = "",
                                        Rcpp::Named("path")     = "");
    try {
        H5::Exception::dontPrint();

        auto* dsA_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(ptr_a));
        auto* dsB_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(ptr_b));

        if (dsA_raw == nullptr || dsB_raw == nullptr)
            throw std::runtime_error("Invalid external pointer");

        //.. 20260304 ..// const std::string filename = dsA_raw->getFileName();
        const std::string filename = dsA_raw->getFullPath();
        const std::string groupA   = dsA_raw->getGroup();
        const std::string nameA    = dsA_raw->getDatasetName();
        const std::string groupB   = dsB_raw->getGroup();
        const std::string nameB    = dsB_raw->getDatasetName();

        const bool bisSymetric = (groupA == groupB && nameA == nameB);

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();
        
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        if (dsA->nrows_r() != dsB->nrows_r())
            throw std::runtime_error("Non-conformable matrices for crossprod");

        const int iblockfactor = 2;
        const int iblock_size  = BigDataStatMeth::getMaxBlockSize(
            dsA->nrows(), dsA->ncols(), dsB->nrows(), dsB->ncols(),
            iblockfactor, block_size);

        const std::string out_group = outgroup.isNull()
            ? std::string("OUTPUT")
                : Rcpp::as<std::string>(outgroup);
        const std::string result_name = outdataset.isNull()
            ? (bisSymetric ? ("CrossProd_" + nameA)
                   : ("CrossProd_" + nameA + "_x_" + nameB))
            : Rcpp::as<std::string>(outdataset);
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, out_group, result_name, true));
        dsC->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());

        bool bparal = false;
        if (paral.isNotNull()) bparal = Rcpp::as<bool>(paral);

        if (bparal) {
            const int memory_block = iblock_size / 2;
            BigDataStatMeth::crossprod(dsA.get(), dsB.get(), dsC.get(),
                                       bisSymetric, iblock_size, memory_block,
                                       true, true, threads);
        } else {
            BigDataStatMeth::crossprod(dsA.get(), dsB.get(), dsC.get(),
                                       bisSymetric, iblock_size, 0,
                                       false, true, threads);
        }

        lst["filename"] = filename;
        lst["path"] = out_group + "/" + result_name;

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error: %s", e.what());
    }

    return lst;
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_tcrossprod
// ---------------------------------------------------------------------------

//' Transposed cross product for HDF5 datasets (R6 wrapper)
//'
//' @description
//' Computes \code{A \%*\% t(B)} using the dedicated BigDataStatMeth
//' block-wise transposed cross-product algorithm. When A and B refer to the
//' same dataset, the symmetric optimisation is applied automatically.
//'
//' @param ptr_a External pointer (SEXP) for matrix A
//' @param ptr_b External pointer (SEXP) for matrix B
//' @param paral Logical or NULL; enable OpenMP parallelisation
//' @param block_size Integer or NULL; block size (NULL = auto)
//' @param threads Integer or NULL; thread count when \code{paral = TRUE}
//' @param outgroup   Character or NULL. Output group in the HDF5 file.
//'   Default \code{"OUTPUT"}.
//' @param outdataset Character or NULL. Output dataset name.
//'   Default \code{"tCrossProd_A"} (single matrix) or
//'   \code{"tCrossProd_A_x_B"} (two matrices).
//'
//' @return Named list with \code{filename} and \code{path} of the result.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_tcrossprod(SEXP ptr_a,
                                       SEXP ptr_b,
                                       Rcpp::Nullable<bool> paral       = R_NilValue,
                                       Rcpp::Nullable<int>  block_size  = R_NilValue,
                                       Rcpp::Nullable<int>  threads     = R_NilValue,
                                       Rcpp::Nullable<int>  compression = R_NilValue,
                                       Rcpp::Nullable<std::string> outgroup   = R_NilValue,
                                       Rcpp::Nullable<std::string> outdataset = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(Rcpp::Named("filename") = "",
                                        Rcpp::Named("path")     = "");
    try {
        H5::Exception::dontPrint();

        auto* dsA_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(ptr_a));
        auto* dsB_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(ptr_b));

        if (dsA_raw == nullptr || dsB_raw == nullptr)
            throw std::runtime_error("Invalid external pointer");

        //.. 20260304 ..// const std::string filename = dsA_raw->getFileName();
        const std::string filename = dsA_raw->getFullPath();
        const std::string groupA   = dsA_raw->getGroup();
        const std::string nameA    = dsA_raw->getDatasetName();
        const std::string groupB   = dsB_raw->getGroup();
        const std::string nameB    = dsB_raw->getDatasetName();

        const bool bisSymetric = (groupA == groupB && nameA == nameB);

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();
        
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        if (dsA->ncols_r() != dsB->ncols_r())
            throw std::runtime_error("Non-conformable matrices for tcrossprod");

        const int iblockfactor = 2;
        const int iblock_size  = BigDataStatMeth::getMaxBlockSize(
            dsA->nrows(), dsA->ncols(), dsB->nrows(), dsB->ncols(),
            iblockfactor, block_size);

        const std::string out_group = outgroup.isNull()
            ? std::string("OUTPUT")
                : Rcpp::as<std::string>(outgroup);
        const std::string result_name = outdataset.isNull()
            ? (bisSymetric ? ("tCrossProd_" + nameA)
                   : ("tCrossProd_" + nameA + "_x_" + nameB))
            : Rcpp::as<std::string>(outdataset);
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, out_group, result_name, true));
        dsC->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());

        bool bparal = false;
        if (paral.isNotNull()) bparal = Rcpp::as<bool>(paral);

        if (bparal) {
            const int memory_block = iblock_size / 2;
            BigDataStatMeth::tcrossprod(dsA.get(), dsB.get(), dsC.get(),
                                        bisSymetric, iblock_size, memory_block,
                                        true, true, threads);
        } else {
            BigDataStatMeth::tcrossprod(dsA.get(), dsB.get(), dsC.get(),
                                        bisSymetric, iblock_size, 0,
                                        false, true, threads);
        }

        lst["filename"] = filename;
        lst["path"] = out_group + "/" + result_name;

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error: %s", e.what());
    }

    return lst;
}
