/**
 * @file hdf5_r6_cholesky.cpp
 * @brief Rcpp wrappers for Cholesky and inverse-Cholesky decomposition of
 *        HDF5-stored symmetric positive-definite matrices.
 *
 * ARCHITECTURAL RULE (roadmap §1.4):
 *   Functions marked [[Rcpp::export]] (e.g. bdCholesky_hdf5, bdInvCholesky_hdf5)
 *   are compiled as DLL symbols accessible ONLY from R — never from another .cpp.
 *   All wrappers must call the inline header functions directly:
 *     - Cholesky:     BigDataStatMeth::Cholesky_decomposition_hdf5()
 *                     (inst/include/hdf5Algebra/matrixInvCholesky.hpp)
 *     - Inverse:      BigDataStatMeth::Rcpp_InvCholesky_hdf5()
 *                     (inst/include/hdf5Algebra/matrixInvCholesky.hpp)
 *
 * Output layout in the HDF5 file (defaults, same group as input):
 *   rcpp_hdf5dataset_chol  → <group>/chol_<dataset>   (lower-triangular L)
 *   rcpp_hdf5dataset_solve → <group>/inv_<dataset>    (full inverse A⁻¹)
 */

#include <BigDataStatMeth.hpp>

/**
 * @brief Compute the Cholesky decomposition of a symmetric positive-definite
 *        HDF5 matrix.
 *
 * Computes the lower-triangular factor L such that A = L L'.
 * Delegates to \c BigDataStatMeth::Cholesky_decomposition_hdf5() from
 * \c inst/include/hdf5Algebra/matrixInvCholesky.hpp.
 *
 * @param filename    Path to the HDF5 file.
 * @param group       Group path of the input dataset.
 * @param dataset     Dataset name within the group.
 * @param full_matrix If true, return the full symmetric matrix (L + L').
 *                    Default false (lower-triangular L only).
 * @param overwrite   If true, overwrite existing result.  Default false.
 * @param threads     Number of OpenMP threads (-1 = auto).
 * @param block_size  Elements per block for block-wise computation.
 *                    -1 = auto (MAXELEMSINBLOCK).
 * @return Named list:
 *   \c file   — full path to the HDF5 file,
 *   \c path_L — full HDF5 path to the Cholesky factor.
 *
 * @export
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_chol(std::string filename,
                                  std::string group,
                                  std::string dataset,
                                  bool        full_matrix = false,
                                  bool        overwrite   = false,
                                  int         threads     = -1,
                                  long        block_size  = -1,
                                  Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        const std::string out_dataset = "chol_" + dataset;
        const std::string out_path    = group + "/" + out_dataset;

        long dElementsBlock = (block_size < 0) ?
            static_cast<long>(MAXELEMSINBLOCK) : block_size;

        Rcpp::Nullable<int> n_threads =
            (threads < 0) ? R_NilValue : Rcpp::wrap(threads);

        // Open input dataset
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false));
        dsA->openDataset();

        if (dsA->getDatasetptr() == nullptr)
            Rf_error("rcpp_hdf5dataset_chol: cannot open dataset '%s/%s'",
                     group.c_str(), dataset.c_str());

        int nrows = dsA->nrows(), ncols = dsA->ncols();
        if (nrows != ncols)
            Rf_error("rcpp_hdf5dataset_chol: matrix must be square (%d x %d)",
                     nrows, ncols);

        // Create output dataset
        std::unique_ptr<BigDataStatMeth::hdf5DatasetInternal> dstmp(
            new BigDataStatMeth::hdf5DatasetInternal(filename, out_path, overwrite));
        dstmp->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());
        dstmp->createDataset(nrows, ncols, "real");

        if (dstmp->getDatasetptr() == nullptr)
            Rf_error("rcpp_hdf5dataset_chol: cannot create output dataset '%s'",
                     out_path.c_str());

        // Call header function directly (not bdCholesky_hdf5 — Rule §1.4)
        //.. 20260302 ..// int res = BigDataStatMeth::Cholesky_decomposition_hdf5(
        //.. 20260302 ..//     dsA.get(), dstmp.get(), nrows, ncols, dElementsBlock, n_threads);
        int res = BigDataStatMeth::Cholesky_decomposition_hdf5(
            dsA.get(), dstmp.get(), nrows, ncols, dElementsBlock, full_matrix, n_threads);
        
        if (res != 0)
            Rf_error("rcpp_hdf5dataset_chol: Cholesky decomposition failed "
                     "(matrix may not be positive definite)");

        return Rcpp::List::create(
            Rcpp::Named("file")   = dsA->getFullPath(),
            Rcpp::Named("path_L") = out_path
        );

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_chol (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_chol (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_chol: %s", e.what());
    }
    return R_NilValue;
}


/**
 * @brief Compute the inverse of a symmetric positive-definite HDF5 matrix
 *        via Cholesky decomposition.
 *
 * Equivalent to base R \code{solve(A)} for SPD matrices.
 * Delegates to \c BigDataStatMeth::Rcpp_InvCholesky_hdf5() from
 * \c inst/include/hdf5Algebra/matrixInvCholesky.hpp.
 *
 * @param filename    Path to the HDF5 file.
 * @param group       Group path of the input dataset.
 * @param dataset     Dataset name within the group.
 * @param full_matrix If true, return full symmetric inverse.  Default true.
 * @param overwrite   If true, overwrite existing result.  Default false.
 * @param threads     Number of OpenMP threads (-1 = auto).
 * @param block_size  Elements per block.  -1 = auto (MAXELEMSINBLOCK).
 * @return Named list:
 *   \c file     — full path to the HDF5 file,
 *   \c path_inv — full HDF5 path to the inverse matrix.
 *
 * @export
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_solve(std::string filename,
                                   std::string group,
                                   std::string dataset,
                                   bool        full_matrix = true,
                                   bool        overwrite   = false,
                                   int         threads     = -1,
                                   long        block_size  = -1,
                                   Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        const std::string out_dataset = "inv_" + dataset;
        const std::string out_path    = group + "/" + out_dataset;

        long dElementsBlock = (block_size < 0) ?
            static_cast<long>(MAXELEMSINBLOCK) : block_size;

        Rcpp::Nullable<int> n_threads =
            (threads < 0) ? R_NilValue : Rcpp::wrap(threads);

        // Open input dataset
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false));
        dsA->openDataset();

        if (dsA->getDatasetptr() == nullptr)
            Rf_error("rcpp_hdf5dataset_solve: cannot open dataset '%s/%s'",
                     group.c_str(), dataset.c_str());

        int nrows = dsA->nrows(), ncols = dsA->ncols();
        if (nrows != ncols)
            Rf_error("rcpp_hdf5dataset_solve: matrix must be square (%d x %d)",
                     nrows, ncols);

        // Create output dataset
        std::unique_ptr<BigDataStatMeth::hdf5DatasetInternal> dstmp(
            new BigDataStatMeth::hdf5DatasetInternal(filename, out_path, overwrite));
        dstmp->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());
        dstmp->createDataset(nrows, ncols, "real");

        if (dstmp->getDatasetptr() == nullptr)
            Rf_error("rcpp_hdf5dataset_solve: cannot create output dataset '%s'",
                     out_path.c_str());

        // Call header function directly (not bdInvCholesky_hdf5 — Rule §1.4)
        BigDataStatMeth::Rcpp_InvCholesky_hdf5(
            dsA.get(), dstmp.get(), full_matrix, dElementsBlock, n_threads);

        return Rcpp::List::create(
            Rcpp::Named("file")    = dsA->getFullPath(),
            Rcpp::Named("path_inv")= out_path
        );

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_solve (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_solve (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_solve: %s", e.what());
    }
    return R_NilValue;
}
