/**
 * @file hdf5_r6_qr.cpp
 * @brief Rcpp wrapper for block-wise QR decomposition of HDF5-stored matrices.
 *
 * ARCHITECTURAL RULE (roadmap §1.4):
 *   Functions marked [[Rcpp::export]] (e.g. bdQR_hdf5) are compiled as DLL
 *   symbols accessible ONLY from R — never from another .cpp.
 *   This wrapper calls the inline header functions directly:
 *     BigDataStatMeth::RcppQRHdf5_dispatch()   (new unified entry point)
 *     BigDataStatMeth::RcppQRHdf5()            (Householder / LAPACK)
 *     BigDataStatMeth::RcppTSQRHdf5()          (Tall-Skinny QR, parallel)
 *   All defined in inst/include/hdf5Algebra/matrixQR.hpp
 *
 * Output layout in the HDF5 file (default, relative to input group):
 *   <group>/QRDec/Q.<dataset>  — Q factor
 *   <group>/QRDec/R.<dataset>  — R factor
 *
 * Q and R dimensions depend on method and thin flag:
 *   method="lapack", thin=FALSE : Q m×m,   R m×n
 *   method="lapack", thin=TRUE  : Q m×k,   R k×n   (k = min(m,n))
 *   method="tsqr",   thin=FALSE : Q m×m,   R m×n   (basis completion, expensive)
 *   method="tsqr",   thin=TRUE  : Q m×n,   R n×n   (natural TSQR output)
 *   method="auto"               : selects tsqr when m > 5n AND m > 1000,
 *                                 lapack otherwise
 *
 * @note 20260304: Added method parameter and TSQR support.
 */

#include <BigDataStatMeth.hpp>

/**
 * @brief Compute the QR decomposition of an HDF5 matrix.
 *
 * Dispatches to BigDataStatMeth::RcppQRHdf5_dispatch() which selects
 * the Householder (LAPACK) or TSQR algorithm based on the \c method argument.
 *
 * @param filename   Path to the HDF5 file containing the input matrix.
 * @param group      Group path of the input dataset.
 * @param dataset    Dataset name within the group.
 * @param thin       If true, compute thin (economy) QR.  Default false.
 *                   For method="tsqr" this is strongly recommended (full Q
 *                   via TSQR requires an expensive basis-completion step).
 * @param block_size Block size hint for TSQR row-block partitioning. -1 = auto.
 *                   Ignored by method="lapack".
 * @param overwrite  If true, overwrite existing QR results.  Default false.
 * @param threads    Number of OpenMP threads (-1 = auto, CRAN-compliant).
 *                   For method="lapack" this controls Eigen::setNbThreads.
 *                   For method="tsqr" this controls the OpenMP block loop.
 * @param method     Algorithm selection:
 *                   \c "auto"   (default) — TSQR for tall-skinny, LAPACK otherwise.
 *                   \c "lapack" — Eigen HouseholderQR, single-threaded Householder.
 *                   \c "tsqr"   — Parallel TSQR; requires m >= n.
 * @return Named list:
 *   \c file   — full path to the HDF5 file,
 *   \c path_Q — full HDF5 path to Q factor,
 *   \c path_R — full HDF5 path to R factor.
 *
 * @export
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_qr(std::string filename,
                                std::string group,
                                std::string dataset,
                                bool        thin       = false,
                                int         block_size = -1,
                                bool        overwrite  = false,
                                int         threads    = -1,
                                std::string method     = "auto",
                                Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        // Output paths — mirror bdQR_hdf5 convention
        const std::string out_group   = group + "/QRDec";
        const std::string ds_Q        = "Q." + dataset;
        const std::string ds_R        = "R." + dataset;
        const std::string path_Q      = out_group + "/" + ds_Q;
        const std::string path_R      = out_group + "/" + ds_R;

        Rcpp::Nullable<int> n_block   =
            (block_size < 0) ? R_NilValue : Rcpp::wrap(block_size);
        Rcpp::Nullable<int> n_threads =
            (threads    < 0) ? R_NilValue : Rcpp::wrap(threads);

        // Open input dataset
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false));
        dsA->openDataset();

        if (dsA->getDatasetptr() == nullptr)
            Rf_error("rcpp_hdf5dataset_qr: cannot open dataset '%s/%s'",
                     group.c_str(), dataset.c_str());

        // Create output datasets for Q and R
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsQ(
            new BigDataStatMeth::hdf5Dataset(filename, out_group, ds_Q, overwrite));
        dsQ->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsR(
            new BigDataStatMeth::hdf5Dataset(filename, out_group, ds_R, overwrite));
        dsR->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());

        // Dispatch to LAPACK or TSQR (Rule §1.4: call header function, not bdQR_hdf5)
        BigDataStatMeth::RcppQRHdf5_dispatch(
            dsA.get(), dsQ.get(), dsR.get(),
            thin, n_block, n_threads, method);

        return Rcpp::List::create(
            Rcpp::Named("file")   = dsA->getFullPath(),
            Rcpp::Named("path_Q") = path_Q,
            Rcpp::Named("path_R") = path_R
        );

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_qr (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_qr (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_qr: %s", e.what());
    }
    return R_NilValue;
}
