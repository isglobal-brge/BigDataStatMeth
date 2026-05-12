/**
 * @file hdf5_r6_svd.cpp
 * @brief Rcpp wrapper for block-wise SVD of HDF5-stored matrices.
 *
 * Receives all inputs as plain strings / scalars (never SEXP pointers) and
 * delegates entirely to BigDataStatMeth::RcppbdSVD_hdf5() defined in
 * inst/include/hdf5Algebra/matrixSvd.hpp.  All HDF5 handles are managed
 * internally by the header function with unique_ptr — no handle conflicts
 * can arise regardless of how many HDF5Matrix objects the caller has open.
 *
 * Output datasets are always written to the input file under the group
 * "SVD/<dataset>/":
 *   - d : 1 × min(m,n) vector of singular values (non-negative, decreasing)
 *   - u : m × r matrix of left  singular vectors
 *   - v : n × r matrix of right singular vectors
 */

#include <BigDataStatMeth.hpp>

/**
 * @brief Compute the Singular Value Decomposition of an HDF5 matrix.
 *
 * Calls BigDataStatMeth::RcppbdSVD_hdf5() (matrixSvd.hpp) which selects
 * automatically between a direct LAPACK path (small matrices) and a
 * block-wise incremental algorithm (large matrices).
 *
 * @param filename   Path to the HDF5 file containing the input matrix.
 * @param group      Group path of the input dataset.
 * @param dataset    Dataset name within the group.
 * @param k          Number of local SVDs to concatenate per level (default 2).
 * @param q          Number of levels for incremental computation (default 1).
 * @param nev        Number of singular values/vectors to compute (0 = all).
 * @param bcenter    If true, center columns before decomposition.
 * @param bscale     If true, scale columns before decomposition.
 * @param rankthreshold Threshold in [0, 0.1] for rank approximation.
 * @param overwrite  If true, overwrite existing SVD results.
 * @param method     Computation method: "auto" | "blocks" | "full".
 *                   Empty string is treated as "auto".
 * @param threads    Number of OpenMP threads (-1 = auto).
 * @return Named list:
 *   \c file    — path to the HDF5 file (same as input),
 *   \c path_d  — full HDF5 path to singular values  (SVD/<dataset>/d),
 *   \c path_u  — full HDF5 path to left  vectors    (SVD/<dataset>/u),
 *   \c path_v  — full HDF5 path to right vectors    (SVD/<dataset>/v).
 *
 * @export
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_svd(std::string filename,
                                 std::string group,
                                 std::string dataset,
                                 int         k              = 2,
                                 int         q              = 1,
                                 int         nev            = 0,
                                 bool        bcenter        = true,
                                 bool        bscale         = true,
                                 double      rankthreshold  = 0.0,
                                 bool        overwrite      = false,
                                 std::string method         = "auto",
                                 int         threads        = -1)
{
    try {
        H5::Exception::dontPrint();

        // Validate threshold
        if (rankthreshold < 0.0 || rankthreshold > 0.1)
            Rf_error("rcpp_hdf5dataset_svd: rankthreshold must be in [0, 0.1]");

        // Build Nullable wrappers expected by the header function
        Rcpp::Nullable<Rcpp::CharacterVector> n_method =
            method.empty() ? R_NilValue
                           : Rcpp::wrap(Rcpp::CharacterVector::create(method));

        Rcpp::Nullable<int> n_threads =
            (threads < 0) ? R_NilValue : Rcpp::wrap(threads);

        // Delegate — all HDF5 I/O happens inside the header
        BigDataStatMeth::RcppbdSVD_hdf5(
            filename, group, dataset,
            k, q, nev,
            bcenter, bscale,
            rankthreshold,
            overwrite,
            /*asRowMajor=*/false,
            n_method,
            n_threads
        );

        // Build output paths (mirrors hdf5_SVD.cpp convention)
        const std::string svd_root = "SVD/" + dataset + "/";

        return Rcpp::List::create(
            Rcpp::Named("file")   = filename,
            Rcpp::Named("path_d") = svd_root + "d",
            Rcpp::Named("path_u") = svd_root + "u",
            Rcpp::Named("path_v") = svd_root + "v"
        );

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_svd (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_svd (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_svd: %s", e.what());
    }
    return R_NilValue;
}
