/**
 * @file hdf5_r6_eigen.cpp
 * @brief R6 wrapper for eigenvalue decomposition of HDF5-stored matrices.
 *
 * ARCHITECTURAL RULE (roadmap §1.4):
 *   Calls BigDataStatMeth::RcppbdEigen_hdf5() defined in
 *   inst/include/hdf5Algebra/matrixEigenDecomposition.hpp directly.
 *   Never calls the Rcpp-exported bdEigen_hdf5() symbol.
 *
 * Output layout (written by the header function):
 *   EIGEN/<dataset>/values   -- 1 x k real eigenvalues (descending)
 *   EIGEN/<dataset>/vectors  -- n x k real eigenvectors (if compute_vectors)
 *
 * @note Uses Rf_error() in catch blocks (not Rcpp::stop()) to avoid
 *   double-free / memory corruption on macOS ARM64 with Rcpp 1.1.0.
 */

#include <BigDataStatMeth.hpp>

/**
 * @brief Eigenvalue decomposition of an HDF5 matrix.
 *
 * Delegates to BigDataStatMeth::RcppbdEigen_hdf5() (Spectra-based,
 * block-wise for large matrices).
 *
 * @param filename        Path to the HDF5 file.
 * @param group           Group path of the input dataset.
 * @param dataset         Dataset name within the group.
 * @param k               Number of eigenvalues/vectors to compute (0 = all).
 * @param which           Which eigenvalues: "LM" (largest magnitude, default),
 *                        "SM", "LR", "SR", "LI", "SI".
 * @param ncv             Number of Lanczos vectors (0 = auto).
 * @param bcenter         Center columns before decomposition.
 * @param bscale          Scale columns before decomposition.
 * @param tolerance       Convergence tolerance (default 1e-10).
 * @param max_iter        Maximum Lanczos iterations (default 1000).
 * @param compute_vectors Whether to compute eigenvectors (default true).
 * @param overwrite       Overwrite existing results (default false).
 * @param threads         OpenMP threads (-1 = auto).
 * @return Named list:
 *   \c file      -- path to the HDF5 file,
 *   \c path_values  -- full HDF5 path to eigenvalues,
 *   \c path_vectors -- full HDF5 path to eigenvectors ("" if not computed).
 *
 * @export
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_eigen(std::string filename,
                                   std::string group,
                                   std::string dataset,
                                   int         k               = 0,
                                   std::string which           = "LM",
                                   int         ncv             = 0,
                                   bool        bcenter         = false,
                                   bool        bscale          = false,
                                   double      tolerance       = 1e-10,
                                   int         max_iter        = 1000,
                                   bool        compute_vectors = true,
                                   bool        overwrite       = false,
                                   int         threads         = -1)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("file")         = "",
        Rcpp::Named("path_values")  = "",
        Rcpp::Named("path_vectors") = "");

    try {
        H5::Exception::dontPrint();

        // Output paths -- mirror bdEigen_hdf5 convention
        const std::string out_group   = "EIGEN/" + dataset;
        const std::string path_values = out_group + "/values";
        const std::string path_vectors = compute_vectors
                                          ? out_group + "/vectors"
                                          : "";

        Rcpp::Nullable<int> n_threads =
            (threads < 0) ? R_NilValue : Rcpp::wrap(threads);

        // Delegate -- all HDF5 I/O is managed inside the header function
        BigDataStatMeth::RcppbdEigen_hdf5(
            filename, group, dataset,
            k, which, ncv,
            bcenter, bscale,
            tolerance, max_iter,
            compute_vectors, overwrite,
            n_threads);

        // Build canonical filename (full path, consistent with other wrappers)
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> ds_tmp(
            new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false));
        ds_tmp->openDataset();
        const std::string full_path = ds_tmp->getFullPath();

        lst["file"]         = full_path;
        lst["path_values"]  = path_values;
        lst["path_vectors"] = path_vectors;

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_eigen (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_eigen (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_eigen: %s", e.what());
    }

    return lst;
}
