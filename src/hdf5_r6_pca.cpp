/**
 * @file hdf5_r6_pca.cpp
 * @brief Rcpp wrapper for block-wise PCA of HDF5-stored matrices.
 *
 * Receives all inputs as plain strings / scalars (never SEXP pointers) and
 * delegates entirely to BigDataStatMeth::RcppPCAHdf5() defined in
 * inst/include/hdf5Algebra/matrixPCA.hpp.  All HDF5 handles are managed
 * internally by the header function with unique_ptr — no handle conflicts
 * can arise regardless of how many HDF5Matrix objects the caller has open.
 *
 * Output datasets are always written to the input file under the group
 * "PCA/<dataset>/":
 *   - lambda       : eigenvalues
 *   - variance     : variance explained by each PC
 *   - cumvar       : cumulative variance explained
 *   - var.coord    : variable coordinates (loadings / rotation matrix)
 *   - var.cos2     : squared cosines for variables
 *   - ind.dist     : distances of individuals from origin
 *   - components   : principal components (projected data)
 *   - ind.coord    : individual coordinates on PCs
 *   - ind.cos2     : squared cosines for individuals
 *   - ind.contrib  : contributions of individuals to each PC
 */

#include <BigDataStatMeth.hpp>

/**
 * @brief Compute Principal Component Analysis of an HDF5 matrix.
 *
 * Calls BigDataStatMeth::RcppPCAHdf5() (matrixPCA.hpp), which runs
 * block-wise SVD internally and derives all PCA statistics on disk.
 * If SVD results already exist in the file (group \p svdgroup) and
 * \p overwrite is false, the SVD step is skipped.
 *
 * @param filename     Path to the HDF5 file containing the input matrix.
 * @param group        Group path of the input dataset.
 * @param dataset      Dataset name within the group.
 * @param ncomponents  Number of PCs to compute (0 = all).
 * @param bcenter      If true, center columns before PCA.
 * @param bscale       If true, scale columns before PCA.
 * @param k            Number of local SVDs per level (default 2).
 * @param q            Levels for incremental SVD (default 1).
 * @param rankthreshold Threshold in [0, 0.1] for rank approximation.
 * @param svdgroup     HDF5 group where intermediate SVD is stored/reused.
 *                     Empty string defaults to "SVD/".
 * @param overwrite    If true, recompute even if results exist.
 * @param method       Computation method: "auto" | "blocks" | "full".
 * @param threads      Number of OpenMP threads (-1 = auto).
 * @return Named list with keys:
 *   \c file, \c path_lambda, \c path_variance, \c path_cumvar,
 *   \c path_var_coord, \c path_var_cos2, \c path_ind_dist,
 *   \c path_components, \c path_ind_coord, \c path_ind_cos2,
 *   \c path_ind_contrib.
 *   All paths are full HDF5 paths inside \p filename.
 *
 * @export
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_pca(std::string filename,
                                 std::string group,
                                 std::string dataset,
                                 int         ncomponents   = 0,
                                 bool        bcenter       = false,
                                 bool        bscale        = false,
                                 int         k             = 2,
                                 int         q             = 1,
                                 double      rankthreshold = 0.0,
                                 std::string svdgroup      = "SVD/",
                                 bool        overwrite     = false,
                                 std::string method        = "auto",
                                 int         threads       = -1)
{
    try {
        H5::Exception::dontPrint();

        if (rankthreshold < 0.0 || rankthreshold > 0.1)
            Rf_error("rcpp_hdf5dataset_pca: rankthreshold must be in [0, 0.1]");

        // Ensure svdgroup ends with "/"
        if (svdgroup.empty()) svdgroup = "SVD/";
        if (svdgroup.back() != '/') svdgroup += '/';
        
        // Build Nullable wrappers expected by the header
        Rcpp::Nullable<Rcpp::CharacterVector> n_method =
            method.empty() ? R_NilValue
                           : Rcpp::wrap(Rcpp::CharacterVector::create(method));
        
        Rcpp::Nullable<int> n_threads =
            (threads < 0) ? R_NilValue : Rcpp::wrap(threads);
        
        // Delegate — all HDF5 I/O happens inside the header
        BigDataStatMeth::RcppPCAHdf5(
            filename, group, dataset,
            svdgroup,
            k, q,
            ncomponents,
            bcenter, bscale,
            rankthreshold,
            overwrite,
            /*asRowMajor=*/false,
            n_method,
            n_threads
        );
        
        // Build output paths (mirrors hdf5_bdPCA.cpp convention)
        const std::string pca_root = "PCA/" + dataset + "/";
        
        return Rcpp::List::create(
            Rcpp::Named("file")           = filename,
            Rcpp::Named("path_lambda")    = pca_root + "lambda",
            Rcpp::Named("path_variance")  = pca_root + "variance",
            Rcpp::Named("path_cumvar")    = pca_root + "cumvar",
            Rcpp::Named("path_var_coord") = pca_root + "var.coord",
            Rcpp::Named("path_var_cos2")  = pca_root + "var.cos2",
            Rcpp::Named("path_ind_dist")  = pca_root + "ind.dist",
            Rcpp::Named("path_components")= pca_root + "components",
            Rcpp::Named("path_ind_coord") = pca_root + "ind.coord",
            Rcpp::Named("path_ind_cos2")  = pca_root + "ind.cos2",
            Rcpp::Named("path_ind_contrib")= pca_root + "ind.contrib"
        );

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_pca (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_pca (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSpaceIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_pca (DataSpace IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_pca: %s", e.what());
    }
    return R_NilValue;
}
