/**
 * @file hdf5_r6_pseudoinv.cpp
 * @brief R6 wrapper for Moore-Penrose pseudoinverse of HDF5-stored matrices.
 *
 * ARCHITECTURAL RULE (roadmap §1.4):
 *   Calls BigDataStatMeth::RcppPseudoinvHdf5() defined in
 *   inst/include/hdf5Algebra/matrixPseudoinverse.hpp directly.
 *   Never calls the Rcpp-exported bdpseudoinv_hdf5() symbol.
 *
 * Output layout:
 *   PseudoInverse/<dataset>  -- m x n pseudoinverse of an n x m input
 *
 * @note Uses Rf_error() in catch blocks to avoid macOS ARM64 memory corruption.
 */

#include <BigDataStatMeth.hpp>

/**
 * @brief Moore-Penrose pseudoinverse of an HDF5 matrix (SVD-based).
 *
 * For A (n x m), computes A+ (m x n) via thin SVD: A+ = V Σ+ Uᵀ.
 * Delegates to BigDataStatMeth::RcppPseudoinvHdf5() (matrixPseudoinverse.hpp).
 *
 * @param filename    Path to the HDF5 file.
 * @param group       Group path of the input dataset.
 * @param dataset     Dataset name within the group.
 * @param out_group   Output group (default "PseudoInverse").
 * @param out_dataset Output dataset name (default = input dataset name).
 * @param overwrite   Overwrite existing result (default false).
 * @param threads     OpenMP threads (-1 = auto).
 * @param compression gzip level 0-9 (-1 = inherit from input).
 * @return Named list:
 *   \c file -- full HDF5 file path,
 *   \c path -- full HDF5 path to the pseudoinverse dataset.
 *
 * @export
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_pseudoinv(std::string filename,
                                       std::string group,
                                       std::string dataset,
                                       std::string out_group   = "PseudoInverse",
                                       std::string out_dataset = "",
                                       bool        overwrite   = false,
                                       int         threads     = -1,
                                       int         compression = -1)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("file") = "",
        Rcpp::Named("path") = "");

    try {
        H5::Exception::dontPrint();

        if (out_dataset.empty()) out_dataset = dataset;

        Rcpp::Nullable<int> n_threads =
            (threads < 0) ? R_NilValue : Rcpp::wrap(threads);

        // Open input dataset
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false));
        dsA->openDataset();

        if (dsA->getDatasetptr() == nullptr)
            throw std::runtime_error(
                "rcpp_hdf5dataset_pseudoinv: cannot open '" +
                group + "/" + dataset + "'");

        int comp_level = (compression < 0) ? dsA->getCompressionLevel()
                                           : compression;

        // Create output dataset (transposed dimensions: m x n for n x m input)
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsR(
            new BigDataStatMeth::hdf5Dataset(
                filename, out_group, out_dataset, overwrite));
        dsR->setCompressionLevel(comp_level);

        // Delegate to header (Rule §1.4)
        BigDataStatMeth::RcppPseudoinvHdf5(dsA.get(), dsR.get(), n_threads);

        lst["file"] = dsA->getFullPath();
        lst["path"] = out_group + "/" + out_dataset;

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_pseudoinv (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_pseudoinv (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_pseudoinv: %s", e.what());
    }

    return lst;
}
