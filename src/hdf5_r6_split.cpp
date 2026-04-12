/**
 * @file hdf5_r6_split.cpp
 * @brief R6 wrapper for splitting an HDF5 dataset into blocks.
 *
 * ARCHITECTURAL RULE (roadmap section 1.4):
 *   Calls BigDataStatMeth::RcppSplit_matrix_hdf5() defined in
 *   inst/include/hdf5Utilities/hdf5SplitDataset.hpp directly.
 *   Never calls the Rcpp-exported bdSplit_matrix_hdf5() symbol.
 *
 * Output layout (produced by the header function):
 *   out_group/out_dataset.0,  out_group/out_dataset.1, ...
 *
 * @note Uses Rf_error() in catch blocks to avoid macOS ARM64 double-free.
 */

#include <BigDataStatMeth.hpp>


/**
 * @brief Split an HDF5 dataset into blocks (R6 wrapper).
 *
 * Splits the HDF5 dataset referenced by ptr into n_blocks equal (or
 * approximately equal) blocks along rows (bycols = FALSE) or columns
 * (bycols = TRUE). Exactly one of n_blocks or block_size must be provided
 * (not both, and not neither).
 * Delegates to BigDataStatMeth::RcppSplit_matrix_hdf5().
 *
 * Output block i (0-based) is stored at out_group/out_dataset.i in the
 * same HDF5 file.
 *
 * @param ptr        External pointer (SEXP) to the input dataset.
 * @param bycols     If TRUE, split by columns; if FALSE (default), split by rows.
 * @param n_blocks   Number of (approximately equal) blocks. Ignored when -1.
 * @param block_size Exact block size (rows or cols per block). Ignored when -1.
 * @param out_group  Output group path (default "SPLIT").
 * @param out_dataset Base dataset name; actual names will be base.0, base.1, etc.
 *                    Default = input dataset name.
 * @param overwrite  Overwrite existing blocks (passed through). Default false.
 * @return Named list with elements "filename", "out_group", "out_dataset", "n_blocks".
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_split(SEXP        ptr,
                                   bool        bycols      = false,
                                   int         n_blocks    = -1,
                                   int         block_size  = -1,
                                   std::string out_group   = "SPLIT",
                                   std::string out_dataset = "",
                                   bool        overwrite   = false)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename")    = "",
        Rcpp::Named("out_group")   = out_group,
        Rcpp::Named("out_dataset") = "",
        Rcpp::Named("n_blocks")    = 0);

    try {
        H5::Exception::dontPrint();

        auto* raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr));
        if (raw == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!raw->isOpen())
            throw std::runtime_error("Dataset is closed");

        const std::string filename = raw->getFullPath();
        const std::string group    = raw->getGroup();
        const std::string name     = raw->getDatasetName();

        if (out_dataset.empty()) out_dataset = name;

        // Fresh handle
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> ds(
            new BigDataStatMeth::hdf5Dataset(filename, group, name, false));
        ds->openDataset();

        if (ds->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open dataset");

        int irows = static_cast<int>(ds->nrows_r());
        int icols = static_cast<int>(ds->ncols_r());
        int dim   = bycols ? icols : irows;   // dimension being split

        // Resolve block_size
        int iblock;
        if (n_blocks > 0 && block_size > 0)
            throw std::runtime_error(
                "split: provide exactly one of n_blocks or block_size, not both");
        if (n_blocks < 0 && block_size < 0)
            throw std::runtime_error(
                "split: provide exactly one of n_blocks or block_size");
        if (n_blocks > 0)
            iblock = (dim + n_blocks - 1) / n_blocks;  // ceiling division
        else
            iblock = block_size;

        if (iblock <= 0)
            throw std::runtime_error("split: block_size must be > 0");

        int actual_blocks = (dim + iblock - 1) / iblock;

        // Delegate (Rule 1.4)
        BigDataStatMeth::RcppSplit_matrix_hdf5(
            ds.get(), bycols,
            out_group, out_dataset,
            iblock, irows, icols);

        lst["filename"]    = filename;
        lst["out_group"]   = out_group;
        lst["out_dataset"] = out_dataset;
        lst["n_blocks"]    = actual_blocks;

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_split (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_split (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_split: %s", e.what());
    }

    return lst;
}
