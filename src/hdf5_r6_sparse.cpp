/**
 * @file hdf5_r6_sparse.cpp
 * @brief R6 wrapper for sparse matrix multiplication of HDF5-stored matrices.
 *
 * ARCHITECTURAL RULE (roadmap section 1.4):
 *   Calls BigDataStatMeth::multiplicationSparse() defined in
 *   inst/include/hdf5Algebra/multiplicationSparse.hpp directly.
 *   Never calls the Rcpp-exported bdblockmult_sparse_hdf5() symbol.
 *
 * @note Uses Rf_error() in catch blocks to avoid macOS ARM64 double-free.
 */

#include <BigDataStatMeth.hpp>
#include <sstream>
#include <ctime>
#include <cstdlib>

static std::string sparse_temp_name(const char* prefix)
{
    std::ostringstream ss;
    ss << prefix << "_" << std::time(nullptr) << "_" << std::rand();
    return ss.str();
}


/**
 * @brief Sparse-aware matrix multiplication for HDF5 datasets (R6 wrapper).
 *
 * Computes A %*% B using BigDataStatMeth block-wise sparse multiplication,
 * which avoids loading zero-dense blocks.
 * Delegates to BigDataStatMeth::multiplicationSparse().
 *
 * @param ptr_a      External pointer (SEXP) for matrix A.
 * @param ptr_b      External pointer (SEXP) for matrix B.
 * @param block_size Integer block size hint (-1 = auto).
 * @param mix_block  Integer memory block size for parallel path (-1 = auto).
 * @param paral      Logical or NULL; enable OpenMP parallelisation.
 * @param threads    Integer or NULL; thread count.
 * @param compression Integer or NULL; gzip level (NULL = inherit from A).
 * @param outgroup   Character or NULL. Output group. Default \code{"OUTPUT"}.
 * @param outdataset Character or NULL. Output dataset name.
 *      Default \code{"A_sparse_x_B"} where A and B are the input dataset names.
 * @return Named list with elements "filename" and "path" of the result dataset.
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_multiply_sparse(SEXP ptr_a,
                                             SEXP ptr_b,
                                             int  block_size  = -1,
                                             int  mix_block   = -1,
                                             Rcpp::Nullable<bool> paral       = R_NilValue,
                                             Rcpp::Nullable<int>  threads     = R_NilValue,
                                             Rcpp::Nullable<int>  compression = R_NilValue,
                                             Rcpp::Nullable<std::string> outgroup   = R_NilValue,
                                             Rcpp::Nullable<std::string> outdataset = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "",
        Rcpp::Named("path")     = "");

    try {
        H5::Exception::dontPrint();

        auto* rawA = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_a));
        auto* rawB = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_b));

        if (rawA == nullptr || rawB == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!rawA->isOpen() || !rawB->isOpen())
            throw std::runtime_error("Dataset is closed");

        const std::string filename = rawA->getFullPath();
        const std::string groupA   = rawA->getGroup();
        const std::string nameA    = rawA->getDatasetName();
        const std::string groupB   = rawB->getGroup();
        const std::string nameB    = rawB->getDatasetName();

        // Dimension check
        if (rawA->ncols_r() != rawB->nrows_r())
            throw std::runtime_error("Non-conformable matrices for sparse multiply");

        bool bparal = (!paral.isNull()) && Rcpp::as<bool>(paral);
        int comp_level = compression.isNotNull()
                          ? Rcpp::as<int>(compression)
                          : static_cast<int>(rawA->getCompressionLevel());

        // Fresh handles
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        // Auto block size (mirrors bdblockmult_sparse_hdf5 logic)
        Rcpp::Nullable<int> n_block =
            (block_size < 0) ? R_NilValue : Rcpp::wrap(block_size);
        constexpr int iblockfactor = 2;
        hsize_t iblock = static_cast<hsize_t>(
            BigDataStatMeth::getMaxBlockSize(
                static_cast<int>(dsA->nrows()),
                static_cast<int>(dsA->ncols()),
                static_cast<int>(dsB->nrows()),
                static_cast<int>(dsB->ncols()),
                iblockfactor, n_block));

        hsize_t imem = (mix_block < 0)
                        ? (bparal ? iblock / 2 : 0)
                        : static_cast<hsize_t>(mix_block);

        const std::string out_grp  = outgroup.isNull()
            ? std::string("OUTPUT")
                : Rcpp::as<std::string>(outgroup);
        const std::string out_name = outdataset.isNull()
            ? (nameA + "_sparse_x_" + nameB)
            : Rcpp::as<std::string>(outdataset);
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, out_grp, out_name, true));
        
        dsC->setCompressionLevel(comp_level);

        // Delegate (Rule 1.4)
        BigDataStatMeth::multiplicationSparse(
            dsA.get(), dsB.get(), dsC.get(),
            iblock, imem, bparal, /*browmajor=*/true, threads);

        lst["filename"] = filename;
        lst["path"] = out_grp + "/" + out_name;

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_multiply_sparse (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_multiply_sparse (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_multiply_sparse: %s", e.what());
    }

    return lst;
}
