/**
 * @file hdf5_r6_arithmetic.cpp
 * @brief R6 wrappers for elementwise arithmetic operations on HDF5 datasets
 *
 * @details
 * Provides Rcpp wrappers for four block-wise elementwise arithmetic operations
 * used by the HDF5Matrix R6/S3 interface:
 *
 *  - \c rcpp_hdf5dataset_add      : element-wise addition       (A + B)
 *  - \c rcpp_hdf5dataset_subtract : element-wise subtraction    (A - B)
 *  - \c rcpp_hdf5dataset_mul_ew   : element-wise multiplication (A * B, Hadamard)
 *  - \c rcpp_hdf5dataset_div_ew   : element-wise division       (A / B)
 *
 * Add and subtract delegate to the existing headers:
 *  - \c BigDataStatMeth::Rcpp_block_matrix_sum_hdf5        (matrixSum.hpp)
 *  - \c BigDataStatMeth::Rcpp_block_matrix_substract_hdf5  (matrixSubstract.hpp)
 *
 * Multiply and divide follow the identical block-wise pattern (read block A,
 * read block B, std::transform with std::multiplies / std::divides, write
 * block C) implemented directly in this file, as no dedicated header exists.
 *
 * All four functions follow the same conventions as hdf5_r6_multiply.cpp:
 *  - inputs are external pointers (SEXP) to open hdf5Dataset objects
 *  - ALL logic inside try{} block
 *  - validation throws std::runtime_error (never Rcpp::stop inside try{})
 *  - C++ objects managed with std::unique_ptr (exception-safe)
 *  - result written to the same HDF5 file as A, named with a unique temp name
 *  - location returned as List(filename, path)
 *
 * @note
 * CRITICAL: Rcpp::stop() must only be called inside catch{} handlers.
 * Use Rf_error() in catch blocks (not Rcpp::stop()) to avoid memory
 * corruption (double-free / SIGABRT) on macOS ARM64 with Rcpp >= 1.1.0.
 */

#include <BigDataStatMeth.hpp>
#include <functional>
#include <sstream>
#include <ctime>
#include <cstdlib>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

static std::string arith_temp_name(const char* prefix)
{
    std::ostringstream ss;
    ss << prefix << "_" << std::time(nullptr) << "_" << std::rand();
    return ss.str();
}

/**
 * @brief Resolve block size for elementwise operations.
 *
 * Mirrors the logic used in bdblockSum_hdf5: if the user supplies a value
 * use it directly; otherwise auto-calculate via getMatrixBlockSize /
 * getVectorBlockSize depending on whether A is a vector or a matrix.
 */
static int resolve_block_size( BigDataStatMeth::hdf5Dataset* dsA,
                               Rcpp::Nullable<int> block_size)
{
    if (block_size.isNotNull())
        return Rcpp::as<int>(block_size);

    const hsize_t nR = dsA->nrows();
    const hsize_t nC = dsA->ncols();

    if (nR == 1 || nC == 1) {
        return static_cast<int>(BigDataStatMeth::getVectorBlockSize(
            static_cast<int>(nR * nC)));
    }

    std::vector<hsize_t> bs = BigDataStatMeth::getMatrixBlockSize(
        static_cast<int>(nR), static_cast<int>(nC));
    return static_cast<int>(nR < nC ? bs.at(0) : bs.at(1));
}

/**
 * @brief Block-wise elementwise transform for two HDF5 matrix datasets.
 *
 * Reads corresponding blocks from dsA and dsB, applies a binary functor
 * element-wise, and writes the result to dsC.  The same adaptive block
 * strategy as matrixSum.hpp is used (iterate over columns when ncols >= nrows,
 * over rows otherwise).
 *
 * @tparam BinaryOp  A callable compatible with std::transform binary form
 *                   (e.g. std::plus<double>(), std::multiplies<double>())
 */
template<typename BinaryOp>
static void block_elementwise_hdf5(BigDataStatMeth::hdf5Dataset* dsA,
                                    BigDataStatMeth::hdf5Dataset* dsB,
                                    BigDataStatMeth::hdf5Dataset* dsC,
                                    hsize_t hdf5_block,
                                    bool    bparal,
                                    Rcpp::Nullable<int> threads,
                                    BinaryOp op)
{
    const hsize_t K = dsA->nrows();
    const hsize_t N = dsA->ncols();

    if (K != dsB->nrows() || N != dsB->ncols())
        throw std::runtime_error("non-conformable arguments");

    dsC->createDataset(N, K, "real");

    const std::vector<hsize_t> stride = {1, 1};
    const std::vector<hsize_t> blk    = {1, 1};

    std::vector<hsize_t> vstart, vsize;

    if (K <= N) {
        BigDataStatMeth::getBlockPositionsSizes(N, hdf5_block, vstart, vsize);

        #pragma omp parallel num_threads(BigDataStatMeth::get_threads(bparal, threads)) \
                shared(dsA, dsB, dsC)
        {
        #pragma omp for schedule(dynamic)
        for (hsize_t ii = 0; ii < vstart.size(); ii++) {
            std::vector<double> vdA(K * vsize[ii]);
            #pragma omp critical(accessFile)
            { dsA->readDatasetBlock({0, vstart[ii]}, {K, vsize[ii]},
                                    stride, blk, vdA.data()); }

            std::vector<double> vdB(K * vsize[ii]);
            #pragma omp critical(accessFile)
            { dsB->readDatasetBlock({0, vstart[ii]}, {K, vsize[ii]},
                                    stride, blk, vdB.data()); }

            std::transform(vdA.begin(), vdA.end(), vdB.begin(), vdA.begin(), op);

            #pragma omp critical(accessFile)
            { dsC->writeDatasetBlock(vdA, {0, vstart[ii]}, {K, vsize[ii]},
                                     stride, blk); }
        }
        }
    } else {
        BigDataStatMeth::getBlockPositionsSizes(K, hdf5_block, vstart, vsize);

        #pragma omp parallel num_threads(BigDataStatMeth::get_threads(bparal, threads)) \
                shared(dsA, dsB, dsC)
        {
        #pragma omp for schedule(dynamic)
        for (hsize_t ii = 0; ii < vstart.size(); ii++) {
            std::vector<double> vdA(vsize[ii] * N);
            #pragma omp critical(accessFile)
            { dsA->readDatasetBlock({vstart[ii], 0}, {vsize[ii], N},
                                    stride, blk, vdA.data()); }

            std::vector<double> vdB(vsize[ii] * N);
            #pragma omp critical(accessFile)
            { dsB->readDatasetBlock({vstart[ii], 0}, {vsize[ii], N},
                                    stride, blk, vdB.data()); }

            std::transform(vdA.begin(), vdA.end(), vdB.begin(), vdA.begin(), op);

            #pragma omp critical(accessFile)
            { dsC->writeDatasetBlock(vdA, {vstart[ii], 0}, {vsize[ii], N},
                                     stride, blk); }
        }
        }
    }
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_add  (A + B)
// ---------------------------------------------------------------------------

//' Element-wise addition of two HDF5 datasets (R6 wrapper)
//'
//' @description
//' Computes \code{A + B} element-wise for two HDF5 datasets referenced by
//' external pointers, using a block-wise algorithm.
//'
//' @param ptr_a External pointer (SEXP) for matrix A
//' @param ptr_b External pointer (SEXP) for matrix B
//' @param paral Logical or NULL; enable OpenMP parallelisation
//' @param block_size Integer or NULL; block size (NULL = auto)
//' @param threads Integer or NULL; thread count
//'
//' @return Named list with \code{filename} and \code{path} of the result.
//'   The result is stored in group \code{"OUTPUT"} with dataset name
//'   \code{"A_plus_B"} (resp. \code{"A_minus_B"}, \code{"A_times_B"},
//'   \code{"A_div_B"}) where A and B are the input dataset names.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_add(SEXP ptr_a,
                                SEXP ptr_b,
                                Rcpp::Nullable<bool> paral      = R_NilValue,
                                Rcpp::Nullable<int>  block_size = R_NilValue,
                                Rcpp::Nullable<int>  threads    = R_NilValue,
                                Rcpp::Nullable<int>  compression = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "", Rcpp::Named("path") = "");

    try {
        H5::Exception::dontPrint();

        auto* dsA_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_a));
        auto* dsB_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_b));

        if (dsA_raw == nullptr || dsB_raw == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!dsA_raw->isOpen() || !dsB_raw->isOpen())
            throw std::runtime_error("Dataset is closed");

        //.. 20260304 ..//  const std::string filename = dsA_raw->getFileName();
        const std::string filename = dsA_raw->getFullPath();
        const std::string groupA   = dsA_raw->getGroup();
        const std::string nameA    = dsA_raw->getDatasetName();
        const std::string groupB   = dsB_raw->getGroup();
        const std::string nameB    = dsB_raw->getDatasetName();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        const int iblock_size = resolve_block_size(dsA.get(), block_size);
        bool bparal = paral.isNotNull() ? Rcpp::as<bool>(paral) : false;

        const std::string result_name = nameA + "_plus_" + nameB;
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, "OUTPUT", result_name, true));
        dsC->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());

        BigDataStatMeth::Rcpp_block_matrix_sum_hdf5(
            dsA.get(), dsB.get(), dsC.get(),
            static_cast<hsize_t>(iblock_size), bparal, threads);

        lst["filename"] = filename;
        lst["path"] = "OUTPUT/" + result_name;

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error in add: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error in add: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error in add: %s", e.what());
    }

    return lst;
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_subtract  (A - B)
// ---------------------------------------------------------------------------

//' Element-wise subtraction of two HDF5 datasets (R6 wrapper)
//'
//' @description
//' Computes \code{A - B} element-wise for two HDF5 datasets referenced by
//' external pointers, using a block-wise algorithm.
//'
//' @param ptr_a External pointer (SEXP) for matrix A
//' @param ptr_b External pointer (SEXP) for matrix B
//' @param paral Logical or NULL; enable OpenMP parallelisation
//' @param block_size Integer or NULL; block size (NULL = auto)
//' @param threads Integer or NULL; thread count
//'
//' @return Named list with \code{filename} and \code{path} of the result.
//'   The result is stored in group \code{"OUTPUT"} with dataset name
//'   \code{"A_plus_B"} (resp. \code{"A_minus_B"}, \code{"A_times_B"},
//'   \code{"A_div_B"}) where A and B are the input dataset names.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_subtract(SEXP ptr_a,
                                     SEXP ptr_b,
                                     Rcpp::Nullable<bool> paral      = R_NilValue,
                                     Rcpp::Nullable<int>  block_size = R_NilValue,
                                     Rcpp::Nullable<int>  threads    = R_NilValue,
                                Rcpp::Nullable<int>  compression = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "", Rcpp::Named("path") = "");

    try {
        H5::Exception::dontPrint();

        auto* dsA_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_a));
        auto* dsB_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_b));

        if (dsA_raw == nullptr || dsB_raw == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!dsA_raw->isOpen() || !dsB_raw->isOpen())
            throw std::runtime_error("Dataset is closed");

        //.. 20260304 ..//  const std::string filename = dsA_raw->getFileName();
        const std::string filename = dsA_raw->getFullPath();
        const std::string groupA   = dsA_raw->getGroup();
        const std::string nameA    = dsA_raw->getDatasetName();
        const std::string groupB   = dsB_raw->getGroup();
        const std::string nameB    = dsB_raw->getDatasetName();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        const int iblock_size = resolve_block_size(dsA.get(), block_size);
        bool bparal = paral.isNotNull() ? Rcpp::as<bool>(paral) : false;

        const std::string result_name = nameA + "_minus_" + nameB;
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, "OUTPUT", result_name, true));
        dsC->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());

        BigDataStatMeth::Rcpp_block_matrix_substract_hdf5(
            dsA.get(), dsB.get(), dsC.get(),
            static_cast<hsize_t>(iblock_size), bparal, threads);

        lst["filename"] = filename;
        lst["path"] = "OUTPUT/" + result_name;

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error in subtract: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error in subtract: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error in subtract: %s", e.what());
    }

    return lst;
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_mul_ew  (A * B, Hadamard product)
// ---------------------------------------------------------------------------

//' Element-wise multiplication of two HDF5 datasets (R6 wrapper)
//'
//' @description
//' Computes the Hadamard (element-wise) product \code{A * B} for two HDF5
//' datasets referenced by external pointers, using a block-wise algorithm.
//'
//' @param ptr_a External pointer (SEXP) for matrix A
//' @param ptr_b External pointer (SEXP) for matrix B
//' @param paral Logical or NULL; enable OpenMP parallelisation
//' @param block_size Integer or NULL; block size (NULL = auto)
//' @param threads Integer or NULL; thread count
//'
//' @return Named list with \code{filename} and \code{path} of the result.
//'   The result is stored in group \code{"OUTPUT"} with dataset name
//'   \code{"A_plus_B"} (resp. \code{"A_minus_B"}, \code{"A_times_B"},
//'   \code{"A_div_B"}) where A and B are the input dataset names.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_mul_ew(SEXP ptr_a,
                                   SEXP ptr_b,
                                   Rcpp::Nullable<bool> paral      = R_NilValue,
                                   Rcpp::Nullable<int>  block_size = R_NilValue,
                                   Rcpp::Nullable<int>  threads    = R_NilValue,
                                Rcpp::Nullable<int>  compression = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "", Rcpp::Named("path") = "");

    try {
        H5::Exception::dontPrint();

        auto* dsA_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_a));
        auto* dsB_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_b));

        if (dsA_raw == nullptr || dsB_raw == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!dsA_raw->isOpen() || !dsB_raw->isOpen())
            throw std::runtime_error("Dataset is closed");

        //.. 20260304 ..//  const std::string filename = dsA_raw->getFileName();
        const std::string filename = dsA_raw->getFullPath();
        const std::string groupA   = dsA_raw->getGroup();
        const std::string nameA    = dsA_raw->getDatasetName();
        const std::string groupB   = dsB_raw->getGroup();
        const std::string nameB    = dsB_raw->getDatasetName();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        const int iblock_size = resolve_block_size(dsA.get(), block_size);
        bool bparal = paral.isNotNull() ? Rcpp::as<bool>(paral) : false;

        const std::string result_name = nameA + "_times_" + nameB;
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, "OUTPUT", result_name, true));
        dsC->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());

        block_elementwise_hdf5(dsA.get(), dsB.get(), dsC.get(),
                                static_cast<hsize_t>(iblock_size),
                                bparal, threads,
                                std::multiplies<double>());

        lst["filename"] = filename;
        lst["path"] = "OUTPUT/" + result_name;

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error in multiply_ew: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error in multiply_ew: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error in multiply_ew: %s", e.what());
    }

    return lst;
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_div_ew  (A / B)
// ---------------------------------------------------------------------------

//' Element-wise division of two HDF5 datasets (R6 wrapper)
//'
//' @description
//' Computes \code{A / B} element-wise for two HDF5 datasets referenced by
//' external pointers, using a block-wise algorithm.  Division by zero
//' produces \code{NaN} or \code{Inf}, matching base R behaviour.
//'
//' @param ptr_a External pointer (SEXP) for matrix A
//' @param ptr_b External pointer (SEXP) for matrix B
//' @param paral Logical or NULL; enable OpenMP parallelisation
//' @param block_size Integer or NULL; block size (NULL = auto)
//' @param threads Integer or NULL; thread count
//'
//' @return Named list with \code{filename} and \code{path} of the result.
//'   The result is stored in group \code{"OUTPUT"} with dataset name
//'   \code{"A_plus_B"} (resp. \code{"A_minus_B"}, \code{"A_times_B"},
//'   \code{"A_div_B"}) where A and B are the input dataset names.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_div_ew(SEXP ptr_a,
                                   SEXP ptr_b,
                                   Rcpp::Nullable<bool> paral      = R_NilValue,
                                   Rcpp::Nullable<int>  block_size = R_NilValue,
                                   Rcpp::Nullable<int>  threads    = R_NilValue,
                                Rcpp::Nullable<int>  compression = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "", Rcpp::Named("path") = "");

    try {
        H5::Exception::dontPrint();

        auto* dsA_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_a));
        auto* dsB_raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_b));

        if (dsA_raw == nullptr || dsB_raw == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!dsA_raw->isOpen() || !dsB_raw->isOpen())
            throw std::runtime_error("Dataset is closed");

        //.. 20260304 ..//  const std::string filename = dsA_raw->getFileName();
        const std::string filename = dsA_raw->getFullPath();
        const std::string groupA   = dsA_raw->getGroup();
        const std::string nameA    = dsA_raw->getDatasetName();
        const std::string groupB   = dsB_raw->getGroup();
        const std::string nameB    = dsB_raw->getDatasetName();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        const int iblock_size = resolve_block_size(dsA.get(), block_size);
        bool bparal = paral.isNotNull() ? Rcpp::as<bool>(paral) : false;

        const std::string result_name = nameA + "_div_" + nameB;
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, "OUTPUT", result_name, true));
        dsC->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());

        block_elementwise_hdf5(dsA.get(), dsB.get(), dsC.get(),
                                static_cast<hsize_t>(iblock_size),
                                bparal, threads,
                                std::divides<double>());

        lst["filename"] = filename;
        lst["path"] = "OUTPUT/" + result_name;

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error in divide_ew: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error in divide_ew: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error in divide_ew: %s", e.what());
    }

    return lst;
}
