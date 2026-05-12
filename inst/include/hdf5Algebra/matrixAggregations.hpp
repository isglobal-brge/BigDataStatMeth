/**
 * @file matrixAggregations.hpp
 * @brief Block-wise aggregation operations for HDF5 matrices
 *
 * @details
 * Provides block-wise, OpenMP-parallel implementations of common aggregation
 * functions for matrices stored in HDF5 format, following the same design
 * conventions as the rest of the BigDataStatMeth library:
 *
 *  - Column-wise: colSums, colMeans, colMins, colMaxs, colVars, colSds
 *  - Row-wise:    rowSums, rowMeans, rowMins, rowMaxs, rowVars, rowSds
 *  - Scalar:      scalar_sum, scalar_mean, scalar_min, scalar_max,
 *                 scalar_var, scalar_sd  (whole matrix treated as a vector)
 *
 * @section axis_mapping HDF5 ↔ R axis mapping
 *
 * BigDataStatMeth stores R matrices transposed in HDF5:
 *  - `dsA->nrows()` = HDF5 dim[0] = R ncols
 *  - `dsA->ncols()` = HDF5 dim[1] = R nrows
 *  - `dsA->nrows_r()` = R nrows
 *  - `dsA->ncols_r()` = R ncols
 *
 * Column-wise functions (one result per R-column) iterate over HDF5 dim[0]
 * in blocks; row-wise functions iterate over HDF5 dim[1].
 *
 * @section parallelism Parallelism strategy
 *
 * Blocks are precomputed (no off-by-one), then dispatched via
 * `#pragma omp parallel for schedule(dynamic)`.  HDF5 reads are serialised
 * with `#pragma omp critical(accessFile)`.  Aggregation inside each block is
 * done by Eigen and is implicitly vectorised.  For reductions across blocks
 * (scalar operations) OpenMP `reduction` clauses are used.
 *
 * @note
 * This header must be included after BigDataStatMeth.hpp (which defines
 * MAXELEMSINBLOCK and BigDataStatMeth::get_threads).
 */

#ifndef BIGDATASTATMETH_HDF5_MATRIXAGGREGATIONS_HPP
#define BIGDATASTATMETH_HDF5_MATRIXAGGREGATIONS_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"
#include <cmath>
#include <limits>

namespace BigDataStatMeth {


// ===========================================================================
// Internal helpers
// ===========================================================================

/**
 * @brief Compute block start positions and sizes for a single dimension.
 *
 * Fills @p starts and @p sizes so that all blocks together cover [0, dim)
 * exactly, with no off-by-one.  The last block may be smaller than
 * @p block_size if @p dim is not an exact multiple.
 *
 * @param dim         Total size of the dimension to iterate.
 * @param block_size  Nominal block size (>= 1).
 * @param starts      Output: start offsets.
 * @param sizes       Output: block sizes (last may be smaller).
 */
inline void agg_make_blocks(hsize_t dim,
                             hsize_t block_size,
                             std::vector<hsize_t>& starts,
                             std::vector<hsize_t>& sizes)
{
    for (hsize_t off = 0; off < dim; off += block_size) {
        starts.push_back(off);
        sizes.push_back(std::min(block_size, dim - off));
    }
}

/**
 * @brief Choose a block size for aggregation.
 *
 * If the user supplied @p wsize, it is used directly (clamped to
 * [1, iterated]).  Otherwise the block size is set so that one block
 * occupies at most MAXELEMSINBLOCK elements (fixed_dim × block_size).
 *
 * @param wsize     User-supplied block size (Nullable).
 * @param iterated  Size of the dimension being blocked (must be > 0).
 * @param fixed     Size of the dimension that is always read completely.
 * @return          Block size in [1, iterated].
 */
inline hsize_t agg_block_size(Rcpp::Nullable<int> wsize,
                               hsize_t iterated,
                               hsize_t fixed)
{
    hsize_t bs;
    if (wsize.isNotNull()) {
        bs = static_cast<hsize_t>(std::max(1, Rcpp::as<int>(wsize)));
    } else {
        // Default: one block holds at most MAXELEMSINBLOCK elements
        bs = (fixed > 0)
            ? std::max(static_cast<hsize_t>(1),
                       static_cast<hsize_t>(
                           std::ceil(static_cast<double>(MAXELEMSINBLOCK) / fixed)))
            : iterated;
    }
    // Never exceed the actual dimension
    return std::min(bs, iterated);
}

// RowMajor Eigen map alias used throughout this file
using RMMatd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;


// ===========================================================================
// Column-wise aggregations
// Result length = R ncols = dsA->nrows()
// Iterate over HDF5 dim[0]; read all of HDF5 dim[1] each time.
// Block shape: (block_rcols × nrows_R)   → X.rowwise().*  gives one value per R-col
// ===========================================================================

/**
 * @brief Column sums of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to base R `colSums(X)`.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length ncols_R.
 */
inline Eigen::VectorXd get_HDF5_colSums(BigDataStatMeth::hdf5Dataset* dsA,
                                         bool bparal,
                                         Rcpp::Nullable<int> wsize,
                                         Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows(); // = R ncols  (iterated)
        const hsize_t nHDF5cols = dsA->ncols(); // = R nrows  (fixed)

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        Eigen::VectorXd result = Eigen::VectorXd::Zero(nHDF5rows);

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..//{ 
            
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols}, stride, blk, vd.data()); 
            
            //.. 20260325 - remove critical ..//}

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            // rowwise sum: one sum per R-column in this block
            result.segment(starts[bi], sizes[bi]) = X.rowwise().sum();
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colSums (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colSums (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_colSums: ")
                                 + e.what());
    }
}

/**
 * @brief Column means of an HDF5 matrix (block-wise, parallel, single-pass).
 *
 * Equivalent to base R `colMeans(X)`.
 *
 * Single-pass implementation: divides the running sum inside each block
 * immediately, avoiding a second traversal of the dataset.  The division
 * constant (nHDF5cols = nrow_R) is the same for every block, so partial
 * means can simply be summed and the grand mean equals the mean of all
 * element values — which is NOT the case in general, but IS the case here
 * because each block contributes an equal fraction of the rows for every
 * column.  Specifically, the fixed dimension (dim[1] = nrow_R) is always
 * read completely in each block, so X.rowwise().mean() per block already
 * gives the exact column mean for the rows in that block.  Accumulating
 * those with += is correct because every column sees every row exactly once
 * across all blocks (the blocks partition dim[0] = ncol_R, not dim[1]).
 * Therefore we accumulate X.rowwise().mean() over blocks and that is the
 * exact colMeans — no further division needed.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length ncols_R.
 */
inline Eigen::VectorXd get_HDF5_colMeans(BigDataStatMeth::hdf5Dataset* dsA,
                                          bool bparal,
                                          Rcpp::Nullable<int> wsize,
                                          Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows(); // = R ncols (iterated)
        const hsize_t nHDF5cols = dsA->ncols(); // = R nrows (fixed, always complete)

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        Eigen::VectorXd result(nHDF5rows);

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols}, stride, blk, vd.data()); 
            
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            // rowwise().mean() = mean over all R rows for each R column in block.
            // The fixed dimension (nHDF5cols = nrow_R) is always fully loaded, so
            // this is the exact column mean — no accumulation across blocks needed.
            result.segment(starts[bi], sizes[bi]) = X.rowwise().mean();
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colMeans (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colMeans (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_colMeans: ")
                                 + e.what());
    }
}

/**
 * @brief Column minimums of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to `apply(X, 2, min)`.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length ncols_R.
 */
inline Eigen::VectorXd get_HDF5_colMins(BigDataStatMeth::hdf5Dataset* dsA,
                                         bool bparal,
                                         Rcpp::Nullable<int> wsize,
                                         Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();
        const hsize_t nHDF5cols = dsA->ncols();

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        // Initialise to +Inf so the first block always wins
        Eigen::VectorXd result =
            Eigen::VectorXd::Constant(nHDF5rows,
                                      std::numeric_limits<double>::infinity());

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols}, stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            result.segment(starts[bi], sizes[bi]) = X.rowwise().minCoeff();
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colMins (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colMins (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_colMins: ")
                                 + e.what());
    }
}

/**
 * @brief Column maximums of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to `apply(X, 2, max)`.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length ncols_R.
 */
inline Eigen::VectorXd get_HDF5_colMaxs(BigDataStatMeth::hdf5Dataset* dsA,
                                         bool bparal,
                                         Rcpp::Nullable<int> wsize,
                                         Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();
        const hsize_t nHDF5cols = dsA->ncols();

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        // Initialise to -Inf so the first block always wins
        Eigen::VectorXd result =
            Eigen::VectorXd::Constant(nHDF5rows,
                                      -std::numeric_limits<double>::infinity());

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            // #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols},stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            result.segment(starts[bi], sizes[bi]) = X.rowwise().maxCoeff();
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colMaxs (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colMaxs (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_colMaxs: ")
                                 + e.what());
    }
}

/**
 * @brief Column variances of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to `apply(X, 2, var)` — uses Bessel's correction (n-1).
 * If nrow_R == 1 the result is a vector of NAs, matching base R behaviour.
 *
 * Each block loads all nrow_R values for a set of columns, so the variance
 * is computed exactly in a single pass per column (no merging required).
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length ncols_R.
 */
inline Eigen::VectorXd get_HDF5_colVars(BigDataStatMeth::hdf5Dataset* dsA,
                                         bool bparal,
                                         Rcpp::Nullable<int> wsize,
                                         Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();  // R ncols (iterated)
        const hsize_t nHDF5cols = dsA->ncols();  // R nrows (fixed)
        const double n = static_cast<double>(nHDF5cols);

        // var undefined for n < 2 — return NaN vector (same as R)
        if (nHDF5cols < 2) {
            return Eigen::VectorXd::Constant(nHDF5rows,
                                             std::numeric_limits<double>::quiet_NaN());
        }

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        Eigen::VectorXd result(nHDF5rows);

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols},stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            // Computational formula: var = (sum_sq - sum^2/n) / (n-1)
            // Each row of X corresponds to one R-column (all n R-rows present)
            const Eigen::VectorXd colsum   = X.rowwise().sum();
            const Eigen::VectorXd colsumsq = X.rowwise().squaredNorm();

            result.segment(starts[bi], sizes[bi]) =
                (colsumsq.array() - colsum.array().square() / n) / (n - 1.0);
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colVars (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_colVars (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_colVars: ")
                                 + e.what());
    }
}

/**
 * @brief Column standard deviations of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to `apply(X, 2, sd)` — uses Bessel's correction (n-1).
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length ncols_R.
 */
inline Eigen::VectorXd get_HDF5_colSds(BigDataStatMeth::hdf5Dataset* dsA,
                                        bool bparal,
                                        Rcpp::Nullable<int> wsize,
                                        Rcpp::Nullable<int> threads)
{
    try {
        return get_HDF5_colVars(dsA, bparal, wsize, threads).array().sqrt();
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_colSds: ")
                                 + e.what());
    }
}


// ===========================================================================
// Row-wise aggregations
// Result length = R nrows = dsA->ncols()
// Iterate over HDF5 dim[1]; read all of HDF5 dim[0] each time.
// Block shape: (ncols_R × block_rrows)   → X.colwise().*  gives one value per R-row
// ===========================================================================

/**
 * @brief Row sums of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to base R `rowSums(X)`.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length nrows_R.
 */
inline Eigen::VectorXd get_HDF5_rowSums(BigDataStatMeth::hdf5Dataset* dsA,
                                         bool bparal,
                                         Rcpp::Nullable<int> wsize,
                                         Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows(); // = R ncols  (fixed)
        const hsize_t nHDF5cols = dsA->ncols(); // = R nrows  (iterated)

        const hsize_t bs = agg_block_size(wsize, nHDF5cols, nHDF5rows);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5cols, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        Eigen::VectorXd result = Eigen::VectorXd::Zero(nHDF5cols);

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(nHDF5rows * sizes[bi]);
            // #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({0, starts[bi]}, {nHDF5rows, sizes[bi]},stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            // Map as (ncols_R × block_rrows) RowMajor
            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(nHDF5rows),
                static_cast<Eigen::Index>(sizes[bi]));

            // colwise sum: one sum per R-row in this block
            result.segment(starts[bi], sizes[bi]) = X.colwise().sum().transpose();
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowSums (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowSums (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_rowSums: ")
                                 + e.what());
    }
}

/**
 * @brief Row means of an HDF5 matrix (block-wise, parallel, single-pass).
 *
 * Equivalent to base R `rowMeans(X)`.
 *
 * Single-pass: each block loads all ncol_R values for a contiguous set of
 * rows (the fixed dimension nHDF5rows = ncol_R is always read complete).
 * Therefore X.colwise().mean() gives the exact row mean for every R-row in
 * the block in one shot — no second pass required.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length nrows_R.
 */
inline Eigen::VectorXd get_HDF5_rowMeans(BigDataStatMeth::hdf5Dataset* dsA,
                                          bool bparal,
                                          Rcpp::Nullable<int> wsize,
                                          Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows(); // = R ncols (fixed)
        const hsize_t nHDF5cols = dsA->ncols(); // = R nrows (iterated)

        const hsize_t bs = agg_block_size(wsize, nHDF5cols, nHDF5rows);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5cols, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        Eigen::VectorXd result(nHDF5cols);

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(nHDF5rows * sizes[bi]);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({0, starts[bi]}, {nHDF5rows, sizes[bi]},stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            // Map as (ncols_R × block_rrows) RowMajor
            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(nHDF5rows),
                static_cast<Eigen::Index>(sizes[bi]));

            // colwise().mean() = mean over all R cols for each R row in this block.
            // nHDF5rows = ncol_R is always fully loaded → exact row mean.
            result.segment(starts[bi], sizes[bi]) =
                X.colwise().mean().transpose();
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowMeans (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowMeans (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_rowMeans: ")
                                 + e.what());
    }
}

/**
 * @brief Row minimums of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to `apply(X, 1, min)`.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length nrows_R.
 */
inline Eigen::VectorXd get_HDF5_rowMins(BigDataStatMeth::hdf5Dataset* dsA,
                                         bool bparal,
                                         Rcpp::Nullable<int> wsize,
                                         Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();
        const hsize_t nHDF5cols = dsA->ncols();

        const hsize_t bs = agg_block_size(wsize, nHDF5cols, nHDF5rows);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5cols, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        Eigen::VectorXd result =
            Eigen::VectorXd::Constant(nHDF5cols,
                                      std::numeric_limits<double>::infinity());

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(nHDF5rows * sizes[bi]);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({0, starts[bi]}, {nHDF5rows, sizes[bi]}, stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(nHDF5rows),
                static_cast<Eigen::Index>(sizes[bi]));

            // colwise min: compares across ncols_R (HDF5 rows), one min per R-row
            Eigen::RowVectorXd block_mins = X.colwise().minCoeff();
            // Take element-wise min with running result to handle multi-block case
            // (This case only occurs when ncols_R > 1 block, which shouldn't
            //  happen for row-wise iteration where ncols_R is the fixed dim.
            //  But defensive code here for correctness.)
            for (hsize_t k = 0; k < sizes[bi]; k++) {
                result[starts[bi] + k] = std::min(result[starts[bi] + k],
                                                   block_mins[k]);
            }
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowMins (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowMins (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_rowMins: ")
                                 + e.what());
    }
}

/**
 * @brief Row maximums of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to `apply(X, 1, max)`.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length nrows_R.
 */
inline Eigen::VectorXd get_HDF5_rowMaxs(BigDataStatMeth::hdf5Dataset* dsA,
                                         bool bparal,
                                         Rcpp::Nullable<int> wsize,
                                         Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();
        const hsize_t nHDF5cols = dsA->ncols();

        const hsize_t bs = agg_block_size(wsize, nHDF5cols, nHDF5rows);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5cols, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        Eigen::VectorXd result =
            Eigen::VectorXd::Constant(nHDF5cols,
                                      -std::numeric_limits<double>::infinity());

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(nHDF5rows * sizes[bi]);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({0, starts[bi]}, {nHDF5rows, sizes[bi]},stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(nHDF5rows),
                static_cast<Eigen::Index>(sizes[bi]));

            Eigen::RowVectorXd block_maxs = X.colwise().maxCoeff();
            for (hsize_t k = 0; k < sizes[bi]; k++) {
                result[starts[bi] + k] = std::max(result[starts[bi] + k],
                                                   block_maxs[k]);
            }
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowMaxs (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowMaxs (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_rowMaxs: ")
                                 + e.what());
    }
}

/**
 * @brief Row variances of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to `apply(X, 1, var)` — uses Bessel's correction (n-1).
 * If ncol_R == 1 the result is a vector of NaNs, matching base R behaviour.
 *
 * Each block loads all ncol_R values for a set of rows, so the variance
 * is computed exactly in a single pass per row.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length nrows_R.
 */
inline Eigen::VectorXd get_HDF5_rowVars(BigDataStatMeth::hdf5Dataset* dsA,
                                         bool bparal,
                                         Rcpp::Nullable<int> wsize,
                                         Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows(); // = R ncols  (fixed)
        const hsize_t nHDF5cols = dsA->ncols(); // = R nrows  (iterated)
        const double n = static_cast<double>(nHDF5rows);

        // var undefined for n < 2
        if (nHDF5rows < 2) {
            return Eigen::VectorXd::Constant(nHDF5cols,
                                             std::numeric_limits<double>::quiet_NaN());
        }

        const hsize_t bs = agg_block_size(wsize, nHDF5cols, nHDF5rows);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5cols, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        Eigen::VectorXd result(nHDF5cols);

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes, result)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(nHDF5rows * sizes[bi]);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
                dsA->readDatasetBlock({0, starts[bi]}, {nHDF5rows, sizes[bi]}, stride, blk, vd.data()); 
                //.. 20260325 - remove critical ..// }

            // Map as (ncols_R × block_rrows) RowMajor
            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(nHDF5rows),
                static_cast<Eigen::Index>(sizes[bi]));

            // Computational formula over R-rows in this block:
            //   var_row = (sum_sq_row - sum_row^2 / n) / (n - 1)
            // where n = ncols_R (HDF5 nrows, all loaded)
            const Eigen::RowVectorXd rowsum   = X.colwise().sum();
            const Eigen::RowVectorXd rowsumsq = X.colwise().squaredNorm();

            result.segment(starts[bi], sizes[bi]) =
                ((rowsumsq.array() - rowsum.array().square() / n) /
                 (n - 1.0)).transpose();
        }

        return result;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowVars (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_rowVars (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_rowVars: ")
                                 + e.what());
    }
}

/**
 * @brief Row standard deviations of an HDF5 matrix (block-wise, parallel).
 *
 * Equivalent to `apply(X, 1, sd)` — uses Bessel's correction (n-1).
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Vector of length nrows_R.
 */
inline Eigen::VectorXd get_HDF5_rowSds(BigDataStatMeth::hdf5Dataset* dsA,
                                        bool bparal,
                                        Rcpp::Nullable<int> wsize,
                                        Rcpp::Nullable<int> threads)
{
    try {
        return get_HDF5_rowVars(dsA, bparal, wsize, threads).array().sqrt();
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_rowSds: ")
                                 + e.what());
    }
}


// ===========================================================================
// Scalar aggregations  (whole matrix treated as a flat vector)
// Iterate column-wise (over HDF5 dim[0]) for memory efficiency.
// ===========================================================================

/**
 * @brief Sum of all elements of an HDF5 matrix.
 *
 * Equivalent to `sum(X)`.  Uses OMP reduction over column blocks.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Scalar sum.
 */
inline double get_HDF5_scalar_sum(BigDataStatMeth::hdf5Dataset* dsA,
                                   bool bparal,
                                   Rcpp::Nullable<int> wsize,
                                   Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();
        const hsize_t nHDF5cols = dsA->ncols();

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        double total_sum = 0.0;

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes) reduction(+:total_sum)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols}, stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            total_sum += X.sum();
        }

        return total_sum;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_scalar_sum (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_scalar_sum (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_scalar_sum: ")
                                 + e.what());
    }
}

/**
 * @brief Mean of all elements of an HDF5 matrix.
 *
 * Equivalent to `mean(X)`.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Scalar mean.
 */
inline double get_HDF5_scalar_mean(BigDataStatMeth::hdf5Dataset* dsA,
                                    bool bparal,
                                    Rcpp::Nullable<int> wsize,
                                    Rcpp::Nullable<int> threads)
{
    try {
        const double N = static_cast<double>(dsA->nrows()) *
                         static_cast<double>(dsA->ncols());
        return get_HDF5_scalar_sum(dsA, bparal, wsize, threads) / N;
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_scalar_mean: ")
                                 + e.what());
    }
}

/**
 * @brief Minimum of all elements of an HDF5 matrix.
 *
 * Equivalent to `min(X)`.  Uses OMP reduction (min clause).
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Scalar minimum.
 */
inline double get_HDF5_scalar_min(BigDataStatMeth::hdf5Dataset* dsA,
                                   bool bparal,
                                   Rcpp::Nullable<int> wsize,
                                   Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();
        const hsize_t nHDF5cols = dsA->ncols();

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        double global_min = std::numeric_limits<double>::infinity();

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes) reduction(min:global_min)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols}, stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            global_min = std::min(global_min, X.minCoeff());
        }

        return global_min;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_scalar_min (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_scalar_min (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_scalar_min: ")
                                 + e.what());
    }
}

/**
 * @brief Maximum of all elements of an HDF5 matrix.
 *
 * Equivalent to `max(X)`.  Uses OMP reduction (max clause).
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Scalar maximum.
 */
inline double get_HDF5_scalar_max(BigDataStatMeth::hdf5Dataset* dsA,
                                   bool bparal,
                                   Rcpp::Nullable<int> wsize,
                                   Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();
        const hsize_t nHDF5cols = dsA->ncols();

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        double global_max = -std::numeric_limits<double>::infinity();

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes) reduction(max:global_max)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols}, stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            global_max = std::max(global_max, X.maxCoeff());
        }

        return global_max;

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_scalar_max (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_scalar_max (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_scalar_max: ")
                                 + e.what());
    }
}

/**
 * @brief Variance of all elements of an HDF5 matrix (treated as a flat vector).
 *
 * Equivalent to `var(as.vector(X))` — uses Bessel's correction (N-1).
 *
 * Algorithm: single-pass computational formula
 *   var = (sum_sq - sum^2 / N) / (N - 1)
 * accumulated over column blocks using OMP reductions.
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Scalar variance.
 */
inline double get_HDF5_scalar_var(BigDataStatMeth::hdf5Dataset* dsA,
                                   bool bparal,
                                   Rcpp::Nullable<int> wsize,
                                   Rcpp::Nullable<int> threads)
{
    try {
        const hsize_t nHDF5rows = dsA->nrows();
        const hsize_t nHDF5cols = dsA->ncols();
        const double N = static_cast<double>(nHDF5rows) *
                         static_cast<double>(nHDF5cols);

        if (N < 2.0)
            return std::numeric_limits<double>::quiet_NaN();

        const hsize_t bs = agg_block_size(wsize, nHDF5rows, nHDF5cols);

        std::vector<hsize_t> starts, sizes;
        agg_make_blocks(nHDF5rows, bs, starts, sizes);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};
        const int nthreads = static_cast<int>(
            BigDataStatMeth::get_threads(bparal, threads));

        double total_sum   = 0.0;
        double total_sumsq = 0.0;

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
                shared(dsA, starts, sizes) \
                reduction(+:total_sum, total_sumsq)
        for (hsize_t bi = 0; bi < starts.size(); bi++) {
            std::vector<double> vd(sizes[bi] * nHDF5cols);
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            //.. 20260325 - remove critical ..// { 
            dsA->readDatasetBlock({starts[bi], 0}, {sizes[bi], nHDF5cols}, stride, blk, vd.data()); 
            //.. 20260325 - remove critical ..// }

            Eigen::Map<const RMMatd> X(vd.data(),
                static_cast<Eigen::Index>(sizes[bi]),
                static_cast<Eigen::Index>(nHDF5cols));

            total_sum   += X.sum();
            total_sumsq += X.array().square().sum();
        }

        // Computational formula: var = (sum_sq - sum^2/N) / (N-1)
        return (total_sumsq - total_sum * total_sum / N) / (N - 1.0);

    } catch (H5::FileIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_scalar_var (File IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (H5::DataSetIException& e) {
        throw std::runtime_error("c++ exception get_HDF5_scalar_var (DataSet IException): "
                                 + std::string(e.getDetailMsg()));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_scalar_var: ")
                                 + e.what());
    }
}

/**
 * @brief Standard deviation of all elements of an HDF5 matrix.
 *
 * Equivalent to `sd(as.vector(X))` — uses Bessel's correction (N-1).
 *
 * @param dsA     Open HDF5 dataset.
 * @param bparal  Enable OpenMP parallelism.
 * @param wsize   Block size (NULL = auto).
 * @param threads Thread count (NULL = auto).
 * @return        Scalar standard deviation.
 */
inline double get_HDF5_scalar_sd(BigDataStatMeth::hdf5Dataset* dsA,
                                  bool bparal,
                                  Rcpp::Nullable<int> wsize,
                                  Rcpp::Nullable<int> threads)
{
    try {
        return std::sqrt(get_HDF5_scalar_var(dsA, bparal, wsize, threads));
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("c++ exception get_HDF5_scalar_sd: ")
                                 + e.what());
    }
}

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_HDF5_MATRIXAGGREGATIONS_HPP
