/**
 * @file crossprod.hpp
 * @brief Implementation of matrix cross-product operations using HDF5
 * @note 2026-03-07 Output datasets now inherit compression level from input datasets
 *         via setCompressionLevel() called before every createDataset() invocation.
 *
 * This file provides functionality for computing matrix cross-products (A^T * B)
 * for large matrices stored in HDF5 format. The implementation uses block-wise
 * operations to handle matrices that don't fit in memory, with optimizations
 * for performance and memory usage.
 *
 * Key features:
 * - Block-wise matrix cross-product computation
 * - Memory-efficient implementation for large matrices
 * - Optimized block size selection
 * - Support for row-major and column-major storage
 * - Optional parallel processing
 *
 * @note This implementation is particularly suited for large matrices
 * where traditional in-memory operations are not feasible.
 *
 * @see BigDataStatMeth::hdf5Dataset
 */

#ifndef BIGDATASTATMETH_ALGEBRA_CROSSPROD_HPP
#define BIGDATASTATMETH_ALGEBRA_CROSSPROD_HPP

// #include <RcppEigen.h>
// #include "H5Cpp.h"

namespace BigDataStatMeth {

    /**
     * @brief Computes the cross-product of two matrices stored in HDF5 format
     *
     * @param dsA First input matrix dataset (A)
     * @param dsB Second input matrix dataset (B)
     * @param dsC Output matrix dataset for result
     * @param hdf5_block Block size for HDF5 operations
     * @param mem_block_size Memory block size for computations
     * @param bparal Whether to use parallel processing
     * @param browmajor Whether matrices are stored in row-major order
     * @param threads Number of threads for parallel processing (optional)
     * @return BigDataStatMeth::hdf5Dataset* Pointer to the result dataset
     *
     * @details Computes C = A^T * B using block-wise operations:
     * - Divides matrices into blocks for memory-efficient processing
     * - Processes blocks using optimized Eigen operations
     * - Accumulates results in the output dataset
     *
     * Implementation notes:
     * - Block sizes are automatically adjusted for matrix boundaries
     * - Uses Eigen for efficient block matrix operations
     * - Handles row-major and column-major storage formats
     *
     * Performance considerations:
     * - Time complexity: O(N*M*K) for matrices of sizes N×K and M×K
     * - Space complexity: O(block_size²) for block operations
     * - I/O complexity depends on block size and matrix dimensions
     *
     * Optimization parameters:
     * - hdf5_block: Controls HDF5 read/write block size
     * - mem_block_size: Controls in-memory block processing size
     * - Block sizes are adjusted for optimal performance
     *
     * Example usage:
     * @code
     * auto* A = hdf5Dataset(...);  // N×K matrix
     * auto* B = hdf5Dataset(...);  // M×K matrix
     * auto* C = hdf5Dataset(...);  // N×M result matrix
     * crossprod(A, B, C, 1000, 1000, true, true);
     * @endcode
     *
     * @throws std::range_error if matrix dimensions are incompatible
     * @throws std::exception on HDF5 operations errors
     *
     * @see getOptimBlockSize for block size optimization
     */
    inline BigDataStatMeth::hdf5Dataset* crossprod( 
            BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, 
            BigDataStatMeth::hdf5Dataset* dsC, bool isSymmetric, 
            hsize_t hdf5_block, hsize_t mem_block_size, bool bparal, 
            bool browmajor, Rcpp::Nullable<int> threads = R_NilValue) 

    {
        
        try {
            
            hsize_t N = dsA->nrows();
            hsize_t K = dsA->ncols();
            hsize_t M = dsB->nrows();
            hsize_t L = dsB->ncols();
            
            if( K == L)
            {
                std::vector<hsize_t> stride = {1, 1};
                std::vector<hsize_t> block = {1, 1};
                
                if (isSymmetric) {
                    if (N != M) {
                        throw std::range_error("Symmetric crossprod requires square result matrix");
                    }
                    //.. 20260304 ..// if (dsA->getFileName() != dsB->getFileName() || 
                    if (dsA->getFullPath() != dsB->getFullPath() || 
                        dsA->getGroup() != dsB->getGroup() || 
                        dsA->getDatasetName() != dsB->getDatasetName()) {
                        Rcpp::warning("isSymmetric=TRUE but different datasets provided. Results may be incorrect.");
                    }
                }
                
                dsC->inheritCompressionLevel(dsA->getCompressionLevel());
                dsC->createDataset( N, M, "real");

                // ── Preload strategy ─────────────────────────────────────────
                // Mirrors the approach in multiplication.hpp.
                // If A (and B) fit within ~20% of available RAM, read them
                // fully into memory and compute a single Eigen/BLAS multiply.
                // This avoids per-block decompression overhead and lets BLAS
                // manage its own thread parallelism over the full matrix.
                // For matrices that exceed the threshold, fall through to the
                // block-wise streaming path (PATH 2) below.
                const double mem_A_MB  = static_cast<double>(N * K) * 8.0 / (1024.0 * 1024.0);
                const double mem_B_MB  = static_cast<double>(M * L) * 8.0 / (1024.0 * 1024.0);
                const double avail_MB  = std::max(512.0, static_cast<double>(getAvailableMemoryMB()));
                const double thresh_MB = avail_MB * 0.20;

                const bool preload_A = (mem_A_MB <= thresh_MB);
                // For the symmetric case A == B, so preload_B follows preload_A.
                const bool preload_B = isSymmetric ? preload_A : (mem_B_MB <= thresh_MB);

                if (preload_A && preload_B) {
                    // ── PATH 1: Preload ───────────────────────────────────────
                    // Strategy: 1 HDF5 read per input + 1 BLAS multiply + blocked writes.
                    //
                    // HDF5/R coordinate convention:
                    //   A_R (nrows_R × ncols_R) is stored in HDF5 as [ncols_R × nrows_R].
                    //   Here N = dsA->nrows() = ncols_R, K = dsA->ncols() = nrows_R.
                    //   Reading {N, K} gives A_eigen (N×K, RowMajor) = t(A_R) mathematically.
                    //
                    // crossprod(A_R) = t(A_R) * A_R = A_eigen * A_eigen^T  (N×N symmetric)
                    // crossprod(A_R, B_R) = t(A_R) * B_R = A_eigen * B_eigen^T  (N×M general)
                    std::vector<hsize_t> stride = {1, 1};
                    std::vector<hsize_t> block  = {1, 1};

                    std::vector<double> vdA_full(N * K);
                    dsA->readDatasetBlock({0, 0}, {N, K}, stride, block, vdA_full.data());
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic,
                                            Eigen::Dynamic, Eigen::RowMajor>>
                        A_full(vdA_full.data(), static_cast<int>(N), static_cast<int>(K));

                    Eigen::MatrixXd C_full;
                    if (isSymmetric) {
                        // t(A_R) * A_R = A_eigen * A_eigen^T  (N×N)
                        C_full = A_full * A_full.transpose();
                    } else {
                        std::vector<double> vdB_full(M * L);
                        dsB->readDatasetBlock({0, 0}, {M, L}, stride, block, vdB_full.data());
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic,
                                                Eigen::Dynamic, Eigen::RowMajor>>
                            B_full(vdB_full.data(), static_cast<int>(M), static_cast<int>(L));
                        // t(A_R) * B_R = A_eigen * B_eigen^T  (N×M)
                        C_full = A_full * B_full.transpose();
                    }

                    // Write C to HDF5 in blocks aligned to hdf5_block.
                    // Same offset convention as PATH 2: offset = {jj, ii}, count = {mj, ni}.
                    for (hsize_t ii = 0; ii < N; ii += hdf5_block) {
                        hsize_t ni = std::min(hdf5_block, N - ii);
                        for (hsize_t jj = 0; jj < M; jj += hdf5_block) {
                            hsize_t mj = std::min(hdf5_block, M - jj);
                            Eigen::MatrixXd C_blk = C_full.block(
                                static_cast<int>(ii), static_cast<int>(jj),
                                static_cast<int>(ni), static_cast<int>(mj));
                            dsC->writeDatasetBlock(Rcpp::wrap(C_blk),
                                                   {jj, ii}, {mj, ni},
                                                   stride, block, false);
                        }
                    }

                } else {
                    // ── PATH 2: Block-wise streaming ──────────────────────────
                    // Matrices do not fit in the preload threshold.
                    // HDF5 I/O is serialised with #pragma omp critical(accessFile)
                    // to prevent concurrent file access; Eigen compute runs fully
                    // parallel between I/O calls.

#ifdef _OPENMP // Configure parallel processing
                int num_threads = 1;
                if (bparal) {
                    num_threads = get_number_threads(threads, Rcpp::wrap(bparal));  
                    omp_set_num_threads(num_threads);
                }
#endif                
                
                // Calculate total blocks for parallelization
                hsize_t blocks_i = (N + hdf5_block - 1) / hdf5_block;
                hsize_t blocks_j = (M + hdf5_block - 1) / hdf5_block;
                hsize_t total_blocks;
                
                if (isSymmetric) {
                    // For symmetric case: only upper triangle blocks
                    total_blocks = (blocks_i * (blocks_i + 1)) / 2;
                } else {
                    // For general case: all blocks
                    total_blocks = blocks_i * blocks_j;
                }
                
#ifdef _OPENMP
#pragma omp parallel for if(bparal) schedule(dynamic)
#endif
                for (hsize_t block_idx = 0; block_idx < total_blocks; ++block_idx)
                {
                    // Convert linear index to (ii_idx, jj_idx)
                    hsize_t ii_idx = 0, 
                            jj_idx = 0;
                    
                    if (isSymmetric) {
                        // Convert to upper triangle coordinates
                        hsize_t remaining = block_idx;
                        for (hsize_t i = 0; i < blocks_i; ++i) {
                            hsize_t blocks_in_row = blocks_i - i;
                            if (remaining < blocks_in_row) {
                                ii_idx = i;
                                jj_idx = i + remaining;
                                break;
                            }
                            remaining -= blocks_in_row;
                        }
                    } else {
                        // Convert to regular coordinates
                        ii_idx = block_idx / blocks_j;
                        jj_idx = block_idx % blocks_j;
                    }
                    
                    // Convert block indices to matrix indices
                    hsize_t ii = ii_idx * hdf5_block;
                    hsize_t jj = jj_idx * hdf5_block;
                    
                    // Boundary checks
                    if (ii >= N || jj >= M) continue;
                    
                    // Thread-local variables (each thread needs its own)
                    hsize_t local_isize = hdf5_block + 1;
                    hsize_t local_jsize = hdf5_block + 1;
                    hsize_t local_ksize = hdf5_block + 1;
                    
                    hsize_t iRowsA = getOptimBlockSize( N, hdf5_block, ii, local_isize);
                    if( ii + hdf5_block > N ) local_isize = N - ii;
                    
                    hsize_t iRowsB = getOptimBlockSize( M, hdf5_block, jj, local_jsize);
                    if( jj + hdf5_block > M) local_jsize = M - jj;
                    
                    Eigen::MatrixXd C_accum = Eigen::MatrixXd::Zero(iRowsA, iRowsB); 
                    
                    std::vector<hsize_t> offset = {jj, ii};        
                    std::vector<hsize_t> count = {iRowsB, iRowsA}; 
                    
                    for(hsize_t kk = 0; kk < K; kk += hdf5_block)
                    {
                        if( kk + hdf5_block > K ) local_ksize = K - kk;
                        
                        hsize_t iColsA = getOptimBlockSize( K, hdf5_block, kk, local_ksize),
                            iColsB = iColsA;
                        
                        // Pre-allocate data vectors outside critical section
                        std::vector<double> vdA(iRowsA * iColsA);
                        std::vector<double> vdB(iRowsB * iColsB);
                        
                        // // HDF5 read operations (inside critical section)
                        // dsA->readDatasetBlock( {ii, kk}, {iRowsA,iColsA}, stride, block, vdA.data() );
                        // dsB->readDatasetBlock( {jj, kk}, {iRowsB, iColsB}, stride, block, vdB.data() );
                        
                        // Serialise HDF5 file I/O across threads.
                        // Concurrent reads of the same B[jj,kk] block by threads sharing
                        // the same jj value cause race conditions in HDF5 without
                        // HDF5_ENABLE_THREADSAFE. The Eigen computation below this block
                        // remains fully parallel (operates on private thread-local data).
                        #ifdef _OPENMP
                        #pragma omp critical(accessFile)
                        #endif
                        {
                            dsA->readDatasetBlock( {ii, kk}, {iRowsA,iColsA}, stride, block, vdA.data() );
                            dsB->readDatasetBlock( {jj, kk}, {iRowsB, iColsB}, stride, block, vdB.data() );
                        }
                        
                        // Parallel computation (outside critical section)
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), iRowsA, iColsA );
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_B (vdB.data(), iRowsB, iColsB);   
                        Eigen::MatrixXd B = tmp_B.transpose();
                        
                        C_accum.noalias() += A * B;
                        
                        // Readjust counters
                        if( kk + hdf5_block > K ) local_ksize = hdf5_block + 1;
                        if( iColsA > hdf5_block ) {
                            kk = kk - hdf5_block + iColsA; 
                        }
                    }
                    
                    // dsC->writeDatasetBlock(Rcpp::wrap(C_accum), offset, count, stride, block, false);
                    // 
                    // if (isSymmetric && ii != jj) {
                    //     std::vector<hsize_t> offset_sym = {ii, jj};
                    //     std::vector<hsize_t> count_sym = {iRowsA, iRowsB};
                    //     dsC->writeDatasetBlock(Rcpp::wrap(C_accum.transpose()), offset_sym, count_sym, stride, block, false);
                    // }
                    
                    #ifdef _OPENMP
                    #pragma omp critical(accessFile)
                    #endif
                    {
                        dsC->writeDatasetBlock(Rcpp::wrap(C_accum), offset, count, stride, block, false);
                        if (isSymmetric && ii != jj) {
                            std::vector<hsize_t> offset_sym = {ii, jj};
                            std::vector<hsize_t> count_sym = {iRowsA, iRowsB};
                            dsC->writeDatasetBlock(Rcpp::wrap(C_accum.transpose()), offset_sym, count_sym, stride, block, false);
                        }
                    }
                    
                }
                
                } // end PATH 2 block-wise streaming (else of preload)
                
            } else {
                throw std::range_error("non-conformable arguments");
            }
            
        } catch(std::exception& ex) {
            throw std::runtime_error(std::string("c++ exception crossprod: ") + ex.what());
        }
        
        return(dsC);
    }

}

#endif // BIGDATASTATMETH_ALGEBRA_CROSSPROD_HPP
