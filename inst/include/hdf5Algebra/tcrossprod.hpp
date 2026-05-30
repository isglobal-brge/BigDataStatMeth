/**
 * @file tcrossprod.hpp
 * @brief Transposed cross-product operations for HDF5 matrices
 * @note 2026-03-07 Output datasets now inherit compression level from input datasets
 *         via setCompressionLevel() called before every createDataset() invocation.
 * @details This header file provides implementations for transposed cross-product
 * operations on matrices stored in HDF5 format. The implementation includes:
 * 
 * Key features:
 * - Transposed matrix multiplication
 * - Block-based computation
 * - Memory-efficient algorithms
 * - Parallel processing support
 * - Error handling and validation
 * 
 * Supported operations:
 * - A'A computation
 * - AA' computation
 * - Block cross-products
 * - Weighted cross-products
 * - Symmetric results
 * 
 * Performance features:
 * - Cache-friendly algorithms
 * - Dynamic block sizing
 * - Multi-threaded processing
 * - I/O optimization
 * - Memory management
 * 
 * The implementation uses:
 * - BLAS Level 3 operations
 * - Block algorithms
 * - HDF5 chunked storage
 * - Parallel I/O
 * - Symmetric optimizations
 */

#ifndef BIGDATASTATMETH_ALGEBRA_TCROSSPROD_HPP
#define BIGDATASTATMETH_ALGEBRA_TCROSSPROD_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

namespace BigDataStatMeth {
    
    /**
     * @brief Transposed cross-product for HDF5 matrices
     * @details Computes the transposed cross-product of matrices stored in HDF5 format.
     * Supports both A'A and AA' computations with block-based processing.
     * 
     * @param dsA Input matrix dataset
     * @param dsB Input matrix dataset
     * @param dsC Output matrix dataset
     * @param hdf5_block Block size for HDF5 I/O operations
     * @param mem_block_size Block size for in-memory operations
     * @param bparal Whether to use parallel processing
     * @param browmajor Whether matrices are stored in row-major order
     * @param threads Number of threads for parallel processing
     */
    inline BigDataStatMeth::hdf5Dataset* tcrossprod( 
            BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, 
            BigDataStatMeth::hdf5Dataset* dsC, bool isSymmetric, hsize_t hdf5_block, 
            hsize_t mem_block_size, bool bparal, bool browmajor, 
            Rcpp::Nullable<int> threads = R_NilValue) 

    {
        
        try {
            
            hsize_t N = dsA->ncols();  
            hsize_t K = dsA->nrows();  
            hsize_t M = dsB->ncols();  
            hsize_t L = dsB->nrows();  
            
            if (K != L) {
                throw std::range_error("non-conformable arguments");
            }
            
            if( K == L)
            {
                if (isSymmetric) {
                    if (N != M) {
                        throw std::range_error("Symmetric tcrossprod requires square result matrix");
                    }
                    //.. 20260304 ..// if (dsA->getFileName() != dsB->getFileName() || 
                    if (dsA->getFullPath() != dsB->getFullPath() || 
                        dsA->getGroup() != dsB->getGroup() || 
                        dsA->getDatasetName() != dsB->getDatasetName()) {
                        Rcpp::warning("isSymmetric=TRUE but different datasets provided. Results may be incorrect.");
                    }
                }

#ifdef _OPENMP  // Configure parallel processing
                int num_threads = 1;
                if (bparal) {
                    num_threads = get_number_threads(threads, Rcpp::wrap(bparal));
                    omp_set_num_threads(num_threads);
                }
#endif
                
                // Calculate total blocks for parallelization
                hsize_t blocks_i = (N + hdf5_block - 1) / hdf5_block;
                hsize_t blocks_j = isSymmetric ? blocks_i : (M + hdf5_block - 1) / hdf5_block;
                hsize_t total_blocks;
                
                if (isSymmetric) {
                    // For symmetric case: only upper triangle blocks
                    total_blocks = (blocks_i * (blocks_i + 1)) / 2;
                } else {
                    // For general case: all blocks
                    total_blocks = blocks_i * blocks_j;
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
                    // HDF5/R coordinate convention for tcrossprod:
                    //   A_R (nrows_R × ncols_R) stored in HDF5 as [ncols_R × nrows_R].
                    //   Here K = dsA->nrows() = ncols_R, N = dsA->ncols() = nrows_R.
                    //   Reading {K, N} gives A_eigen (K×N, RowMajor) = t(A_R) mathematically.
                    //   A_R = A_eigen^T  (N×K in math).
                    //
                    // tcrossprod(A_R) = A_R * t(A_R) = A_eigen^T * A_eigen  (N×N symmetric)
                    // tcrossprod(A_R, B_R) = A_R * t(B_R) = A_eigen^T * B_eigen  (N×M general)
                    std::vector<hsize_t> stride = {1, 1};
                    std::vector<hsize_t> block  = {1, 1};

                    std::vector<double> vdA_full(K * N);
                    dsA->readDatasetBlock({0, 0}, {K, N}, stride, block, vdA_full.data());
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic,
                                            Eigen::Dynamic, Eigen::RowMajor>>
                        A_full(vdA_full.data(), static_cast<int>(K), static_cast<int>(N));

                    Eigen::MatrixXd C_full;
                    if (isSymmetric) {
                        // A_R * t(A_R) = A_eigen^T * A_eigen  (N×N)
                        C_full = A_full.transpose() * A_full;
                    } else {
                        std::vector<double> vdB_full(L * M);
                        dsB->readDatasetBlock({0, 0}, {L, M}, stride, block, vdB_full.data());
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic,
                                                Eigen::Dynamic, Eigen::RowMajor>>
                            B_full(vdB_full.data(), static_cast<int>(L), static_cast<int>(M));
                        // A_R * t(B_R) = A_eigen^T * B_eigen  (N×M)
                        C_full = A_full.transpose() * B_full;
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
                            std::vector<hsize_t> stride_w = {1, 1};
                            std::vector<hsize_t> block_w  = {1, 1};
                            dsC->writeDatasetBlock(Rcpp::wrap(C_blk),
                                                   {jj, ii}, {mj, ni},
                                                   stride_w, block_w, false);
                        }
                    }

                } else {
                    // ── PATH 2: Block-wise streaming ──────────────────────────
                    // Matrices do not fit in the preload threshold.
                    // HDF5 I/O is serialised with #pragma omp critical(accessFile)
                    // to prevent concurrent file access; Eigen compute runs fully
                    // parallel between I/O calls.

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
                    
                    hsize_t Nii = std::min(hdf5_block, N - ii);
                    hsize_t Mjj = std::min(hdf5_block, M - jj);
                    
                    // Thread-local variables
                    std::vector<hsize_t> stride = {1, 1};
                    std::vector<hsize_t> block = {1, 1};
                    
                    Eigen::MatrixXd C_accum = Eigen::MatrixXd::Zero(Nii, Mjj);
                    
                    std::vector<hsize_t> offset = {jj, ii};
                    std::vector<hsize_t> count = {Mjj, Nii};
                    
                    for(hsize_t kk = 0; kk < K; kk += hdf5_block)
                    {
                        hsize_t Kkk = std::min(hdf5_block, K - kk); 
                        
                        Eigen::MatrixXd A;
                        
                        // {
                        //     std::vector<double> vdA( Nii * Kkk);
                        //     dsA->readDatasetBlock( {kk, ii}, {Kkk, Nii}, stride, block, vdA.data() );
                        //     Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_A (vdA.data(), Kkk, Nii );    
                        //     A = tmp_A.transpose();
                        // }
                        // 
                        // std::vector<double> vdB( Mjj * Kkk);
                        // dsB->readDatasetBlock( {kk, jj}, { Kkk, Mjj}, stride, block, vdB.data() );
                        
                        {
                            std::vector<double> vdA( Nii * Kkk);
                            // Serialise HDF5 reads: concurrent reads of the same
                            // B[kk,jj] block by threads sharing the same jj value
                            // are a race condition without HDF5_ENABLE_THREADSAFE.
                        #ifdef _OPENMP
                        #pragma omp critical(accessFile)
                        #endif
                        {
                            dsA->readDatasetBlock( {kk, ii}, {Kkk, Nii}, stride, block, vdA.data() );
                        }
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_A (vdA.data(), Kkk, Nii );    
                        A = tmp_A.transpose();
                                                }
                                                
                                                std::vector<double> vdB( Mjj * Kkk);
                        #ifdef _OPENMP
                        #pragma omp critical(accessFile)
                        #endif
                        {
                            dsB->readDatasetBlock( {kk, jj}, { Kkk, Mjj}, stride, block, vdB.data() );
                        }
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), Kkk, Mjj);   
                        
                        C_accum.noalias() += A * B;
                    }
                    
                    // dsC->writeDatasetBlock(Rcpp::wrap(C_accum), offset, count, stride, block, false);
                    // 
                    // // Symetric matrices
                    // if (isSymmetric && ii != jj) {
                    //     std::vector<hsize_t> offset_sym = {ii, jj};
                    //     std::vector<hsize_t> count_sym = {Nii, Mjj};
                    //     dsC->writeDatasetBlock(Rcpp::wrap(C_accum.transpose()), offset_sym, count_sym, stride, block, false);
                    // }
                
                    // Serialise writes: regions are non-overlapping per thread
                    // but HDF5 file state is shared across threads.
                    #ifdef _OPENMP
                    #pragma omp critical(accessFile)
                    #endif
                    {
                        dsC->writeDatasetBlock(Rcpp::wrap(C_accum), offset, count, stride, block, false);
                        // Symmetric matrices: write the transposed block
                        if (isSymmetric && ii != jj) {
                            std::vector<hsize_t> offset_sym = {ii, jj};
                            std::vector<hsize_t> count_sym = {Nii, Mjj};
                            dsC->writeDatasetBlock(Rcpp::wrap(C_accum.transpose()), offset_sym, count_sym, stride, block, false);
                        }
                    }
                
                } // end PATH 2 block-wise streaming (else of preload)
                
                } // end if(K == L)
                
            } else {
                throw std::range_error("non-conformable arguments");
            }
            
        } catch(std::exception& ex) {
            throw std::runtime_error(std::string("c++ exception tcrossprod: ") + ex.what());
        }
        
        return(dsC);
    }

}

#endif // BIGDATASTATMETH_ALGEBRA_TCROSSPROD_HPP
