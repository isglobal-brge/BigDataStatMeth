/**
 * @file multiplication.hpp
 * @brief Matrix multiplication operations for HDF5 matrices
 * @note 2026-03-07 Output datasets now inherit compression level from input datasets
 *         via setCompressionLevel() called before every createDataset() invocation.
 * @details This header file provides implementations for matrix multiplication
 * operations on matrices stored in HDF5 format. The implementation includes:
 * 
 * Key features:
 * - Dense matrix multiplication
 * - Block-based multiplication
 * - Parallel processing support
 * - Memory-efficient algorithms
 * - Error handling and validation
 * 
 * Supported operations:
 * - Matrix-matrix multiplication
 * - Block matrix multiplication
 * - Transposed multiplication
 * - Multi-threaded multiplication
 * - Out-of-core processing
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
 * - Cache blocking
 */

#ifndef BIGDATASTATMETH_HDF5_MULTIPLICATION_HPP
#define BIGDATASTATMETH_HDF5_MULTIPLICATION_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"
#include "Utilities/system-utils.hpp"

namespace BigDataStatMeth {


/**
 * @brief Main matrix multiplication function for HDF5 matrices
 * @details Performs matrix multiplication C = A * B where A, B, and C are HDF5 datasets.
 * Supports parallel processing and block-based computation for memory efficiency.
 * 
 * @param dsA First input matrix dataset
 * @param dsB Second input matrix dataset
 * @param dsC Output matrix dataset
 * @param transpose_A Whether to transpose matrix A
 * @param transpose_B Whether to transpose matrix B
 * @param bparal Whether to use parallel processing
 * @param hdf5_block Block size for HDF5 I/O operations
 * @param threads Number of threads for parallel processing
 */
// inline void multiplication( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
//                                    Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> hdf5_block, Rcpp::Nullable<int> threads);

inline void multiplication( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
                            bool transpose_A, bool transpose_B, Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> hdf5_block, Rcpp::Nullable<int> threads);


/**
 * @brief Calculate block positions and sizes for HDF5 matrix operations
 * @details Determines optimal block positions and sizes for block-based matrix
 * operations on HDF5 datasets.
 * 
 * @param maxPosition Maximum position to process
 * @param blockSize Size of each block
 * @param[out] starts Vector to store starting positions of blocks
 * @param[out] sizes Vector to store sizes of blocks
 */
inline void getBlockPositionsSizes_hdf5( hsize_t maxPosition, hsize_t blockSize, std::vector<hsize_t>& starts, std::vector<hsize_t>& sizes ){

        hsize_t isize = blockSize + 1;

        for (hsize_t ii = 0; ii < maxPosition; ii += blockSize)
        {
            if( ii + blockSize > maxPosition ) {
                isize = maxPosition - ii; }

            hsize_t sizetoRead = getOptimBlockSize( maxPosition, blockSize, ii, isize);

            starts.push_back(ii);
            sizes.push_back(sizetoRead);

            // if( ii + blockSize > maxPosition ) {
            //     isize = blockSize + 1; }
            if( sizetoRead > blockSize ) {
                ii = ii - blockSize + sizetoRead; }
        }

    }


    /**
     * @brief Calculate optimal block sizes specifically for multiplication
     * @param N Number of result rows (dsA->ncols)
     * @param M Number of result cols (dsB->nrows)  
     * @param K Inner dimension (dsA->nrows = dsB->ncols)
     */
    inline BlockSizes calculate_multiplication_blocks(hsize_t N, hsize_t M, hsize_t K) {
        const hsize_t MEMORY_BUDGET = (hsize_t)(1.5 * 1024 * 1024 * 1024);  // 1.5GB
        const hsize_t bytes_per_element = sizeof(double);
        const hsize_t min_block = 256;
        
        BlockSizes result;
        result.type = classify_matrix_type(N, M, K);
        
        switch (result.type) {
        case RECTANGULAR_EXTREME: {
            // Para big-omics: K >> N,M - priorizar bloques grandes en K
            hsize_t small_dim = std::min(N, M);
            result.output_block = std::min(small_dim, (hsize_t)1024);
            
            // Maximizar bloque K con memoria restante
            // Memory = A_block + B_block = K×(N+M) aprox
            hsize_t remaining_budget = (hsize_t)(0.8 * MEMORY_BUDGET);
            result.inner_block = remaining_budget / ((N + M) * bytes_per_element);
            result.inner_block = std::min({result.inner_block, K, (hsize_t)32768});
            break;
        }
            
        case SQUARE_SMALL: {
            hsize_t block = sqrt(MEMORY_BUDGET / (3 * bytes_per_element));
            // block = std::min({block, N/2, M/2, K/2});
            block = std::min({block, N, M, K});
            result.inner_block = result.output_block = block;
            break;
        }
            
        case SQUARE_LARGE: {
            hsize_t target_blocks = 2;  // 2×2 output blocks
            result.output_block = std::max({N / target_blocks, M / target_blocks, min_block});
            
            hsize_t output_memory = result.output_block * result.output_block * bytes_per_element;
            hsize_t remaining_budget = MEMORY_BUDGET - output_memory;
            result.inner_block = remaining_budget / (2 * result.output_block * bytes_per_element);
            result.inner_block = std::min({result.inner_block, K, (hsize_t)16384});
            break;
        }
            
        case SQUARE_EXTREME: {
            hsize_t block = sqrt(MEMORY_BUDGET / (3 * bytes_per_element));
            block = std::min({block, (hsize_t)8192});
            result.inner_block = result.output_block = block;
            break;
        }
        }
        
        // Aplicar límites y redondear
        result.inner_block = std::max(result.inner_block, min_block);
        result.output_block = std::max(result.output_block, min_block);
        
        result.inner_block = ((result.inner_block + 63) / 64) * 64;
        result.output_block = ((result.output_block + 63) / 64) * 64;
        
        return result;
    }

    // In-memory execution - Parallel version
    // 
    //  IMPORTANT : FUNCIÓ MODIFICADA EL 2024/04/06  I NO TESTEJADA !!!!
    // 

    /**
     * @brief Parallel block-based matrix multiplication
     * @details Implements parallel block-based matrix multiplication for in-memory matrices.
     * 
     * @param A First input matrix
     * @param B Second input matrix
     * @param block_size Size of blocks for computation
     * @param threads Number of threads for parallel processing
     * @return Result of matrix multiplication
     */
    inline Eigen::MatrixXd Bblock_matrix_mul_parallel( Eigen::MatrixXd A, Eigen::MatrixXd B, 
                                                             int block_size, Rcpp::Nullable<int> threads  = R_NilValue)
    {
        
        // unsigned int ithreads;
        Eigen::MatrixXd C;
        
        try {

            int M = A.rows();
            int K = A.cols();
            int N = B.cols();
            
            C = Eigen::MatrixXd::Zero(M,N) ;
            
            if( A.cols() == B.rows())   // inner dimensions must match for A*B
            {
                
                if(block_size > std::min( N, std::min(M,K)) )
                    block_size = std::min( N, std::min(M,K)); 
                
                std::vector<hsize_t> vsizetoRead, vstart,
                                     vsizetoReadM, vstartM,
                                     vsizetoReadK, vstartK;
                
                getBlockPositionsSizes_hdf5( N, block_size, vstart, vsizetoRead );
                getBlockPositionsSizes_hdf5( M, block_size, vstartM, vsizetoReadM );
                getBlockPositionsSizes_hdf5( K, block_size, vstartK, vsizetoReadK );
                
                // int ithreads = get_number_threads(threads, R_NilValue);
                // int chunks = vstart.size()/ithreads;
                
                #pragma omp parallel num_threads( get_number_threads(threads, R_NilValue) ) shared(A, B, C) //..// , chunk) private(tid ) 
                {
                    
                    #pragma omp for schedule (static) // collapse(3)
                    for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                    {
                        for (hsize_t jj = 0; jj < vstartM.size(); jj++)
                        {
                            for (hsize_t kk = 0; kk < vstartK.size(); kk++)
                            {
                                C.block(vstart[ii], vstartM[jj], vsizetoRead[ii], vsizetoReadM[jj]) = 
                                    C.block(vstart[ii], vstartM[jj], vsizetoRead[ii], vsizetoReadM[jj]) + 
                                    ( A.block(vstart[ii], vstartK[kk], vsizetoRead[ii], vsizetoReadK[kk]) * 
                                      B.block(vstartK[kk], vstartM[jj], vsizetoReadK[kk], vsizetoReadM[jj]) );
                            }
                        }
                    }
                }
                
            } else {
                throw std::range_error("multiplication error: non-conformable arguments");
            }
            
        } catch(std::exception& ex) {
            throw std::runtime_error(std::string("c++ exception Bblock_matrix_mul_parallel: ") + ex.what());
        }
        
        return(C);
        
    }
    
    
    inline void multiplication( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
                                bool transpose_A, bool transpose_B, Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> hdf5_block, Rcpp::Nullable<int> threads = R_NilValue) 
    {
        
        try {
            
            hsize_t K, N, L, M;
            
            if (transpose_A) {
                K = dsA->ncols();  // inner dim for t(A)*B = R's nrows(A) [HDF5 second dim]
                N = dsA->nrows();  // output rows of t(A)*B = R's ncols(A) [HDF5 first dim]
            } else {
                K = dsA->nrows();  // inner dim for A*B = R's ncols(A) [HDF5 first dim]
                N = dsA->ncols();  // output rows of A*B = R's nrows(A) [HDF5 second dim]
            }
            
            if (transpose_B) {
                L = dsB->nrows();  // inner dim check: R's ncols(B) [HDF5 first dim]
                M = dsB->ncols();  // output cols of A*t(B) = R's nrows(B) [HDF5 second dim]
            } else {
                L = dsB->ncols();  // inner dim check: R's nrows(B) [HDF5 second dim]
                M = dsB->nrows();  // output cols of A*B = R's ncols(B) [HDF5 first dim]
            }
            
            int ihdf5_block_N, ihdf5_block_M, ihdf5_block_K;
            
            if( hdf5_block.isNotNull()) {
                ihdf5_block_N = ihdf5_block_M = ihdf5_block_K = Rcpp::as<int>(hdf5_block);
            } else {
                BlockSizes blocks = calculate_multiplication_blocks(N, M, K);
                ihdf5_block_N = ihdf5_block_M = blocks.output_block;
                ihdf5_block_K = blocks.inner_block;
            }
            
            if( K == L )
            {
                std::vector<hsize_t> stride = {1, 1},
                    block  = {1, 1},
                    vsizetoRead, vstart,
                    vsizetoReadM, vstartM,
                    vsizetoReadK, vstartK;
                
                dsC->inheritCompressionLevel(dsA->getCompressionLevel());
                dsC->createDataset( N, M, "real");
                
                if( dsC->getDatasetptr() != nullptr) 
                {
                    getBlockPositionsSizes_hdf5( N, ihdf5_block_N, vstart,  vsizetoRead  );
                    getBlockPositionsSizes_hdf5( M, ihdf5_block_M, vstartM, vsizetoReadM );
                    getBlockPositionsSizes_hdf5( K, ihdf5_block_K, vstartK, vsizetoReadK );
                    
                    // --- Preload strategy ---
                    // Estimate memory footprint of A and B (in MB).
                    // If a matrix fits within 20% of available RAM it is read from HDF5 once
                    // and reused entirely from RAM, eliminating redundant disk reads.
                    // Threshold is conservative (20%) to leave room for C accumulators and
                    // OS/HDF5 caches.
                    const double mem_A_MB  = static_cast<double>(N * K) * 8.0 / (1024.0 * 1024.0);
                    const double mem_B_MB  = static_cast<double>(M * K) * 8.0 / (1024.0 * 1024.0);
                    const double avail_MB  = std::max(512.0, static_cast<double>(getAvailableMemoryMB()));
                    const double thresh_MB = avail_MB * 0.20;
                    
                    const bool preload_A = (mem_A_MB <= thresh_MB);
                    const bool preload_B = (mem_B_MB <= thresh_MB);
                    
                    // ═══════════════════════════════════════════════════════════════════
                    // PATH 1 — both A and B fit in RAM
                    //
                    // Strategy: 2 HDF5 reads (full A, full B) + 1 BLAS multiply + block writes.
                    // Zero disk reads inside the compute loop.
                    //
                    // HDF5 physical layout (BigDataStatMeth convention, transposed vs R):
                    //   A stored as [N × K] when transpose_A, else [K × N]
                    //   B stored as [K × M] when transpose_B, else [M × K]
                    //   C stored as [M × N]  (createDataset(N,M) → HDF5 [M,N])
                    // ═══════════════════════════════════════════════════════════════════
                    if (preload_A && preload_B)
                    {
                        // Read full A into RAM
                        std::vector<double> vdA_full(N * K);
                        if (transpose_A)
                            dsA->readDatasetBlock({0, 0}, {N, K}, stride, block, vdA_full.data());
                        else
                            dsA->readDatasetBlock({0, 0}, {K, N}, stride, block, vdA_full.data());
                        Eigen::MatrixXd A_full = Eigen::Map<
                            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                                    vdA_full.data(),
                                    transpose_A ? static_cast<int>(N) : static_cast<int>(K),
                                    transpose_A ? static_cast<int>(K) : static_cast<int>(N));
                        
                        // Read full B into RAM
                        std::vector<double> vdB_full(M * K);
                        if (transpose_B)
                            dsB->readDatasetBlock({0, 0}, {K, M}, stride, block, vdB_full.data());
                        else
                            dsB->readDatasetBlock({0, 0}, {M, K}, stride, block, vdB_full.data());
                        Eigen::MatrixXd B_full = Eigen::Map<
                            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                                    vdB_full.data(),
                                    transpose_B ? static_cast<int>(K) : static_cast<int>(M),
                                    transpose_B ? static_cast<int>(M) : static_cast<int>(K));
                        
                        // Single BLAS-optimised multiply — C is M×N
                        // Operación según flags de transposición
                        Eigen::MatrixXd C_full;
                        if (!transpose_A && !transpose_B)      C_full = B_full * A_full;                           // A * B
                        else if (transpose_A && !transpose_B)  C_full = B_full * A_full.transpose();               // t(A) * B
                        else if (!transpose_A && transpose_B)  C_full = B_full.transpose() * A_full;               // A * t(B)
                        else                                   C_full = B_full.transpose() * A_full.transpose();   // t(A) * t(B)
                        
                        // Write C block-by-block to HDF5 (C_full is M×N)
                        for (hsize_t ii = 0; ii < vstart.size(); ii++) {
                            for (hsize_t jj = 0; jj < vstartM.size(); jj++) {
                                std::vector<double> vdC_final(vsizetoReadM[jj] * vsizetoRead[ii]);
                                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
                                    C_final_map(vdC_final.data(), vsizetoReadM[jj], vsizetoRead[ii]);
                                C_final_map = C_full.block(vstartM[jj], vstart[ii], vsizetoReadM[jj], vsizetoRead[ii]);
                                
                                std::vector<hsize_t> offset = {vstartM[jj], vstart[ii]};
                                std::vector<hsize_t> count  = {vsizetoReadM[jj], vsizetoRead[ii]};
                                dsC->writeDatasetBlock(vdC_final, offset, count, stride, block);
                            }
                        }
                        
                    } else {
                        
                        // ═══════════════════════════════════════════════════════════════════
                        // PATH 2 / 3 — streaming, with optional partial preload
                        //
                        // Loop order: parallel(ii) → kk → jj
                        //
                        // Why kk before jj (vs old order jj before kk):
                        //   Old order  parallel(ii) → jj → kk:
                        //     · A[ii,kk] read M_blocks × K_blocks times per ii  (M_blocks redundant)
                        //     · B[jj,kk] read K_blocks times per (ii,jj)
                        //   New order  parallel(ii) → kk → jj:
                        //     · A[ii,kk] read K_blocks times per ii  (once per kk, reused for all jj)
                        //     · B[jj,kk] read from RAM (preload_B) or disk once per (ii,jj,kk)
                        //     · C_acc[jj] accumulated across all kk; single write per (ii,jj)
                        //
                        // Preload_A: A_full_mat in RAM → zero A reads in the loop.
                        // Preload_B: B_full_mat in RAM → zero B reads in the loop.
                        // Neither:   A read once per (ii,kk); B still read once per (ii,jj,kk).
                        // ═══════════════════════════════════════════════════════════════════
                        
                        // Preload whichever matrix fits in RAM (declared before parallel region)
                        Eigen::MatrixXd A_full_mat, B_full_mat;
                        
                        if (preload_A) {
                            // A_full_mat layout: (N×K) if transpose_A, else (K×N)
                            std::vector<double> vdA_full(N * K);
                            if (transpose_A)
                                dsA->readDatasetBlock({0, 0}, {N, K}, stride, block, vdA_full.data());
                            else
                                dsA->readDatasetBlock({0, 0}, {K, N}, stride, block, vdA_full.data());
                            A_full_mat = Eigen::Map<
                                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                                        vdA_full.data(),
                                        transpose_A ? static_cast<int>(N) : static_cast<int>(K),
                                        transpose_A ? static_cast<int>(K) : static_cast<int>(N));
                        }
                        
                        if (preload_B) {
                            // B_full_mat layout: (K×M) if transpose_B, else (M×K)
                            std::vector<double> vdB_full(M * K);
                            if (transpose_B)
                                dsB->readDatasetBlock({0, 0}, {K, M}, stride, block, vdB_full.data());
                            else
                                dsB->readDatasetBlock({0, 0}, {M, K}, stride, block, vdB_full.data());
                            B_full_mat = Eigen::Map<
                                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                                        vdB_full.data(),
                                        transpose_B ? static_cast<int>(K) : static_cast<int>(M),
                                        transpose_B ? static_cast<int>(M) : static_cast<int>(K));
                        }
                        
#pragma omp parallel num_threads( get_number_threads(threads, R_NilValue) ) shared(dsA, dsB, dsC, A_full_mat, B_full_mat, vstart, vsizetoRead, vstartM, vsizetoReadM, vstartK, vsizetoReadK)
{
#pragma omp for schedule(dynamic) nowait
    for (hsize_t ii = 0; ii < vstart.size(); ii++)
    {
        // One C accumulator per output-column block, zeroed once per ii.
        // Accumulates across all kk before a single write — avoids
        // re-reading C from disk.
        std::vector<Eigen::MatrixXd> C_acc(vstartM.size());
        for (hsize_t jj = 0; jj < vstartM.size(); jj++)
            C_acc[jj] = Eigen::MatrixXd::Zero(vsizetoReadM[jj], vsizetoRead[ii]);
        
        for (hsize_t kk = 0; kk < vstartK.size(); kk++)
        {
            // --- Read or fetch block of A (once per (ii,kk), reused for all jj) ---
            hsize_t rowsA = transpose_A ? vsizetoRead[ii]  : vsizetoReadK[kk];
            hsize_t colsA = transpose_A ? vsizetoReadK[kk] : vsizetoRead[ii];
            Eigen::MatrixXd A_block;
            
            if (preload_A) {
                // No I/O: extract sub-block from preloaded A_full_mat
                // A_full_mat is (N×K) if transpose_A, else (K×N)
                if (transpose_A)
                    A_block = A_full_mat.block(vstart[ii],   vstartK[kk], vsizetoRead[ii],  vsizetoReadK[kk]);
                else
                    A_block = A_full_mat.block(vstartK[kk],  vstart[ii],  vsizetoReadK[kk], vsizetoRead[ii]);
            } else {
                // Read A block from HDF5
                //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
                std::vector<double> vdA(rowsA * colsA);
                if (transpose_A) {
                    // HDF5 dim1=N, dim2=K  →  read {N_offset, K_offset}
                    dsA->readDatasetBlock( {vstart[ii], vstartK[kk]}, {vsizetoRead[ii], vsizetoReadK[kk]}, stride, block, vdA.data() );
                } else {
                    // HDF5 dim1=K, dim2=N  →  read {K_offset, N_offset}
                    dsA->readDatasetBlock( {vstartK[kk], vstart[ii]}, {vsizetoReadK[kk], vsizetoRead[ii]}, stride, block, vdA.data() );
                }
                A_block = Eigen::Map<
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                            vdA.data(), rowsA, colsA);
            }
            
            // --- Loop over output-column blocks, read or fetch B ---
            for (hsize_t jj = 0; jj < vstartM.size(); jj++)
            {
                hsize_t rowsB = transpose_B ? vsizetoReadK[kk] : vsizetoReadM[jj];
                hsize_t colsB = transpose_B ? vsizetoReadM[jj] : vsizetoReadK[kk];
                Eigen::MatrixXd B_block;
                
                if (preload_B) {
                    // No I/O: extract sub-block from preloaded B_full_mat
                    // B_full_mat is (K×M) if transpose_B, else (M×K)
                    if (transpose_B)
                        B_block = B_full_mat.block(vstartK[kk],  vstartM[jj], vsizetoReadK[kk], vsizetoReadM[jj]);
                    else
                        B_block = B_full_mat.block(vstartM[jj],  vstartK[kk], vsizetoReadM[jj], vsizetoReadK[kk]);
                } else {
                    // Read B block from HDF5
                    //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
                    std::vector<double> vdB(rowsB * colsB);
                    if (transpose_B) {
                        // HDF5 dim1=K, dim2=M  →  read {K_offset, M_offset}
                        dsB->readDatasetBlock( {vstartK[kk], vstartM[jj]}, {vsizetoReadK[kk], vsizetoReadM[jj]}, stride, block, vdB.data() );
                    } else {
                        // HDF5 dim1=M, dim2=K  →  read {M_offset, K_offset}
                        dsB->readDatasetBlock( {vstartM[jj], vstartK[kk]}, {vsizetoReadM[jj], vsizetoReadK[kk]}, stride, block, vdB.data() );
                    }
                    B_block = Eigen::Map<
                        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                                vdB.data(), rowsB, colsB);
                }
                
                // C_accumulator += B * A;
                // Operación según flags de transposición
                if (!transpose_A && !transpose_B) {
                    C_acc[jj] += B_block * A_block;                           // A * B
                } else if (transpose_A && !transpose_B) {
                    C_acc[jj] += B_block * A_block.transpose();               // t(A) * B
                } else if (!transpose_A && transpose_B) {
                    C_acc[jj] += B_block.transpose() * A_block;               // A * t(B)
                } else {
                    C_acc[jj] += B_block.transpose() * A_block.transpose();   // t(A) * t(B)
                }
            }
        }
        
        // Write all output-column blocks for this ii row-block (once per ii,
        // after accumulating all kk — avoids re-reading C)
        for (hsize_t jj = 0; jj < vstartM.size(); jj++) {
            std::vector<double> vdC_final(vsizetoReadM[jj] * vsizetoRead[ii]);
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
                C_final_map(vdC_final.data(), vsizetoReadM[jj], vsizetoRead[ii]);
            C_final_map = C_acc[jj];
            
            std::vector<hsize_t> offset = {vstartM[jj], vstart[ii]};
            std::vector<hsize_t> count  = {vsizetoReadM[jj], vsizetoRead[ii]};
            //.. 20260325 - remove critical ..// #pragma omp critical(accessFile)
            dsC->writeDatasetBlock(vdC_final, offset, count, stride, block);
        }
    }
}
                    }
                }
                
            } else {
                throw std::range_error("multiplication error: non-conformable arguments");
            }
            
        }  catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            // checkClose_file(dsA, dsB, dsC);
            throw std::runtime_error("c++ c++ exception multiplication (File IException)");
            // return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            // checkClose_file(dsA, dsB, dsC);
            throw std::runtime_error("c++ exception multiplication (DataSet IException)");
            // return void();
        } catch(std::exception &ex) {
            // checkClose_file(dsA, dsB, dsC);
            throw std::runtime_error(std::string("c++ exception multiplication: ") + ex.what());
            // return void();
        }  catch (...) {
            // checkClose_file(dsA, dsB, dsC);
            throw std::runtime_error("C++ exception multiplication (unknown reason)");
            // return void();
        }
        
        return void();
    }

    // inline void multiplication( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
    //                             bool transpose_A, bool transpose_B, Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> hdf5_block, Rcpp::Nullable<int> threads = R_NilValue) 
    // {
    //     
    //     try {
    //               
    //          hsize_t K, N, L, M;
    //         
    //         if (transpose_A) {
    //             K = dsA->ncols();  // inner dim for t(A)*B = R's nrows(A) [HDF5 second dim]
    //             N = dsA->nrows();  // output rows of t(A)*B = R's ncols(A) [HDF5 first dim]
    //         } else {
    //             K = dsA->nrows();  // inner dim for A*B = R's ncols(A) [HDF5 first dim]
    //             N = dsA->ncols();  // output rows of A*B = R's nrows(A) [HDF5 second dim]
    //         }
    //         
    //         if (transpose_B) {
    //             L = dsB->nrows();  // inner dim check: R's ncols(B) [HDF5 first dim]
    //             M = dsB->ncols();  // output cols of A*t(B) = R's nrows(B) [HDF5 second dim]
    //         } else {
    //             L = dsB->ncols();  // inner dim check: R's nrows(B) [HDF5 second dim]
    //             M = dsB->nrows();  // output cols of A*B = R's ncols(B) [HDF5 first dim]
    //         }
    //          
    //          int ihdf5_block_N, ihdf5_block_M, ihdf5_block_K;
    //          
    //          if( hdf5_block.isNotNull()) {
    //              ihdf5_block_N = ihdf5_block_M = ihdf5_block_K = Rcpp::as<int>(hdf5_block);
    //          } else {
    //              BlockSizes blocks = calculate_multiplication_blocks(N, M, K);
    //              ihdf5_block_N = ihdf5_block_M = blocks.output_block;
    //              ihdf5_block_K = blocks.inner_block;
    //          }
    //          
    //         if( K == L )
    //         {
    //             std::vector<hsize_t> stride = {1, 1},
    //                                  block = {1, 1},
    //                                  vsizetoRead, vstart,
    //                                  vsizetoReadM, vstartM,
    //                                  vsizetoReadK, vstartK;
    //             
    //             dsC->inheritCompressionLevel(dsA->getCompressionLevel());
    //             dsC->createDataset( N, M, "real");
    //             
    //             if( dsC->getDatasetptr() != nullptr) 
    //             {
    //                 getBlockPositionsSizes_hdf5( N, ihdf5_block_N, vstart, vsizetoRead );
    //                 getBlockPositionsSizes_hdf5( M, ihdf5_block_M, vstartM, vsizetoReadM );
    //                 getBlockPositionsSizes_hdf5( K, ihdf5_block_K, vstartK, vsizetoReadK );
    //                 
    //                 #pragma omp parallel num_threads( get_number_threads(threads, R_NilValue) ) shared(dsA, dsB, dsC, vstart, vsizetoRead) // chunks
    //                 {
    //                     
    //                     #pragma omp for schedule(dynamic) nowait
    //                     for (hsize_t ii = 0; ii < vstart.size(); ii++)
    //                     {
    //                         for (hsize_t jj = 0; jj < vstartM.size(); jj++)
    //                         {
    //                             
    //                             Eigen::MatrixXd C_accumulator = Eigen::MatrixXd::Zero(vsizetoReadM[jj], vsizetoRead[ii]);
    //                             
    //                             for (hsize_t kk = 0; kk < vstartK.size(); kk++)
    //                             {
    //                                 hsize_t iColsA = vsizetoReadK[kk],
    //                                         iRowsA = vsizetoRead[ii],
    //                                         iColsB = vsizetoReadM[jj],
    //                                         iRowsB = vsizetoReadK[kk];
    //                                 
    //                                 std::vector<double> vdA( iRowsA * iColsA );
    //                                 if (transpose_A) {
    //                                     // HDF5 dim1=N, dim2=K  →  read {N_offset, K_offset}
    //                                     dsA->readDatasetBlock( {vstart[ii], vstartK[kk]}, {vsizetoRead[ii], vsizetoReadK[kk]}, stride, block, vdA.data() );
    //                                 } else {
    //                                     // HDF5 dim1=K, dim2=N  →  read {K_offset, N_offset}
    //                                     dsA->readDatasetBlock( {vstartK[kk], vstart[ii]}, {vsizetoReadK[kk], vsizetoRead[ii]}, stride, block, vdA.data() );
    //                                 }
    //                                 Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(),
    //                                      transpose_A ? vsizetoRead[ii]    : vsizetoReadK[kk],   // rows del bloque en memoria
    //                                      transpose_A ? vsizetoReadK[kk]   : vsizetoRead[ii]);   // cols del bloque en memoria
    //                                 
    //                                 std::vector<double> vdB( iRowsB * iColsB );
    //                                 if (transpose_B) {
    //                                     // HDF5 dim1=K, dim2=M  →  read {K_offset, M_offset}
    //                                     dsB->readDatasetBlock( {vstartK[kk], vstartM[jj]}, {vsizetoReadK[kk], vsizetoReadM[jj]}, stride, block, vdB.data() );
    //                                 } else {
    //                                     // HDF5 dim1=M, dim2=K  →  read {M_offset, K_offset}
    //                                     dsB->readDatasetBlock( {vstartM[jj], vstartK[kk]}, {vsizetoReadM[jj], vsizetoReadK[kk]}, stride, block, vdB.data() );
    //                                 }
    //                                 Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(),
    //                                      transpose_B ? vsizetoReadK[kk]   : vsizetoReadM[jj],   // rows del bloque en memoria
    //                                      transpose_B ? vsizetoReadM[jj]   : vsizetoReadK[kk]);  // cols del bloque en memoria
    //                                     
    //                                 // C_accumulator += B * A;
    //                                 // Operación según flags de transposición
    //                                 if (!transpose_A && !transpose_B) {
    //                                     C_accumulator += B * A;                           // A * B
    //                                 } else if (transpose_A && !transpose_B) {
    //                                     C_accumulator += B * A.transpose();               // t(A) * B
    //                                 } else if (!transpose_A && transpose_B) {
    //                                     C_accumulator += B.transpose() * A;               // A * t(B)
    //                                 } else {
    //                                     C_accumulator += B.transpose() * A.transpose();   // t(A) * t(B)
    //                                 }
    //                             }
    //                         
    //                             std::vector<double> vdC_final(vsizetoReadM[jj] * vsizetoRead[ii]);
    //                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> C_final_map(vdC_final.data(), vsizetoReadM[jj], vsizetoRead[ii]);
    //                             C_final_map = C_accumulator;
    //                             
    //                             std::vector<hsize_t> offset = {vstartM[jj], vstart[ii]};
    //                             std::vector<hsize_t> count = {vsizetoReadM[jj], vsizetoRead[ii]};
    //                                                             
    //                             dsC->writeDatasetBlock(vdC_final, offset, count, stride, block);
    //                         }
    //                     }
    //                 }
    //             } 
    //             
    //         } else {
    //             throw std::range_error("multiplication error: non-conformable arguments");
    //         }
    // 
    //     }  catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
    //         // checkClose_file(dsA, dsB, dsC);
    //         throw std::runtime_error("c++ c++ exception multiplication (File IException)");
    //         // return void();
    //     } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
    //         // checkClose_file(dsA, dsB, dsC);
    //         throw std::runtime_error("c++ exception multiplication (DataSet IException)");
    //         // return void();
    //     } catch(std::exception &ex) {
    //         // checkClose_file(dsA, dsB, dsC);
    //         throw std::runtime_error(std::string("c++ exception multiplication: ") + ex.what());
    //         // return void();
    //     }  catch (...) {
    //         // checkClose_file(dsA, dsB, dsC);
    //         throw std::runtime_error("C++ exception multiplication (unknown reason)");
    //         // return void();
    //     }
    // 
    //     return void();
    // }
    
}

#endif // BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_HPP
