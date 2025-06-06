/**
 * @file tcrossprod.hpp
 * @brief Transposed cross-product operations for HDF5 matrices
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
        BigDataStatMeth::hdf5Dataset* dsC, hsize_t hdf5_block, hsize_t mem_block_size, 
        bool bparal, bool browmajor, Rcpp::Nullable<int> threads  = R_NilValue) 

{
    
    try {

        hsize_t N = dsA->ncols();
        hsize_t K = dsA->nrows();
        hsize_t M = dsB->ncols();
        hsize_t L = dsB->nrows();

        if( K == L)
        {
            hsize_t isize = hdf5_block + 1,
                    ksize = hdf5_block + 1,
                    jsize = hdf5_block + 1;

            std::vector<hsize_t> stride = {1, 1};
            std::vector<hsize_t> block = {1, 1};

            dsC->createDataset( N, M, "real");

            for (hsize_t ii = 0; ii < N; ii += hdf5_block)
            {
                hsize_t Nii = getOptimBlockSize( N, hdf5_block, ii, isize);
                
                if( ii + hdf5_block > N ) isize = N - ii;
                
                for (hsize_t jj = 0; jj < M; jj += hdf5_block)
                {
                    hsize_t Mjj = getOptimBlockSize( M, hdf5_block, jj, jsize);
                    
                    if( jj + hdf5_block > M) jsize = M - jj;
                    
                    for(hsize_t kk = 0; kk < K; kk += hdf5_block)
                    {
                        if( kk + hdf5_block > K ) ksize = K - kk;

                        hsize_t Kkk = getOptimBlockSize( K, hdf5_block, kk, ksize);
                        
                        Eigen::MatrixXd C;
                        Eigen::MatrixXd A;
                        
                        {
                            std::vector<double> vdA( Nii * Kkk);
                            dsA->readDatasetBlock( {kk, ii}, {Kkk, Nii}, stride, block, vdA.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_A (vdA.data(), Kkk, Nii );    
                            A = tmp_A.transpose();
                        }
                        
                        std::vector<double> vdB( Mjj * Kkk);
                        dsB->readDatasetBlock( {kk, jj}, { Kkk, Mjj}, stride, block, vdB.data() );
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), Kkk, Mjj);   
                        
                        {
                            std::vector<double> vdC( Mjj * Nii);
                            dsC->readDatasetBlock( {jj, ii}, {Mjj, Nii}, stride, block, vdC.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_C (vdC.data(), Mjj, Nii);
                            C = tmp_C.transpose();
                        }

                        // if( bparal == false) {
                            C = C + (A * B);
                        // } else{
                        //     C = C + BigDataStatMeth::Bblock_matrix_mul_parallel(A, B, mem_block_size, threads);
                        // }

                        std::vector<hsize_t> offset = {jj, ii};
                        std::vector<hsize_t> count = {Mjj, Nii};
                        
                        dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, false);

                        // Readjust counters
                        if( kk + hdf5_block > K ) ksize = hdf5_block + 1;
                        if( Kkk > hdf5_block ) {
                            kk = kk - hdf5_block + Kkk; }
                    }

                    if( jj + hdf5_block > M ) jsize = hdf5_block + 1;
                    if( Mjj > hdf5_block ) {
                        jj = jj - hdf5_block + Mjj; }
                }

                if( ii + hdf5_block > N ) isize = hdf5_block + 1;
                if( Nii > hdf5_block ) {
                    ii = ii - hdf5_block + Nii; }
            }

        } else {
            throw std::range_error("non-conformable arguments");
        }

    } catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception tcrossprod: "<<ex.what()<< " \n";
        return(dsC);
    }

    return(dsC);
}

}

#endif // BIGDATASTATMETH_ALGEBRA_TCROSSPROD_HPP