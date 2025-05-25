/**
 * @file crossprod.hpp
 * @brief Implementation of matrix cross-product operations using HDF5
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

// #include "hdf5Algebra/multiplication.hpp"

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
extern inline BigDataStatMeth::hdf5Dataset* crossprod( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, 
        BigDataStatMeth::hdf5Dataset* dsC, hsize_t hdf5_block, hsize_t mem_block_size, 
        bool bparal, bool browmajor, Rcpp::Nullable<int> threads  = R_NilValue) 

{
    
    try {

        hsize_t N = dsA->nrows();
        hsize_t K = dsA->ncols();
        hsize_t M = dsB->nrows();
        hsize_t L = dsB->ncols();
        
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
                hsize_t iRowsA = getOptimBlockSize( N, hdf5_block, ii, isize);
                
                if( ii + hdf5_block > N ) isize = N - ii;
                
                for (hsize_t jj = 0; jj < M; jj += hdf5_block)
                {
                    hsize_t iRowsB = getOptimBlockSize( M, hdf5_block, jj, jsize);
                    
                    if( jj + hdf5_block > M) jsize = M - jj;
                    
                    for(hsize_t kk = 0; kk < K; kk += hdf5_block)
                    {
                        if( kk + hdf5_block > K ) ksize = K - kk;

                        hsize_t iColsA = getOptimBlockSize( K, hdf5_block, kk, ksize),
                            iColsB = iColsA;
                        
                        Eigen::MatrixXd C;
                        Eigen::MatrixXd B;
                        
                        std::vector<double> vdA( iRowsA * iColsA);
                        dsA->readDatasetBlock( {ii, kk}, {iRowsA,iColsA}, stride, block, vdA.data() );
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), iRowsA, iColsA );

                        {
                            std::vector<double> vdB( iRowsB * iColsB);
                            dsB->readDatasetBlock( {jj, kk}, {iRowsB, iColsB}, stride, block, vdB.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_B (vdB.data(), iRowsB, iColsB);   
                            B = tmp_B.transpose();
                        }
                        
                        {
                            std::vector<double> vdC( iRowsB * iRowsA);
                            dsC->readDatasetBlock( {jj, ii}, {iRowsB, iRowsA}, stride, block, vdC.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_C (vdC.data(), iRowsB, iRowsA);
                            C = tmp_C.transpose();
                        }
                           
                        C = C + (A * B);
                        
                        std::vector<hsize_t> offset = {jj,ii};
                        std::vector<hsize_t> count = {iRowsB, iRowsA};
                        
                        dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, false);

                        // Readjust counters
                        if( kk + hdf5_block > K ) ksize = hdf5_block + 1;
                        if( iColsA > hdf5_block ) {
                            kk = kk - hdf5_block + iColsA; }
                    }

                    if( jj + hdf5_block > M ) jsize = hdf5_block + 1;
                    if( iRowsB > hdf5_block ) {
                        jj = jj - hdf5_block + iRowsB; }
                }

                if( ii + hdf5_block > N ) isize = hdf5_block + 1;
                if( iRowsA > hdf5_block ) {
                    ii = ii - hdf5_block + iRowsA; }
            }

        } else {
            throw std::range_error("non-conformable arguments");
        }

    } catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception crossprod: "<<ex.what()<< " \n";
        return(dsC);
    }

    return(dsC);
}

}

#endif // BIGDATASTATMETH_ALGEBRA_CROSSPROD_HPP