/**
 * @file matrixDiagonal.hpp
 * @brief Functions for manipulating matrix diagonals in HDF5 datasets
 *
 * This file provides functionality for extracting and setting diagonal elements
 * of matrices stored in HDF5 format. The implementation is optimized for
 * efficient access to diagonal elements without loading the entire matrix
 * into memory.
 *
 * Key features:
 * - Diagonal extraction from HDF5 matrices
 * - Diagonal element setting for HDF5 matrices
 * - Memory-efficient block-wise operations
 * - Exception-safe implementation
 *
 * @note These operations are particularly useful for large matrices where
 * loading the entire matrix into memory is not feasible.
 *
 * @see BigDataStatMeth::hdf5Dataset
 */

#ifndef BIGDATASTATMETH_HDF5_MATRIXDIAGONAL_HPP
#define BIGDATASTATMETH_HDF5_MATRIXDIAGONAL_HPP

// #include <RcppEigen.h>
// #include "H5Cpp.h"

namespace BigDataStatMeth {
    
    /**
     * @brief Extracts the diagonal elements from a matrix stored in HDF5 format
     *
     * @param dsMat Input HDF5 dataset containing the matrix
     * @return Rcpp::NumericVector Vector containing the diagonal elements
     *
     * @details Implementation details:
     * - Reads diagonal elements one at a time to minimize memory usage
     * - Uses HDF5 block reading for efficient access
     * - Returns empty vector on error
     *
     * @note This function is optimized for matrices where reading the entire
     * matrix into memory would be impractical.
     *
     * @throws std::exception on HDF5 read errors or memory allocation failures
     */
    inline Rcpp::NumericVector getDiagonalfromMatrix(BigDataStatMeth::hdf5Dataset* dsMat)
    {
        hsize_t matrix_size = dsMat->nrows();
        Rcpp::NumericVector diagonal(matrix_size);
        
        try {
            const hsize_t DIAG_BLOCK_SIZE = 256;
            std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
            
            for (hsize_t block_start = 0; block_start < matrix_size; block_start += DIAG_BLOCK_SIZE) {
                hsize_t current_block_size = std::min(DIAG_BLOCK_SIZE, matrix_size - block_start);
                
                // Read square block starting from diagonal position  
                std::vector<hsize_t> offset = {block_start, block_start};
                std::vector<hsize_t> count = {current_block_size, current_block_size};
                
                std::vector<double> block_data(current_block_size * current_block_size);
                dsMat->readDatasetBlock(offset, count, stride, block, block_data.data());
                
                // Map to Eigen for correct R/HDF5 layout handling
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> 
                    block_matrix(block_data.data(), current_block_size, current_block_size);
                
                // Use Eigen diagonal extraction - NO LOOP NEEDED
                Eigen::Map<Eigen::VectorXd> diagonal_segment(REAL(diagonal) + block_start, current_block_size);
                diagonal_segment = block_matrix.diagonal();
            }
            
        } catch(std::exception& ex) {
            Rcpp::stop ("c++ exception getDiagonalfromMatrix: %s", ex.what());
            // Rcpp::Rcout << "c++ exception getDiagonalfromMatrix: " << ex.what();
        }
        
        return diagonal;
    }
    
    
    
    /**
     * @brief Sets the diagonal elements of a matrix stored in HDF5 format
     *
     * @param dsMat Target HDF5 dataset containing the matrix
     * @param intNewDiagonal Vector of new diagonal values to set
     *
     * @details Implementation approach:
     * - Writes diagonal elements one at a time
     * - Uses HDF5 block writing for efficient access
     * - Preserves existing non-diagonal elements
     *
     * Usage example:
     * @code
     * BigDataStatMeth::hdf5Dataset* matrix = ...;
     * Rcpp::NumericVector newDiag = ...;
     * setDiagonalMatrix(matrix, newDiag);
     * @endcode
     *
     * @throws std::exception on HDF5 write errors or dimension mismatch
     */
    inline void setDiagonalMatrix(BigDataStatMeth::hdf5Dataset* dsMat, Rcpp::NumericVector intNewDiagonal)
    {
        try {
            const hsize_t DIAG_BLOCK_SIZE = 256;
            hsize_t matrix_size = intNewDiagonal.size();
            std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
            
            for (hsize_t block_start = 0; block_start < matrix_size; block_start += DIAG_BLOCK_SIZE) {
                hsize_t N = std::min(DIAG_BLOCK_SIZE, matrix_size - block_start);
                
                std::vector<hsize_t> offset = {block_start, block_start};
                std::vector<hsize_t> count  = {N, N};
                
                // Read block as raw HDF5 row-major bytes.
                // block_data[i*N + j] == HDF5[block_start+i, block_start+j]
                std::vector<double> block_data(N * N);
                dsMat->readDatasetBlock(offset, count, stride, block, block_data.data());
                
                // Diagonal elements sit at block_data[k*N + k] (same index in
                // HDF5 row-major and R column-major for the diagonal).
                // Write new diagonal values directly — no transposition needed.
                for (hsize_t k = 0; k < N; ++k)
                    block_data[k * N + k] = intNewDiagonal[block_start + k];
                
                // Write back using the vector overload, which passes the raw
                // row-major bytes straight to HDF5 without any R/Rcpp wrapping.
                dsMat->writeDatasetBlock(block_data, offset, count, stride, block);
            }
            
        } catch(std::exception& ex) {
            Rcpp::stop("c++ exception setDiagonalMatrix: "+ std::string(ex.what()));
        }
    }

}

#endif // BIGDATASTATMETH_HDF5_MATRIXDIAGONAL_HPP

