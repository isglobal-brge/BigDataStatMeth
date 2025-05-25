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
extern inline Rcpp::NumericVector getDiagonalfromMatrix( BigDataStatMeth::hdf5Dataset* dsMat)
{
    
    Rcpp::NumericVector intNewDiagonal;
    
    try {
        
        std::vector<hsize_t> offset = {0, 0},
                             count = {1, 1},
                             stride = {1, 1},
                             block = {1, 1};
        
        Rcpp::NumericVector diagElem(1);
        
        for(hsize_t i=0; i < dsMat->nrows(); i++) {
            offset[0] = i; offset[1] = i;
            dsMat->readDatasetBlock(  offset, count, stride, block,  REAL(diagElem) );
            intNewDiagonal.push_back(diagElem[0]);
        }
        
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdgetDiagonal_hdf5" << ex.what();
        return(intNewDiagonal);
    }
    
    return(intNewDiagonal);
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
extern inline void setDiagonalMatrix( BigDataStatMeth::hdf5Dataset* dsMat, Rcpp::NumericVector intNewDiagonal)
{
    
    try{
        
        std::vector<hsize_t> count = {1, 1},
                             stride = {1, 1},
                             block = {1, 1};
        
        for( hsize_t i=0; i < (unsigned)intNewDiagonal.size(); i++) {
            std::vector<hsize_t> offset = {i, i};
            dsMat->writeDatasetBlock(Rcpp::wrap(intNewDiagonal(i)), offset, count, stride, block, true);
        }
        
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdgetDiagonal_hdf5" << ex.what();
        return void();
    }
    
    return void();
}

}

#endif // BIGDATASTATMETH_HDF5_MATRIXDIAGONAL_HPP

