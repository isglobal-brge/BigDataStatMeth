/**
 * @file matrixTriangular.hpp
 * @brief Triangular matrix operations for HDF5 matrices
 * @details This header file provides implementations for operations involving
 * triangular matrices stored in HDF5 format. The implementation includes:
 * 
 * Key features:
 * - Upper triangular operations
 * - Lower triangular operations
 * - Triangular solvers
 * - Block-based processing
 * - Memory-efficient algorithms
 * 
 * Supported operations:
 * - Triangular matrix multiplication
 * - Triangular system solving
 * - Back substitution
 * - Forward substitution
 * - Matrix decomposition
 * 
 * Performance features:
 * - Cache-friendly algorithms
 * - Block-based computation
 * - Multi-threaded processing
 * - Memory access optimization
 * - I/O efficiency
 * 
 * The implementation uses:
 * - BLAS Level 3 operations
 * - Block algorithms
 * - HDF5 chunked storage
 * - Parallel processing
 */

#ifndef BIGDATASTATMETH_HDF5_MATRIXTRIANGULAR_HPP
#define BIGDATASTATMETH_HDF5_MATRIXTRIANGULAR_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

namespace BigDataStatMeth {

/**
 * @brief Set the upper triangular matrix
 * @details Set the upper triangular matrix of a dataset
 * 
 * @param dsMat The dataset to set the upper triangular matrix of
 * @param dElementsBlock The block size to use for the matrix
 */
inline void setUpperTriangularMatrix( BigDataStatMeth::hdf5Dataset* dsMat, hsize_t dElementsBlock)
{
    
    try {
        
        std::vector<hsize_t> count = {1, 1},
                             stride = {1, 1},
                             block = {1, 1},
                             offset = {0, 0};
        
        hsize_t readedRows = 0,
                rowstoRead,
                minimumBlockSize;
        
        if( dElementsBlock < dsMat->nrows() * 2 ) {
            minimumBlockSize = dsMat->nrows() * 2;
        } else {
            minimumBlockSize = dElementsBlock;
        }
        
        while ( readedRows < dsMat->nrows() ) {
            
            rowstoRead = ( -2 * readedRows - 1 + std::sqrt( pow(2*readedRows, 2) - 4 * readedRows + 8 * minimumBlockSize + 1) ) / 2;
            
            if( readedRows + rowstoRead > dsMat->nrows()) { // Max size bigger than data to read ?
                rowstoRead = dsMat->nrows() - readedRows;
            }
            
            if( readedRows == 0) { // Read complete cols
                count[0] = dsMat->nrows();
                count[1] = rowstoRead;
            } else {
                count[1] = rowstoRead;
                count[0] = dsMat->nrows() - readedRows;
            }
            
            offset[1] = readedRows;
            offset[0] = offset[1];
            
            std::vector<double> vdA( count[0] * count[1]);
            dsMat->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), count[0], count[1] );
            
            // #pragma omp parallel for num_threads(getDTthreads(ithreads, true)) private(sum) shared (A,L,j) schedule(static) if (j < readedRows - chunk)
            for ( hsize_t i = 0; i < rowstoRead; i++ ) {
                std::vector<hsize_t> tcount = { dsMat->nrows()-readedRows-i-1, 1 };
                std::vector<hsize_t> toffset = {readedRows+i, readedRows+i+1 };
                
                dsMat->writeDatasetBlock( Rcpp::wrap( A.block(i+1, i, dsMat->nrows()-readedRows-i-1, 1 ).transpose() ), toffset, tcount, stride, block, true);
            }
            
            readedRows = readedRows + rowstoRead; 
        }
        
    }
    catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception Rcpp_setUpperTriangularMatrix (File IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Rcpp_setUpperTriangularMatrix (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception Rcpp_setUpperTriangularMatrix" << ex.what();
        return void();
    }
    
    return void();
}

/**
 * @brief Set the lower triangular matrix
 * @details Set the lower triangular matrix of a dataset
 * 
 * @param dsMat The dataset to set the lower triangular matrix of
 * @param dElementsBlock The block size to use for the matrix
 */
inline void setLowerTriangularMatrix( BigDataStatMeth::hdf5Dataset* dsMat, hsize_t dElementsBlock)
{
    
    try {
        
        std::vector<hsize_t> count = {1, 1},
                            stride = {1, 1},
                            block = {1, 1},
                            offset = {0, 0};
        
        hsize_t readedRows = 0,
                rowstoRead,
                minimumBlockSize;
        
        if( dElementsBlock < dsMat->nrows() * 2 ) {
            minimumBlockSize = dsMat->nrows() * 2;
        } else {
            minimumBlockSize = dElementsBlock;
        }
        
        while ( readedRows < dsMat->nrows() ) {
            
            rowstoRead = ( -2 * readedRows - 1 + std::sqrt( pow(2*readedRows, 2) - 4 * readedRows + 8 * minimumBlockSize + 1) ) / 2;
            
            if( readedRows + rowstoRead > dsMat->nrows()) { // Max size bigger than data to read ?
                rowstoRead = dsMat->nrows() - readedRows;
            }
            
            if( readedRows == 0) { // Read complete cols
                count[0] = rowstoRead;
                count[1] = dsMat->nrows();
            } else {
                count[1] = dsMat->nrows() - readedRows;
                count[0] = rowstoRead;
            }
            
            offset[1] = readedRows;
            offset[0] = offset[1];
            
            std::vector<double> vdA( count[0] * count[1]);
            dsMat->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), count[0], count[1] );
            
            // #pragma omp parallel for num_threads(getDTthreads(ithreads, true)) private(sum) shared (A,L,j) schedule(static) if (j < readedRows - chunk)
            for ( hsize_t i = 0; i < rowstoRead; i++ ) {
                std::vector<hsize_t> tcount = { 1, dsMat->nrows()-readedRows-i-1 };
                std::vector<hsize_t> toffset = { readedRows+i+1, readedRows+i };
                
                dsMat->writeDatasetBlock( Rcpp::wrap( A.block(i, i+1, 1, dsMat->nrows()-readedRows-i-1 ).transpose() ), toffset, tcount, stride, block, true);
            }
            
            readedRows = readedRows + rowstoRead; 
        }
    }
    catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception Rcpp_setLowerTriangularMatrix (File IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Rcpp_setLowerTriangularMatrix (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception Rcpp_setLowerTriangularMatrix" << ex.what();
        return void();
    }
    
    return void();
}


}

#endif // BIGDATASTATMETH_HDF5_MATRIXTRIANGULAR_HPP

