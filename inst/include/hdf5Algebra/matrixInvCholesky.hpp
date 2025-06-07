/**
 * @file matrixInvCholesky.hpp
 * @brief Implementation of Cholesky decomposition and matrix inversion using HDF5
 *
 * This file provides functionality for computing matrix inverses using the Cholesky
 * decomposition method, specifically designed for large matrices stored in HDF5 format.
 * The implementation includes parallel processing capabilities and memory-efficient
 * block operations.
 *
 * Key features:
 * - Cholesky decomposition for positive definite matrices
 * - Matrix inversion using Cholesky decomposition
 * - Block-wise processing for memory efficiency
 * - Parallel computation support
 * - HDF5 integration for large matrix handling
 *
 * @note This implementation is particularly efficient for large, symmetric,
 * positive-definite matrices that don't fit in memory.
 *
 * @see BigDataStatMeth::hdf5Dataset
 * @see BigDataStatMeth::hdf5DatasetInternal
 */

#ifndef BIGDATASTATMETH_HDF5_INVCHOLESKY_HPP
#define BIGDATASTATMETH_HDF5_INVCHOLESKY_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

namespace BigDataStatMeth {

/**
 * @brief Computes matrix inverse using Cholesky decomposition with HDF5 storage
 *
 * @param inDataset Input matrix dataset (must be symmetric positive-definite)
 * @param outDataset Output dataset for the inverse matrix
 * @param bfull If true, computes full matrix inverse; if false, only lower triangular part
 * @param dElementsBlock Block size for processing (minimum 2 * matrix dimension)
 * @param threads Number of threads for parallel processing (optional)
 *
 * @details This function performs matrix inversion in three steps:
 * 1. Cholesky decomposition (A = LL^T)
 * 2. Inverse of the Cholesky factor (L^-1)
 * 3. Computation of full inverse (A^-1 = L^-T L^-1)
 *
 * Performance considerations:
 * - Time complexity: O(n³) for n×n matrix
 * - Space complexity: O(b²) where b is the block size
 * - Parallel efficiency depends on matrix and block sizes
 *
 * @throws H5::FileIException on HDF5 file operation errors
 * @throws H5::GroupIException on HDF5 group operation errors
 * @throws H5::DataSetIException on HDF5 dataset operation errors
 * @throws std::exception on general errors
 */
inline void Rcpp_InvCholesky_hdf5 ( BigDataStatMeth::hdf5Dataset* inDataset,  BigDataStatMeth::hdf5DatasetInternal* outDataset, bool bfull, long dElementsBlock, Rcpp::Nullable<int> threads );

/**
 * @brief Performs Cholesky decomposition on a matrix stored in HDF5 format
 *
 * @param inDataset Input matrix dataset
 * @param outDataset Output dataset for Cholesky factor L
 * @param idim0 Number of rows
 * @param idim1 Number of columns
 * @param dElementsBlock Block size for processing
 * @param threads Number of threads for parallel processing (optional)
 * @return int 0 on success, 1 if matrix is not positive definite, 2 on HDF5 errors
 *
 * @details Implements block-wise Cholesky decomposition algorithm:
 * - Processes matrix in blocks to manage memory usage
 * - Computes lower triangular factor L where A = LL^T
 * - Uses parallel processing for row computations
 *
 * Implementation notes:
 * - Minimum block size is 2 * matrix dimension
 * - Checks for positive definiteness during computation
 * - Optimized for symmetric matrices
 *
 * @throws H5::FileIException on HDF5 file operation errors
 * @throws H5::GroupIException on HDF5 group operation errors
 * @throws H5::DataSetIException on HDF5 dataset operation errors
 */
inline int Cholesky_decomposition_hdf5( BigDataStatMeth::hdf5Dataset* inDataset,  
                                             BigDataStatMeth::hdf5Dataset* outDataset, 
                                             int idim0, int idim1, long dElementsBlock, 
                                             Rcpp::Nullable<int> threads );

/**
 * @brief Computes inverse of Cholesky factor in-place
 *
 * @param InOutDataset Dataset containing Cholesky factor L (will be overwritten with L^-1)
 * @param idim0 Number of rows
 * @param idim1 Number of columns
 * @param dElementsBlock Block size for processing
 * @param threads Number of threads for parallel processing (optional)
 *
 * @details Implements block-wise inversion of lower triangular matrix:
 * - Processes matrix in blocks for memory efficiency
 * - Overwrites input with its inverse
 * - Uses parallel processing for independent computations
 *
 * @note This function assumes input is a lower triangular matrix from Cholesky decomposition
 *
 * @throws H5::FileIException on HDF5 file operation errors
 * @throws H5::GroupIException on HDF5 group operation errors
 * @throws H5::DataSetIException on HDF5 dataset operation errors
 */
inline void Inverse_of_Cholesky_decomposition_hdf5(  BigDataStatMeth::hdf5Dataset* InOutDataset, 
                                                         int idim0, int idim1, long dElementsBlock, 
                                                         Rcpp::Nullable<int> threads);

/**
 * @brief Computes final matrix inverse using inverted Cholesky factors
 *
 * @param InOutDataset Dataset containing L^-1 (will be overwritten with A^-1)
 * @param idim0 Number of rows
 * @param idim1 Number of columns
 * @param dElementsBlock Block size for processing
 * @param threads Number of threads for parallel processing (optional)
 *
 * @details Computes A^-1 = L^-T L^-1 where L^-1 is the inverse of Cholesky factor:
 * - Uses block matrix multiplication for memory efficiency
 * - Implements parallel processing for block operations
 * - Result is symmetric but only lower triangular part is computed
 *
 * @note Upper triangular part can be filled using setUpperTriangularMatrix if needed
 *
 * @throws H5::FileIException on HDF5 file operation errors
 * @throws H5::GroupIException on HDF5 group operation errors
 * @throws H5::DataSetIException on HDF5 dataset operation errors
 */
inline void Inverse_Matrix_Cholesky_parallel( BigDataStatMeth::hdf5Dataset* InOutDataset, 
                                                   int idim0, int idim1, long dElementsBlock, 
                                                   Rcpp::Nullable<int> threads );



inline void Rcpp_InvCholesky_hdf5 ( BigDataStatMeth::hdf5Dataset* inDataset, 
                           BigDataStatMeth::hdf5DatasetInternal* outDataset, 
                           bool bfull, long dElementsBlock, 
                           Rcpp::Nullable<int> threads = R_NilValue )
{
    
    try{
        
        int nrows = inDataset->nrows();
        int ncols = inDataset->ncols();
        
        int res = Cholesky_decomposition_hdf5(inDataset, outDataset, nrows, ncols, dElementsBlock, threads);
        
        if(res == 0)
        {
            Inverse_of_Cholesky_decomposition_hdf5( outDataset, nrows, ncols, dElementsBlock, threads);
            Inverse_Matrix_Cholesky_parallel( outDataset, nrows, ncols, dElementsBlock, threads);
            
            if( bfull==true ) {
                setUpperTriangularMatrix( outDataset, dElementsBlock);
            }
        }
        
    } catch( H5::FileIException& error ) { 
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr<<"c++ exception Rcpp_InvCholesky_hdf5 (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { 
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr << "c++ exception Rcpp_InvCholesky_hdf5 (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { 
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr << "c++ exception Rcpp_InvCholesky_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr << "c++ exception Rcpp_InvCholesky_hdf5" << ex.what();
        return void();
    } catch (...) {
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr<<"\nC++ exception Rcpp_InvCholesky_hdf5 (unknown reason)";
        return void();
    }
    
    return void();
    
}



inline int Cholesky_decomposition_hdf5( BigDataStatMeth::hdf5Dataset* inDataset,  
                           BigDataStatMeth::hdf5Dataset* outDataset, 
                           int idim0,  int idim1,  long dElementsBlock, 
                           Rcpp::Nullable<int> threads  = R_NilValue)
{
    
    try {
        
        int dimensionSize = idim0,
            readedRows = 0,
            rowstoRead,
            minimumBlockSize;
            // chunk = 1,
            
        bool bcancel = false;
        double sum = 0;
        
        std::vector<hsize_t> offset = {0,0},
                             count = {1,1},
                             stride = {1,1},
                             block = {1,1};
        
        
        // Set minimum elements in block (mandatory : minimum = 2 * longest line)
        if( dElementsBlock < dimensionSize * 2 ) {
            minimumBlockSize = dimensionSize * 2;
        } else {
            minimumBlockSize = dElementsBlock;
        }
        
        // Poso el codi per llegir els blocks aquí i desprès a dins hauria d'anar-hi la j
        if( idim0 == idim1)
        {
            while ( readedRows < dimensionSize ) {
                
                rowstoRead = ( -2 * readedRows - 1 + std::sqrt( pow(2*readedRows, 2) - 4 * readedRows + 8 * minimumBlockSize + 1) ) / 2;
                
                if( readedRows + rowstoRead > idim0) { // Max size bigger than data to read ?
                    rowstoRead = idim0 - readedRows;
                }
                
                offset[0] = readedRows;
                
                if( readedRows == 0) {
                    count[0] = dimensionSize - readedRows;
                    count[1] = rowstoRead;
                } else {
                    offset[0] = offset[0] - 1; // We need results from previous line to compute next line results
                    count[1] = readedRows + rowstoRead;
                    count[0] = dimensionSize - readedRows + 1;
                }
                
                readedRows = readedRows + rowstoRead; // Ho preparem perquè desprès necessitarem llegir a partir de la línea anterior
                
                Eigen::MatrixXd A, L;
                
                std::vector<double> vdA( count[0] * count[1] ); 
                inDataset->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
                A = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdA.data(), count[0], count[1] );
                
                std::vector<double> vdL( count[0] * count[1] ); 
                outDataset->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdL.data() );
                L = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdL.data(), count[0], count[1] );
                
                if( offset[0] != 0 )
                    rowstoRead = rowstoRead + 1;
                
                for (int j = 0; j < rowstoRead; j++)  {
                    
                    if( j + offset[0] == 0) {
                        L(j, j) = std::sqrt(A(j,j));
                    } else {
                        L(j, j + offset[0]) = std::sqrt(A(j,j + offset[0]) - (L.row(j).head(j + offset[0]).array().pow(2).sum() ));    
                    }
                    
                    
#pragma omp parallel for num_threads( get_number_threads(threads, R_NilValue) ) private(sum) shared (A,L,j) schedule(dynamic) if (j < readedRows - chunk)
                    for ( int i = j + 1; i < dimensionSize - offset[0]  ; i++ )
                    {
                        if(bcancel == false) {
                            if( j + offset[0] > 0) {
                                sum = (L.block(i, 0, 1, j + offset[0]).array() * L.block(j, 0, 1, j + offset[0]).array()).array().sum();
                                if( sum != sum ) {
                                    Rcpp::Rcout<<"\n Can't get inverse matrix using Cholesky decomposition matrix is not positive definite\n";
                                    bcancel = true;
                                }
                            } else {
                                sum = 0;
                            }
                            if(!bcancel){
                                L(i,j + offset[0]) =  (1/L(j,j + offset[0])*(A(i,j + offset[0]) - sum));    
                            }    
                        }
                    }
                }
                if(bcancel == true) {
                    return(1);
                }
                
                if( offset[0] != 0) {
                    offset[0] = offset[0] + 1;
                    count[0] = count[0] - 1;
                    
                    outDataset->writeDatasetBlock( Rcpp::wrap(L.block(1, 0, L.rows()-1, L.cols())), offset, count, stride, block, false);
                    
                } else {
                    outDataset->writeDatasetBlock( Rcpp::wrap(L), offset, count, stride, block, false);
                }
                
                offset[0] = offset[0] + count[0] - 1;
                
            }
            
        } else {
            throw std::range_error("non-conformable arguments");
        }
        
    } catch( H5::FileIException& error ) { 
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr<<"c++ exception Cholesky_decomposition_hdf5 (File IException)";
        return(2);
    } catch( H5::GroupIException & error ) { 
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr << "c++ exception Cholesky_decomposition_hdf5 (Group IException)";
        return(2);
    } catch( H5::DataSetIException& error ) { 
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr << "c++ exception Cholesky_decomposition_hdf5 (DataSet IException)";
        return(2);
    } catch(std::exception& ex) {
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr << "c++ exception Cholesky_decomposition_hdf5" << ex.what();
        return(2);
    } catch (...) {
        checkClose_file(inDataset, outDataset);
        Rcpp::Rcerr<<"\nC++ exception Rcpp_InvCholesky_hdf5 (unknown reason)";
        return(2);
    }
    return(0);
    
}




inline void Inverse_of_Cholesky_decomposition_hdf5( BigDataStatMeth::hdf5Dataset* InOutDataset, 
                                    int idim0, int idim1, long dElementsBlock, 
                                    Rcpp::Nullable<int> threads = R_NilValue)
{
    
    try{
        
        int dimensionSize = idim0, 
            readedCols = 0,
            colstoRead,
            minimumBlockSize;
        
        std::vector<hsize_t> offset = {0,0},
                             count = {1,1},
                             stride = {1,1},
                             block = {1,1};

        
        // Set minimum elements in block (mandatory : minimum = 2 * longest line)
        if( dElementsBlock < dimensionSize * 2 ) {
            minimumBlockSize = dimensionSize * 2;
        } else {
            minimumBlockSize = dElementsBlock;
        }
        
        Eigen::VectorXd Diagonal = Rcpp::as<Eigen::VectorXd>(getDiagonalfromMatrix(InOutDataset));
        // Set the new diagonal in the result matrix
        setDiagonalMatrix( InOutDataset, Rcpp::wrap(Diagonal.cwiseInverse()) );
        
        if( idim0 == idim1) {
            
            while ( readedCols < dimensionSize ) {
                
                colstoRead = ( -2 * readedCols - 1 + std::sqrt( pow(2*readedCols, 2) - 4 * readedCols + 8 * minimumBlockSize + 1) ) / 2;
                
                if( readedCols + colstoRead > idim0) { // Max size bigger than data to read ?
                    colstoRead = idim0 - readedCols;
                }
                
                offset[0] = readedCols;
                count[0] =  dimensionSize - offset[0];
                count[1] = readedCols + colstoRead;
                
                Eigen::MatrixXd verticalData;
                
                std::vector<double> vverticalData( count[0] * count[1] ); 
                InOutDataset->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vverticalData.data() );
                verticalData = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vverticalData.data(), count[0], count[1] );
                
                for (int j = 1; j < dimensionSize - offset[0]; j++)
                {
                    
                    Eigen::VectorXd vR = Eigen::VectorXd::Zero(j+offset[0]);
                    Eigen::ArrayXd ar_j;
                    
                    int size_j;
                    if( j < colstoRead) {
                        size_j = j;
                    } else {
                        size_j = colstoRead;
                    }
                    
                    ar_j = verticalData.block( j, offset[0], 1,  size_j).transpose().array();
                    
                    vR = Eigen::VectorXd::Zero(offset[0] + ar_j.size());
                    
#pragma omp parallel for num_threads( get_number_threads(threads, R_NilValue) ) shared (ar_j, j, verticalData, offset, colstoRead, vR) schedule(dynamic) 
                    for (int i = 0; i < offset[0] + ar_j.size() ; i++) {
                        
                        Eigen::ArrayXd ar_i = verticalData.block( 0, i, size_j, 1).array();
                        
                        if( offset[0] > 0 ) {
                            if( j  <= colstoRead ) {
                                if( i < offset[0] ){
                                    vR(i) =  (( verticalData.coeff(j, i) + ((ar_j.transpose() * ar_i.transpose()).sum())) * (-1)) / Diagonal[j+offset[0]];
                                } else {
                                    vR(i) =   ((ar_j.transpose() * ar_i.transpose()) * (-1)).sum() / Diagonal[j+offset[0]];
                                    
                                }
                            } else {
                                if( i < offset[0] ){
                                    vR(i) =  ( verticalData.coeff(j, i) + ((ar_j.transpose() * ar_i.transpose()).sum()));
                                } else {
                                    vR(i) =   (ar_j.transpose() * ar_i.transpose()).sum();
                                }
                            }
                        } else {
                            if( j <= colstoRead ) {
                                vR(i) =   ((ar_j.transpose() * ar_i.transpose()) * (-1)).sum() / Diagonal[j+offset[0]];
                            } else {
                                vR(i) =   (ar_j.transpose() * ar_i.transpose()).sum();
                            }
                        }
                    }
                    
                    verticalData.block(j, 0, 1, vR.size()) = vR.transpose();
                    
                }
                
                InOutDataset->writeDatasetBlock( Rcpp::wrap(verticalData), offset, count, stride, block, false);
                
                readedCols = readedCols + colstoRead; // Ho preparem perquè desprès necessitarem llegir a partir de la línea anterior
                
            }
        } else {
            throw std::range_error("non-conformable arguments");
        }
        
    } catch( H5::FileIException& error ) { 
        checkClose_file(InOutDataset);
        Rcpp::Rcerr<<"c++ exception Inverse_of_Cholesky_decomposition_hdf5 (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { 
        checkClose_file(InOutDataset);
        Rcpp::Rcerr << "c++ exception Inverse_of_Cholesky_decomposition_hdf5 (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { 
        checkClose_file(InOutDataset);
        Rcpp::Rcerr << "c++ exception Inverse_of_Cholesky_decomposition_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        checkClose_file(InOutDataset);
        Rcpp::Rcerr << "c++ exception Inverse_of_Cholesky_decomposition_hdf5" << ex.what();
        return void();
    } catch (...) {
        checkClose_file(InOutDataset);
        Rcpp::Rcerr<<"\nC++ exception Inverse_of_Cholesky_decomposition_hdf5 (unknown reason)";
        return void();
    }
    
    return void();
    
}



inline void Inverse_Matrix_Cholesky_parallel( BigDataStatMeth::hdf5Dataset* InOutDataset, 
                                     int idim0, int idim1, long dElementsBlock, 
                                     Rcpp::Nullable<int> threads = R_NilValue)
{
    
    try {
        
        int dimensionSize = idim0,
            readedCols = 0,
            colstoRead,
            minimumBlockSize;
        
        Eigen::VectorXd newDiag(idim0);
        
        std::vector<hsize_t> offset = {0,0},
                             count = {1,1},
                             stride = {1,1},
                             block = {1,1};
        
        
        // Set minimum elements in block (mandatory : minimum = 2 * longest line)
        if( dElementsBlock < dimensionSize * 2 ) {
            minimumBlockSize = dimensionSize * 2;
        } else {
            minimumBlockSize = dElementsBlock;
        }
        
        if( idim0 == idim1 )
        {
            while ( readedCols < dimensionSize ) {
                
                colstoRead = ( -2 * readedCols - 1 + std::sqrt( pow(2*readedCols, 2) - 4 * readedCols + 8 * minimumBlockSize + 1) ) / 2;
                if(colstoRead ==1) {
                    colstoRead = 2;
                }
                
                if( readedCols + colstoRead > idim0) { // Max size bigger than data to read ?
                    colstoRead = idim0 - readedCols;
                }
                
                offset[0] = readedCols;
                count[0] =  dimensionSize - readedCols;
                count[1] = readedCols + colstoRead;
                
                Eigen::MatrixXd verticalData;
                
                std::vector<double> vverticalData( count[0] * count[1] ); 
                InOutDataset->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vverticalData.data() );
                verticalData = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vverticalData.data(), count[0], count[1] );
                
                
// #pragma omp parallel for num_threads(ithreads) shared (verticalData, colstoRead, offset) schedule(dynamic)
#pragma omp parallel for num_threads( get_number_threads(threads, R_NilValue) ) shared (verticalData, colstoRead, offset) schedule(static) ordered
                for ( int i = 0; i < colstoRead + offset[0]; i++)   // Columns
                {
                    int init;
                    
                    if(offset[0] == 0) {
                        init = i + 1;
                        newDiag(i) = verticalData.block(i, i, idim0-i, 1 ).array().pow(2).sum();
                    } else {
                        if(  i < offset[0]) {
                            init = 0;
                        } else {
                            newDiag(i) = verticalData.block( i-offset[0], i, idim0-i, 1 ).array().pow(2).sum();
                            if( i < offset[0] + colstoRead - 1) {
                                init = i - offset[0] + 1;
                            } else {
                                init = colstoRead; // force end
                            }
                        }
                    }
                    
                    #pragma omp ordered
                    {
                        for ( int j = init; j < colstoRead ; j++) { // Rows
                            
                            if( offset[0] + j < verticalData.cols()) {
                                verticalData(j,i) = (verticalData.block( j, i , verticalData.rows() - j, 1).array() * verticalData.block( j, j + offset[0],  verticalData.rows() - j, 1).array()).sum();
                            }
                        }
                    }
                   
                }
                
                InOutDataset->writeDatasetBlock( Rcpp::wrap(verticalData), offset, count, stride, block, false);
                
                readedCols = readedCols + colstoRead; 
            }
            
            setDiagonalMatrix( InOutDataset, Rcpp::wrap(newDiag) );
            
        } else {
            throw std::range_error("non-conformable arguments");
        }
        
    } catch( H5::FileIException& error ) { 
        checkClose_file(InOutDataset);
        Rcpp::Rcerr<<"c++ exception Inverse_Matrix_Cholesky_parallel (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { 
        checkClose_file(InOutDataset);
        Rcpp::Rcerr << "c++ exception Inverse_Matrix_Cholesky_parallel (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { 
        checkClose_file(InOutDataset);
        Rcpp::Rcerr << "c++ exception Inverse_Matrix_Cholesky_parallel (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        checkClose_file(InOutDataset);
        Rcpp::Rcerr << "c++ exception Inverse_Matrix_Cholesky_parallel: " << ex.what();
        return void();
    } catch (...) {
        checkClose_file(InOutDataset);
        Rcpp::Rcerr<<"\nC++ exception Inverse_Matrix_Cholesky_parallel (unknown reason)";
        return void();
    }
    
    return void();
}



}

#endif // BIGDATASTATMETH_HDF5_INVCHOLESKY_HPP

