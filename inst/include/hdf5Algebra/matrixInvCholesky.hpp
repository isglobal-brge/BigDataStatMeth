#ifndef BIGDATASTATMETH_HDF5_INVCHOLESKY_HPP
#define BIGDATASTATMETH_HDF5_INVCHOLESKY_HPP

// #include <RcppEigen.h>
// #include "H5Cpp.h"

#include "Utilities/openme-utils.hpp"
#include "hdf5Algebra/matrixDiagonal.hpp"
#include "hdf5Algebra/matrixTriangular.hpp"


namespace BigDataStatMeth {

extern inline void Rcpp_InvCholesky_hdf5 ( BigDataStatMeth::hdf5Dataset* inDataset,  BigDataStatMeth::hdf5DatasetInternal* outDataset, bool bfull, long dElementsBlock, Rcpp::Nullable<int> threads );
extern inline int Cholesky_decomposition_hdf5( BigDataStatMeth::hdf5Dataset* inDataset,  BigDataStatMeth::hdf5Dataset* outDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads );
extern inline void Inverse_of_Cholesky_decomposition_hdf5(  BigDataStatMeth::hdf5Dataset* InOutDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads);
extern inline void Inverse_Matrix_Cholesky_parallel( BigDataStatMeth::hdf5Dataset* InOutDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads);

extern inline void Rcpp_InvCholesky_hdf5 ( BigDataStatMeth::hdf5Dataset* inDataset, BigDataStatMeth::hdf5DatasetInternal* outDataset, bool bfull, long dElementsBlock, Rcpp::Nullable<int> threads = R_NilValue )
{
    
    try{
        
        int nrows = inDataset->nrows();
        int ncols = inDataset->ncols();
        
        int res = Cholesky_decomposition_hdf5(inDataset, outDataset, nrows, ncols, dElementsBlock, threads);
        // delete dsA;
        
        if(res == 0)
        {
            Inverse_of_Cholesky_decomposition_hdf5( outDataset, nrows, ncols, dElementsBlock, threads);
            Inverse_Matrix_Cholesky_parallel( outDataset, nrows, ncols, dElementsBlock, threads); 
            
            if( bfull==true ) {
                setUpperTriangularMatrix( outDataset, dElementsBlock);
            }
        }
        
    }catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception Rcpp_InvCholesky_hdf5 (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Rcpp_InvCholesky_hdf5 (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Rcpp_InvCholesky_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception Rcpp_InvCholesky_hdf5" << ex.what();
        return void();
    }
    
    return void();
    
}



//..// void Cholesky_decomposition_hdf5( H5File* file, DataSet* inDataset, DataSet* outDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads )
extern inline int Cholesky_decomposition_hdf5( BigDataStatMeth::hdf5Dataset* inDataset,  BigDataStatMeth::hdf5Dataset* outDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads )
{
    
    try {
        
        int dimensionSize = idim0,
            readedRows = 0,
            chunk = 1,
            rowstoRead,
            minimumBlockSize;
        bool bcancel = false;
        // i,k;
        double sum = 0;
        unsigned int ithreads;
        
        std::vector<hsize_t> offset = {0,0},
                             count = {1,1},
                             stride = {1,1},
                             block = {1,1};
        
        if(threads.isNotNull()) {
            if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
                ithreads = Rcpp::as<int> (threads);
            } else {
                ithreads = getDTthreads(0, true);
            }
        } else {
            ithreads = getDTthreads(0, true);
        }
        
        // Set minimum elements in block (mandatory : minimum = 2 * longest line)
        if( dElementsBlock < dimensionSize * 2 ) {
            minimumBlockSize = dimensionSize * 2;
        } else {
            minimumBlockSize = dElementsBlock;
        }
        
        // Poso el codi per llegir els blocks aquí i desprès a dins hauria d'anar-hi la j
        if( idim0 == idim1)
        {
            // llegir files : minimumBlockSize <=
            //      - nfiles ??
            //      - columnes : mínim j, máxim ??
            //
            //  Primer block :
            //      - nfiles^2 = maxelements   =>   nfiles = sqrt(maxelements)
            //  Blocs successius :
            //      - nfiles = maxelements (nfiles anterior +  )
            //      ((nfiles + (nfiles + x)) /2) * (x + 1)
            
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
                    
#pragma omp parallel for num_threads(ithreads) private(sum) shared (A,L,j) schedule(dynamic) if (j < readedRows - chunk)
                    for ( int i = j + 1; i < (dimensionSize - offset[0] && bcancel == false)  ; i++ )
                    {
                        if( j + offset[0] > 0) {
                            sum = (L.block(i, 0, 1, j + offset[0]).array() * L.block(j, 0, 1, j + offset[0]).array()).array().sum();
                            if( sum != sum ) {
                                Rcpp::Rcout<<"\n Can't get inverse matrix using Cholesky decomposition matrix is not positive definite\n";
                                //..// return(1);
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
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception Cholesky_decomposition_hdf5 (File IException)";
        return(2);
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Cholesky_decomposition_hdf5 (Group IException)";
        return(2);
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Cholesky_decomposition_hdf5 (DataSet IException)";
        return(2);
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception Cholesky_decomposition_hdf5" << ex.what();
        return(2);
    }
    return(0);
    
}






//..// void Inverse_of_Cholesky_decomposition_hdf5(  H5File* file, DataSet* InOutDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads = R_NilValue)
extern inline void Inverse_of_Cholesky_decomposition_hdf5(  BigDataStatMeth::hdf5Dataset* InOutDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads = R_NilValue)
{
    
    try{
        
        unsigned int ithreads;
        int dimensionSize = idim0, 
            readedCols = 0,
            // chunk = 1,
            colstoRead,
            minimumBlockSize;
        
        std::vector<hsize_t> offset = {0,0},
                             count = {1,1},
                             stride = {1,1},
                             block = {1,1};
        
        if(threads.isNotNull()) {
            if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()) {
                ithreads = Rcpp::as<int> (threads);
            } else {
                ithreads = getDTthreads(0, true);
            }
        } else {
            ithreads = getDTthreads(0, true);
        }
        
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
                //..// verticalData = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vverticalData.data(), count[1], count[0] );
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
                    
#pragma omp parallel for num_threads(ithreads) shared (ar_j, j, verticalData, offset, colstoRead, vR) schedule(dynamic) 
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
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception Inverse_of_Cholesky_decomposition_hdf5 (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5 (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5" << ex.what();
        return void();
    }
    
    return void();
    
}



extern inline void Inverse_Matrix_Cholesky_parallel( BigDataStatMeth::hdf5Dataset* InOutDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads = R_NilValue)
{
    
    try {
        
        int dimensionSize = idim0,
            readedCols = 0,
            colstoRead,
            minimumBlockSize;
        
        Eigen::VectorXd newDiag(idim0);
        unsigned int ithreads;
        
        std::vector<hsize_t> offset = {0,0},
                             count = {1,1},
                             stride = {1,1},
                             block = {1,1};
        
        if(threads.isNotNull()) {
            if (Rcpp::as<unsigned int> (threads) <= std::thread::hardware_concurrency()){
                ithreads = Rcpp::as<int> (threads);
            } else {
                ithreads = getDTthreads(0, true);
            }
        } else {
            ithreads = getDTthreads(0, true);
        }
        
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
                
                //..// Eigen::MatrixXd verticalData = GetCurrentBlock_hdf5(file, InOutDataset, offset[0], offset[1], count[0], count[1]);
                
                Eigen::MatrixXd verticalData;
                
                std::vector<double> vverticalData( count[0] * count[1] ); 
                InOutDataset->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vverticalData.data() );
                verticalData = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vverticalData.data(), count[0], count[1] );
                
                
#pragma omp parallel for num_threads(ithreads) shared (verticalData, colstoRead, offset) schedule(dynamic)
                for ( int i = 0; i < colstoRead + offset[0]; i++)   // Columnes
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
                    
                    for ( int j = init; j < colstoRead ; j++) {     // Files
                        
                        if( offset[0] + j < verticalData.cols()) {
                            verticalData(j,i) = (verticalData.block( j, i , verticalData.rows() - j, 1).array() * verticalData.block( j, j + offset[0],  verticalData.rows() - j, 1).array()).sum();
                        }
                    }
                }
                
                InOutDataset->writeDatasetBlock( Rcpp::wrap(verticalData), offset, count, stride, block, false);
                
                readedCols = readedCols + colstoRead; // Ho preparem perquè desprès necessitarem llegir a partir de la línea anterior
            }
            
            setDiagonalMatrix( InOutDataset, Rcpp::wrap(newDiag) );
            
        } else {
            throw std::range_error("non-conformable arguments");
        }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception Inverse_Matrix_Cholesky_parallel (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Inverse_Matrix_Cholesky_parallel (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception Inverse_Matrix_Cholesky_parallel (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception Inverse_Matrix_Cholesky_parallel" << ex.what();
        return void();
    }
    
    return void();
}



}

#endif // BIGDATASTATMETH_HDF5_INVCHOLESKY_HPP

