#include "include/hdf5_invCholesky.h"


// Suma :
//  Primer element = a
//  Differència entre elements = 1
//  Nombre de termes = n 
//          sum = n/2[2a + (n-1)*d]
//  Tening en compte que d = 1 quedaria : 
//          sum = n/2[2a + (n-1)]
//  a mi m'interessa saber n sabent sum = size
int get_rowsinBlock( int start, int size) {
    
    int iRowsCols;
    
    iRowsCols = std::ceil(-2 * start + 1 * std::sqrt( 4*std::pow(start,2) - 4 * start + 1 + 16 * size));
    iRowsCols = std::ceil(iRowsCols / 2);
    return(iRowsCols);
    
}



void Cholesky_decomposition_hdf5( H5File* file, DataSet* inDataset, DataSet* outDataset, int idim0, int idim1, double dElementsBlock, Rcpp::Nullable<int> threads )
{

    try {

        int dimensionSize = idim0,
            readedRows = 0,
            chunk = 1,
            rowstoRead,
            minimumBlockSize,
            i,k;
        double sum = 0;
        unsigned int ithreads;
        
        Rcpp::IntegerVector offset = Rcpp::IntegerVector::create(0, 0),
                            count = Rcpp::IntegerVector::create(1, 1),
                            stride = Rcpp::IntegerVector::create(1, 1),
                            block = Rcpp::IntegerVector::create(1, 1);

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

                Eigen::MatrixXd A = GetCurrentBlock_hdf5(file, inDataset, offset[0], offset[1], count[0], count[1]);
                Eigen::MatrixXd L = GetCurrentBlock_hdf5(file, outDataset, offset[0], offset[1], count[0], count[1]);

                if( offset[0] != 0 )
                    rowstoRead = rowstoRead + 1;
                
                for (int j = 0; j < rowstoRead; j++)  {
                    if( j + offset[0] == 0) {
                        L(j, j) = std::sqrt(A(j,j));
                    } else {
                        L(j, j + offset[0]) = std::sqrt(A(j,j + offset[0]) - (L.row(j).head(j + offset[0]).array().pow(2).sum() ));    
                    }

    #pragma omp parallel for num_threads(getDTthreads(ithreads, true)) private(i,sum) shared (A,L,j) schedule(static) if (j < readedRows - chunk)
                    for ( int i = j + 1; i < dimensionSize - offset[0]  ; i++ )
                    {
                        if( j + offset[0] > 0) {
                            sum = (L.block(i, 0, 1, j + offset[0]).array() * L.block(j, 0, 1, j + offset[0]).array()).array().sum();
                        } else {
                            sum = 0;
                        }
                        L(i,j + offset[0]) =  (1/L(j,j + offset[0])*(A(i,j + offset[0]) - sum));    
                    }
                }
                            
                if( offset[0] != 0) {
                    offset[0] = offset[0] + 1;
                    count[0] = count[0] - 1;
                    write_HDF5_matrix_subset_v2( file, outDataset, offset, count, stride, block, Rcpp::wrap( L.block(1, 0, L.rows()-1, L.cols()) ) );
                    
                } else {
                    write_HDF5_matrix_subset_v2( file, outDataset, offset, count, stride, block, Rcpp::wrap( L ) );
                }
                    
                offset[0] = offset[0] + count[0] - 1;
                
            }

        } else {
            throw std::range_error("non-conformable arguments");
        }

    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        outDataset->close();
        inDataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception Cholesky_decomposition_hdf5 (File IException)";
        return void();
    } catch( GroupIException & error ) { // catch failure caused by the DataSet operations
        outDataset->close();
        inDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Cholesky_decomposition_hdf5 (Group IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        outDataset->close();
        inDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Cholesky_decomposition_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        outDataset->close();
        inDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Cholesky_decomposition_hdf5" << ex.what();
        return void();
    }
    return void();

}



// void Inverse_of_Cholesky_decomposition_hdf5(  H5File* file, DataSet* InOutDataset, int idim0, int idim1, double dElementsBlock, Rcpp::Nullable<int> threads = R_NilValue)
// {
//     
//     try{
//         
//         unsigned int ithreads;
//         int dimensionSize = idim0, 
//             readedRows = 0,
//             chunk = 1,
//             rowstoRead,
//             minimumBlockSize,
//             i, k;
//         
//         Rcpp::IntegerVector offset = Rcpp::IntegerVector::create(0, 0),
//                             count = Rcpp::IntegerVector::create(1, 1),
//                             stride = Rcpp::IntegerVector::create(1, 1),
//                             block = Rcpp::IntegerVector::create(1, 1);
//         
//         if(threads.isNotNull()) {
//             if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()) {
//                 ithreads = Rcpp::as<int> (threads);
//             } else {
//                 ithreads = getDTthreads(0, true);
//             }
//         } else {
//             ithreads = getDTthreads(0, true);
//         }
//         
//         // Set minimum elements in block (mandatory : minimum = 2 * longest line)
//         if( dElementsBlock < dimensionSize * 2 ) {
//             minimumBlockSize = dimensionSize * 2;
//         } else {
//             minimumBlockSize = dElementsBlock;
//         }
//         
// 
//         Eigen::VectorXd Diagonal = as<Eigen::VectorXd>(Rcpp_getDiagonalfromMatrix(file, InOutDataset));
//         // Set diagonal in the new result matrix
//         Rcpp_setDiagonalMatrix( file, InOutDataset, Rcpp::wrap(Diagonal.cwiseInverse()) );
//         
//         Rcpp::Rcout<<"\n Matriu ORIGINAL : \n"<<GetCurrentBlock_hdf5(file, InOutDataset, 0, 0, idim0, idim1) <<"\n";
//         
//         if( idim0 == idim1) {
// // #pragma omp parallel for num_threads(getDTthreads(ithreads, true)) private(j) shared (L,TT,dimensionSize) schedule(static) if (j < dimensionSize)
//             for (int j = 1; j < dimensionSize; j++)
//             {
//                 Eigen::VectorXd vR = Eigen::VectorXd::Zero(j);
//                 Eigen::ArrayXd ar_j = GetCurrentBlock_hdf5(file, InOutDataset, j, 0, 1, j).transpose().array();
//                 
//                 Rcpp::Rcout<<"\n ========================================================== \n";
//                 Rcpp::Rcout<<"\n ---------------------------------------------------------- \n";
//                 Rcpp::Rcout<<"\n ar_j: "<<ar_j.transpose();
//                 Rcpp::Rcout<<"\n ---------------------------------------------------------- \n";
// 
//                 for (int i = 0; i < j ; i++) {
//                     Eigen::ArrayXd ar_i = GetCurrentBlock_hdf5(file, InOutDataset, 0, i, j, 1).array();
//                     Rcpp::Rcout<<"\n ar_i: "<<ar_i.transpose();
//                     
//                     vR(i) = ((ar_j.transpose() * ar_i.transpose()) * (-1)).sum() / Diagonal[j];
//                 }
//                 
//                 Rcpp::Rcout<<"\n ========================================================== \n";
//                 
//                 count[1] = j;
//                 offset[0] = j;
//                 
//                 write_HDF5_matrix_subset_v2( file, InOutDataset, offset, count, stride, block, Rcpp::wrap( vR.transpose() ) );
//             }
// 
//         } else {
//             throw std::range_error("non-conformable arguments");
//         }
//         
//         
//     } catch( FileIException& error ) { // catch failure caused by the H5File operations
//         InOutDataset->close();
//         file->close();
//         Rcpp::Rcout<<"c++ exception Inverse_of_Cholesky_decomposition_hdf5 (File IException)";
//         return void();
//     } catch( GroupIException & error ) { // catch failure caused by the DataSet operations
//         InOutDataset->close();
//         file->close();
//         Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5 (Group IException)";
//         return void();
//     } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
//         InOutDataset->close();
//         file->close();
//         Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5 (DataSet IException)";
//         return void();
//     } catch(std::exception& ex) {
//         InOutDataset->close();
//         file->close();
//         Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5" << ex.what();
//         return void();
//     }
// 
// return void();
//     
// }
//     






void Inverse_of_Cholesky_decomposition_hdf5(  H5File* file, DataSet* InOutDataset, int idim0, int idim1, double dElementsBlock, Rcpp::Nullable<int> threads = R_NilValue)
{
    
    try{
        
        unsigned int ithreads;
        int dimensionSize = idim0, 
            readedCols = 0,
            chunk = 1,
            colstoRead,
            minimumBlockSize;
        
        Rcpp::IntegerVector offset = Rcpp::IntegerVector::create(0, 0),
            count = Rcpp::IntegerVector::create(1, 1),
            stride = Rcpp::IntegerVector::create(1, 1),
            block = Rcpp::IntegerVector::create(1, 1);
        
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
        
        
        Eigen::VectorXd Diagonal = as<Eigen::VectorXd>(Rcpp_getDiagonalfromMatrix(file, InOutDataset));
        // Set diagonal in the new result matrix
        Rcpp_setDiagonalMatrix( file, InOutDataset, Rcpp::wrap(Diagonal.cwiseInverse()) );
        
        if( idim0 == idim1) {
            
            while ( readedCols < dimensionSize ) {
                
                colstoRead = ( -2 * readedCols - 1 + std::sqrt( pow(2*readedCols, 2) - 4 * readedCols + 8 * minimumBlockSize + 1) ) / 2;
                
                if( readedCols + colstoRead > idim0) { // Max size bigger than data to read ?
                    colstoRead = idim0 - readedCols;
                }
                
                offset[0] = readedCols;
                count[0] =  dimensionSize - offset[0];
                count[1] = readedCols + colstoRead;
                
                Eigen::MatrixXd verticalData = GetCurrentBlock_hdf5(file, InOutDataset, offset[0], offset[1], dimensionSize - offset[0], readedCols + colstoRead);
                
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
                    
#pragma omp parallel for num_threads(getDTthreads(ithreads, true)) shared (ar_j, j, verticalData, offset, colstoRead, vR) schedule(static) 
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
                
                write_HDF5_matrix_subset_v2( file, InOutDataset, offset, count, stride, block, Rcpp::wrap( verticalData ) );
                readedCols = readedCols + colstoRead; // Ho preparem perquè desprès necessitarem llegir a partir de la línea anterior
                    
            }
        } else {
            throw std::range_error("non-conformable arguments");
        }

    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        InOutDataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception Inverse_of_Cholesky_decomposition_hdf5 (File IException)";
        return void();
    } catch( GroupIException & error ) { // catch failure caused by the DataSet operations
        InOutDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5 (Group IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        InOutDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        InOutDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Inverse_of_Cholesky_decomposition_hdf5" << ex.what();
        return void();
    }
    
    return void();
    
}



void Inverse_Matrix_Cholesky_parallel(  H5File* file, DataSet* InOutDataset, int idim0, int idim1, double dElementsBlock, Rcpp::Nullable<int> threads = R_NilValue)
{

    try {

        int dimensionSize = idim0,
            readedCols = 0,
            colstoRead,
            minimumBlockSize;

        Eigen::VectorXd newDiag(idim0);
        unsigned int ithreads;

        Rcpp::IntegerVector offset = Rcpp::IntegerVector::create(0, 0),
                            count = Rcpp::IntegerVector::create(1, 1),
                            stride = Rcpp::IntegerVector::create(1, 1),
                            block = Rcpp::IntegerVector::create(1, 1);

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
                
                if(offset[0] + count[0]> idim0){
                    Rcpp::Rcout<<"\nWWWWHHHHHAAAAAATTTTT !!!!! offset[0] + count[0]> idim0 ------ PAAAAANNNNNIIIIICCCCCC!!!!\n";
                    Rcpp::Rcout<<"\nWWWWHHHHHAAAAAATTTTT !!!!! offset[0] + count[0]> idim0 ------ PAAAAANNNNNIIIIICCCCCC!!!!\n";
                }
                
                if(offset[1] + count[1]> idim1){
                    Rcpp::Rcout<<"\nWWWWHHHHHAAAAAATTTTT !!!!! offset[1] + count[1]> idim1 ------ PAAAAANNNNNIIIIICCCCCC!!!!\n";
                    Rcpp::Rcout<<"\nWWWWHHHHHAAAAAATTTTT !!!!! offset[1] + count[1]> idim1 ------ PAAAAANNNNNIIIIICCCCCC!!!!\n";
                }
                
                Eigen::MatrixXd verticalData = GetCurrentBlock_hdf5(file, InOutDataset, offset[0], offset[1], count[0], count[1]);
                
#pragma omp parallel for num_threads(getDTthreads(ithreads, true)) shared (verticalData, colstoRead, offset) schedule(static)
                for ( int i = 0; i < colstoRead + offset[0]; i++)   // Columnes
                {
                    int init;
                    int end;
                    newDiag( offset[0] + i) = verticalData.block(i, i + offset[0], verticalData.rows()-i, 1 ).array().pow(2).sum();
                    
                    if(offset[0] == 0) {
                        init = i + 1;
                    } else {
                        if(  i < offset[0]) {
                            init = 0;
                        } else {
                            if( i < offset[0] + colstoRead - 1) {
                                init = i - offset[0] + 1;
                            } else {
                                init = colstoRead; // force end
                            }
                        }
                    }
                    
                    // Rcpp::Rcout<<"\n Offset[0] val : "<<offset[0]<<" i val : "<<i<<" init val "<<init;
                    
                    for ( int j = init; j < colstoRead ; j++) {     // Files
                        
                        if( offset[0] + j < verticalData.cols()) {
                            verticalData(j,i) = (verticalData.block( j, i , verticalData.rows() - j, 1).array() * verticalData.block( j, j + offset[0],  verticalData.rows() - j, 1).array()).sum();    
                        }
                        
                    }
                }

                
                if(offset[0] + count[0]> idim0){
                    Rcpp::Rcout<<"\nWWWWHHHHHAAAAAATTTTT !!!!! offset[0] + count[0]> idim0 ------ PAAAAANNNNNIIIIICCCCCC FINAAAAAAL!!!!\n";
                    Rcpp::Rcout<<"\nWWWWHHHHHAAAAAATTTTT !!!!! offset[0] + count[0]> idim0 ------ PAAAAANNNNNIIIIICCCCCC FINAAAAAAL!!!!\n";
                }
                
                if(offset[1] + count[1]> idim1){
                    Rcpp::Rcout<<"\nWWWWHHHHHAAAAAATTTTT !!!!! offset[1] + count[1]> idim1 ------ PAAAAANNNNNIIIIICCCCCC!!!!\n";
                    Rcpp::Rcout<<"\nWWWWHHHHHAAAAAATTTTT !!!!! offset[1] + count[1]> idim1 ------ PAAAAANNNNNIIIIICCCCCC!!!!\n";
                }
                
                write_HDF5_matrix_subset_v2( file, InOutDataset, offset, count, stride, block, Rcpp::wrap( verticalData ) );
                readedCols = readedCols + colstoRead; // Ho preparem perquè desprès necessitarem llegir a partir de la línea anterior
            }

            Rcpp_setDiagonalMatrix(file, InOutDataset, Rcpp::wrap(newDiag));

        } else {
            throw std::range_error("non-conformable arguments");
        }

    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        InOutDataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception Inverse_Matrix_Cholesky_parallel (File IException)";
        return void();
    } catch( GroupIException & error ) { // catch failure caused by the DataSet operations
        InOutDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Inverse_Matrix_Cholesky_parallel (Group IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        InOutDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Inverse_Matrix_Cholesky_parallel (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        InOutDataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Inverse_Matrix_Cholesky_parallel" << ex.what();
        return void();
    }

    return void();
}
    
    


//' Compute inverse cholesky with hdf5 data files
//'
//' Compute inverse cholesky with datasets stored in hdf5 data files
//'
//' @param filename, character array with the name of an existin hdf5 data file containing the dataset to be modified
//' @param group, character array indicating the input group where the data set to be modified. 
//' @param datasets, character array indicating the input dataset to be modified
//' @param outdataset character array with output dataset name where we want to store results
//' @param outgroup optional, character array with output group name where we want to 
//' store results if not provided then results are stored in the same group as original dataset
//' @param threads optional parameter. Integer with numbers of threads to be used
//' @param force, optional boolean if true, previous results in same location inside 
//' hdf5 will be overwritten, by default force = false, data was not overwritten.
//' @param elementsBlock, optional integer defines de maximum number of elements to read from hdf5 data file in each block. 
//' By default this value is set to 10000. If matrix is bigger thant 5000x5000 then block is set to number of rows or columns x 2
//' @return Original hdf5 data file with Inverse of Cholesky
//' @examples
//' 
//' library(BigDataStatMeth)
//' library(rhdf5)
//' 
//' set.seed(1234)
//' A  <- matrix(sample.int(10, 10000, replace = TRUE), ncol = 100)
//' A <- crossprod(A)
//' 
//' # Create hdf5 data file with  data (Y)
//' bdCreate_hdf5_matrix_file("test_file2.hdf5", A, "data", "A", force = T)
//' 
//' # Get Inverse Cholesky
//' res <- bdInvCholesky_hdf5("test_file.hdf5", "data", "A", "results", "InverseA", force = T)
//' 
//' @export
// [[Rcpp::export]]
void bdInvCholesky_hdf5( std::string filename, std::string group, std::string dataset,
                         std::string  outdataset,
                         Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                         Rcpp::Nullable<bool> force = R_NilValue,
                         Rcpp::Nullable<int> threads = 2,
                         Rcpp::Nullable<double> elementsBlock = 1000000)
{
    
    H5File* file = nullptr;
    DataSet* pdataset = nullptr;
    DataSet* poutdataset_tmp = nullptr;
    double dElementsBlock;
    std::string strOutgroup, strIndataset, 
                strOutdataset, strOutdataset_tmp;
    
    int ithreads;
    bool bforce;
    
    try
    {
     
        // Get default values for Nullable variables
        if(force.isNull()) { bforce = false; } 
        else {  bforce = Rcpp::as<bool>(force); }
        
        if(outgroup.isNull()) { strOutgroup = group; } 
        else {   strOutgroup = Rcpp::as<std::string>(outgroup); }
        
        if(threads.isNull()) { ithreads = 2; } 
        else { ithreads = Rcpp::as<int>(threads); }
        
        if(elementsBlock.isNull()) { dElementsBlock = MAXELEMSINBLOCK; } 
        else { dElementsBlock = Rcpp::as<double>(elementsBlock); }
        
        // Rcpp::Rcout<<"\n Maxim elements block : "<<dElementsBlock<<"\n";
        
        strIndataset = group + "/" + dataset;
        strOutdataset = strOutgroup + "/" + outdataset;
        strOutdataset_tmp = "tmp/tmp_L";
        
        // Test if file exists dataset exists
        if(!ResFileExist(filename)) {
            Rcpp::Rcout<<"\nFile not exits, create file before get inverse of Cholesky\n";  
            return void();
        }
        
        file = new H5File( filename, H5F_ACC_RDWR );
        
        // Test if input dataset exists
        if(exists_HDF5_element_ptr(file, group)==0) {
            Rcpp::Rcout<<"\nGroup "<< group<<" not exits, create file and dataset before get inverse of Cholesky\n";
            file->close();
            return void();
        } else {
            
            if(!exists_HDF5_element_ptr(file, group + "/" + dataset)) {
                Rcpp::Rcout<<"\n Dataset "<< strIndataset <<" not exits, create dataset before get inverse of Cholesky\n";
                file->close();
                return void();
            }
        }
        
        // Test if output dataset exists and remove it if force = TRUE
        if( exists_HDF5_element_ptr(file, strOutdataset) && bforce == false) {
            Rcpp::Rcout<<"\n Dataset "<< strOutdataset <<"  also exists, please set force = TRUE to overwrite\n";
            file->close();
            return void();
        } else if(exists_HDF5_element_ptr(file, outdataset) && bforce == true) {
            remove_HDF5_element_ptr(file, strOutdataset); 
        }
        
        pdataset = new DataSet(file->openDataSet(strIndataset));
        
        // Real data set dimension
        IntegerVector dims_out = get_HDF5_dataset_size_ptr(pdataset);
        
        if(dims_out[0] == dims_out[1]) {
            
            // COM VEURE SI LA MATRIU ESTÀ DEFINIDA SEMI-POSITIVA ??
    /***
            Eigen::LLT<Eigen::MatrixXd> lltOfA(A);
            if(lltOfA.info() == Eigen::NumericalIssue) {
                throw std::runtime_error("Possibly non semi-positive definitie matrix!");
                
            } else {
    ***/
            // Test dataset
            if( !exists_HDF5_element_ptr(file, strOutdataset_tmp)) {
                create_HDF5_group_ptr(file, "tmp");
                create_HDF5_dataset_ptr(file, strOutdataset_tmp, dims_out[0], dims_out[1], "real"); 
                
            } else {
                remove_HDF5_element_ptr(file, strOutdataset_tmp);
                create_HDF5_dataset_ptr(file, strOutdataset_tmp, dims_out[0], dims_out[1], "real"); 
            }
            
                poutdataset_tmp = new DataSet(file->openDataSet(strOutdataset_tmp));
                
                Cholesky_decomposition_hdf5(file, pdataset, poutdataset_tmp, dims_out[0], dims_out[1], dElementsBlock, threads);
                // Rcpp::Rcout<<"\n =================>>>>>>>>>>>>< Acabat el primer pas !!! \n";
                pdataset->close();
                Inverse_of_Cholesky_decomposition_hdf5(  file, poutdataset_tmp, dims_out[0], dims_out[1], dElementsBlock, threads); // Resultats emmagatzemats Triangular inferior (menys diagonal)
                // // Rcpp::Rcout<<"\n =================>>>>>>>>>>>>< Acabat el segon pas !!! \n";
                Inverse_Matrix_Cholesky_parallel( file, poutdataset_tmp, dims_out[0], dims_out[1], dElementsBlock, threads); // Resultats emmagatzemats Triangular superior (juntament amb diagonal)
                //  Rcpp::Rcout<<"\n =================>>>>>>>>>>>>< Acabat el tercer pas !!! \n";
                
                
            
        } else {
            pdataset->close();
            file->close();
            Rcpp::Rcout<<"\n Can't get inverse matrix using Cholesky decomposition \n";
            return void();
        }
        
        
    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        // Rcpp::Rcout<<"\nPASSA PER FILE\n";
        // poutdataset->close();
        poutdataset_tmp->close();
        pdataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception bdInvCholesky_hdf5 (File IException)";
        return void();
    } catch( GroupIException & error ) { // catch failure caused by the DataSet operations
        // Rcpp::Rcout<<"\nPASSA PER GROUP\n";
        // poutdataset->close();
        poutdataset_tmp->close();
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdInvCholesky_hdf5 (Group IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        // Rcpp::Rcout<<"\nPASSA PER DATASET\n";
        // poutdataset->close();
        poutdataset_tmp->close();
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdInvCholesky_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        // Rcpp::Rcout<<"\nPASSA PER EXCEPTION\n";
        // poutdataset->close();
        poutdataset_tmp->close();
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdInvCholesky_hdf5" << ex.what();
        return void();
    }
    
    poutdataset_tmp->close();
    file->close();
    return void();
}






/*** R

library(BigDataStatMeth)
library(rhdf5)

#..# setwd("C:/tmp_test/")
setwd("/Volumes/XtraSpace/PhD_Test/BigDataStatMeth")

set.seed(1234)
A  <- matrix(sample.int(10, 1000000, replace = TRUE), ncol = 100)
A <- crossprod(A)
E <- BigDataStatMeth:::inversechol_par(A,threads = 2)
# 
# # Test on-memory matrix
# inversechol_par(A, 2)

# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("test_file22.hdf5", A, "data", "A", force = TRUE)
# bdCreate_hdf5_matrix_file("test_file23.hdf5", A, "data", "A", force = TRUE)

# Get Inverse Cholesky
# res <- bdInvCholesky_hdf5("test_file22.hdf5", "data", "A", "results", "InverseA", force = T)
res <- bdInvCholesky_hdf5("test_file22.hdf5", "data", "A", "results", "InverseA",elementsBlock = 10)

#..# E <- BigDataStatMeth:::inversechol_par(A,threads = 2)


microbenchmark::microbenchmark( T <- solve(A),
                                # res <- bdInvCholesky_hdf5("test_file22.hdf5", "data", "A", "results", "InverseA",elementsBlock = 20),
                                res <- bdInvCholesky_hdf5("test_file22.hdf5", "data", "A", "results", "InverseA", elementsBlock = 100000),
                                # res <- bdInvCholesky_hdf5("test_file23.hdf5", "data", "A", "results", "InverseA", elementsBlock = 10),
                                F <- BigDataStatMeth:::inversechol_par(A),
                                G2 <- BigDataStatMeth:::inversechol_par(A,threads = 4),
                                times = 5
)

microbenchmark::microbenchmark( res <- bdInvCholesky_hdf5("test_file22.hdf5", "data", "A", "results", "InverseA", elementsBlock = 10),
                                times = 1
)


devtools::reload(pkgload::inst("BigDataStatMeth"))
library(BigDataStatMeth)
setwd("/Volumes/XtraSpace/PhD_Test/BigDataStatMeth")
set.seed(1234)
A  <- matrix(sample.int(10, 10000, replace = TRUE), ncol = 100)
A <- crossprod(A)


library(BigDataStatMeth)

setwd("/Volumes/XtraSpace/PhD_Test/BigDataStatMeth")
# Get Inverse Cholesky
res <- bdInvCholesky_hdf5("test_file22.hdf5", "data", "A", "results", "InverseA", force = T,elementsBlock = 20)

*/
