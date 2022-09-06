#include "include/hdf5_writeOppositeTriangular.h"


void Rcpp_setUpperTriangularMatrix( H5File* file, DataSet* pdataset, int dimensionSize, long dElementsBlock)
{
    
    try {
        
        Rcpp::IntegerVector count = Rcpp::IntegerVector::create(1, 1),
                            stride = Rcpp::IntegerVector::create(1, 1),
                            block = Rcpp::IntegerVector::create(1, 1),
                            offset = Rcpp::IntegerVector::create(0, 0);
        
        int readedRows = 0,
            rowstoRead,
            minimumBlockSize;
        
        if( dElementsBlock < dimensionSize * 2 ) {
            minimumBlockSize = dimensionSize * 2;
        } else {
            minimumBlockSize = dElementsBlock;
        }
        
        while ( readedRows < dimensionSize ) {
            
            rowstoRead = ( -2 * readedRows - 1 + std::sqrt( pow(2*readedRows, 2) - 4 * readedRows + 8 * minimumBlockSize + 1) ) / 2;
            
            if( readedRows + rowstoRead > dimensionSize) { // Max size bigger than data to read ?
                rowstoRead = dimensionSize - readedRows;
            }
            
            if( readedRows == 0) { // Read complete cols
                count[0] = dimensionSize;
                count[1] = rowstoRead;
            } else {
                count[1] = rowstoRead;
                count[0] = dimensionSize - readedRows;
            }
            
            offset[1] = readedRows;
            offset[0] = offset[1];
            
            Eigen::MatrixXd A = GetCurrentBlock_hdf5(file, pdataset, offset[0], offset[1], count[0], count[1]);
            
            // #pragma omp parallel for num_threads(getDTthreads(ithreads, true)) private(sum) shared (A,L,j) schedule(static) if (j < readedRows - chunk)
            for ( int i = 0; i < rowstoRead; i++ ) {
                Rcpp::IntegerVector tcount = Rcpp::IntegerVector::create(1, dimensionSize-readedRows-i-1);
                Rcpp::IntegerVector toffset = Rcpp::IntegerVector::create( readedRows+i, readedRows+i+1);
                
                write_HDF5_matrix_subset_v2( file, pdataset, toffset, tcount, stride, block, Rcpp::wrap( A.block(i+1, i, dimensionSize-readedRows-i-1, 1 ).transpose() ) );
            }
            
            readedRows = readedRows + rowstoRead; 
        }
        
    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception Rcpp_setUpperTriangularMatrix (File IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Rcpp_setUpperTriangularMatrix (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Rcpp_setUpperTriangularMatrix" << ex.what();
        return void();
    }
    
    return void();
}


void Rcpp_setLowerTriangularMatrix( H5File* file, DataSet* pdataset, int dimensionSize, long dElementsBlock)
{
    
    try {
        
        Rcpp::IntegerVector count = Rcpp::IntegerVector::create(1, 1),
            stride = Rcpp::IntegerVector::create(1, 1),
            block = Rcpp::IntegerVector::create(1, 1),
            offset = Rcpp::IntegerVector::create(0, 0);
        
        int readedRows = 0,
            rowstoRead,
            minimumBlockSize;
        
        if( dElementsBlock < dimensionSize * 2 ) {
            minimumBlockSize = dimensionSize * 2;
        } else {
            minimumBlockSize = dElementsBlock;
        }
        
        while ( readedRows < dimensionSize ) {
            
            rowstoRead = ( -2 * readedRows - 1 + std::sqrt( pow(2*readedRows, 2) - 4 * readedRows + 8 * minimumBlockSize + 1) ) / 2;
            
            if( readedRows + rowstoRead > dimensionSize) { // Max size bigger than data to read ?
                rowstoRead = dimensionSize - readedRows;
            }
            
            if( readedRows == 0) { // Read complete cols
                count[0] = rowstoRead;
                count[1] = dimensionSize;
            } else {
                count[1] = dimensionSize - readedRows;
                count[0] = rowstoRead;
            }
            
            offset[1] = readedRows;
            offset[0] = offset[1];
            
            Eigen::MatrixXd A = GetCurrentBlock_hdf5(file, pdataset, offset[0], offset[1], count[0], count[1]);
            
            // #pragma omp parallel for num_threads(getDTthreads(ithreads, true)) private(sum) shared (A,L,j) schedule(static) if (j < readedRows - chunk)
            for ( int i = 0; i < rowstoRead; i++ ) {
                Rcpp::IntegerVector tcount = Rcpp::IntegerVector::create( dimensionSize-readedRows-i-1, 1);
                Rcpp::IntegerVector toffset = Rcpp::IntegerVector::create( readedRows+i+1, readedRows+i);
                
                write_HDF5_matrix_subset_v2( file, pdataset, toffset, tcount, stride, block, Rcpp::wrap( A.block(i, i+1, 1, dimensionSize-readedRows-i-1 ).transpose() ) );
            }
            
            readedRows = readedRows + rowstoRead; 
        }
        
    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception Rcpp_setLowerTriangularMatrix (File IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Rcpp_setLowerTriangularMatrix (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception Rcpp_setLowerTriangularMatrix" << ex.what();
        return void();
    }
    
    return void();
    
}


//' Write Upper/Lower triangular matrix
//'
//' Write diagonal matrix to an existing dataset inside hdf5
//'
//' @param filename, character array with the name of an existin hdf5 data file containing the dataset to be modified
//' @param group, character array indicating the input group where the data set to be modified. 
//' @param datasets, character array indicating the input dataset to be modified
//' @param copytolower, boolean with default value = false. If true, sets lower 
//' triangular matrix using upper triangular matrix, if lower=false (default 
//' value) sets upper triangular matrix using lower triangular matrix.
//' @param elementsBlock, optional integer defines de maximum number of elements 
//' to read from hdf5 data file in each block. By default this value is set 
//' to 10000. If matrix is bigger thant 5000x5000 then block is set to number 
//' of rows or columns x 2
//' @examples
//' library(BigDataStatMeth)
//' 
//' # Prepare data and functions
//' X <- matrix(rnorm(150), 10, 10)
//' X.1 <- X
//' X[lower.tri(X)] <- 0
//' # Create hdf5 data file with  data (Y)
//' bdCreate_hdf5_matrix_file("test_file.hdf5", X, "data", "X", force = T)
//' # Update Lower triangular matrix in hdf5
//' bdWriteOppsiteTriangularMatrix_hdf5(filename = "test_file.hdf5", 
//'         group = "data", dataset = "X", copytolower = T, elementsBlock = 10)
//' 
//' X <- X.1
//' X[upper.tri(X)] <- 0
//' # CAdd matrix data to a file
//' bdAdd_hdf5_matrix(X, "test_file.hdf5", "data", "Y", force = T )
//' # Update Upper triangular matrix in hdf5
//' bdWriteOppsiteTriangularMatrix_hdf5(filename = "test_file.hdf5", 
//'         group = "data", dataset = "Y", copytolower = F, elementsBlock = 10)
//' 
//' @export
// [[Rcpp::export]]
void bdWriteOppsiteTriangularMatrix_hdf5(std::string filename, 
                                std::string group, std::string dataset, 
                                Rcpp::Nullable<bool> copytolower = R_NilValue,
                                Rcpp::Nullable<long> elementsBlock = 1000000)
{
    
    H5File* file = nullptr;
    DataSet* pdataset = nullptr;
    Rcpp::NumericVector intNewDiagonal;
    bool blower;
    long dElementsBlock;
    
    try
    {
        IntegerVector dims_out;
        
        
        // Get default values for Nullable variables
        if(copytolower.isNull()) { blower = false; } 
        else {  blower = Rcpp::as<bool>(copytolower); }
        
        if(elementsBlock.isNull()) { dElementsBlock = MAXELEMSINBLOCK; } 
        else { dElementsBlock = Rcpp::as<long>(elementsBlock); }
        
        std::string strDataset = group + "/" + dataset;
        

        // Test file
        if( ResFileExist_filestream(filename) ) {
            file = new H5File( filename, H5F_ACC_RDWR ); 
        } else {
            Rcpp::Rcout<<"\nFile not exits, create file before copy triangular matrix";
            return void();
        }
        
        if( exists_HDF5_element_ptr(file, strDataset ) ) {
            pdataset = new DataSet(file->openDataSet(strDataset));
            // Real data set dimension
            dims_out = get_HDF5_dataset_size_ptr(pdataset);
            if(dims_out[0] != dims_out[1]) {
                pdataset->close();
                file->close();
                Rcpp::Rcout<<"\nCan not write opposite triangular matrix - Non squuare matrix";
                return void();
            }

        } else {
            file -> close();
            ::Rf_error( "c++ exception in bdWriteOppsiteTriangularMatrix_hdf5 (Dataset does not exist!)" );
            return void();
        }
        
        
        
        if(blower == false) {
            Rcpp_setLowerTriangularMatrix( file, pdataset, dims_out[0], dElementsBlock);
        } else {
            Rcpp_setUpperTriangularMatrix( file, pdataset, dims_out[0], dElementsBlock);
        }
        

    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception bdWriteOppsiteTriangularMatrix_hdf5 (File IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdWriteOppsiteTriangularMatrix_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdWriteOppsiteTriangularMatrix_hdf5" << ex.what();
        return void();
    }
    
    pdataset->close();
    file->close();
    Rcpp::Rcout<<dataset<<" Triangular Matrix has been mirrored\n";
    return void();
}



/***R

library(BigDataStatMeth)
library(rhdf5)

# setwd("C:/tmp_test/")
setwd("/Volumes/XtraSpace/PhD_Test/BigDataStatMeth")

# devtools::reload(pkgload::inst("BigDataStatMeth"))

# Prepare data and functions
X <- matrix(rnorm(150), 10, 10)

X.1 <- X
X[lower.tri(X)] <- 0

# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("test_file.hdf5", X, "data", "X", force = T)
# Update diagonal
bdWriteOppsiteTriangularMatrix_hdf5(filename = "test_file.hdf5", group = "data", dataset = "X", copytolower = T, elementsBlock = 10)


X <- X.1
X[upper.tri(X)] <- 0

# Create hdf5 data file with  data (Y)
bdAdd_hdf5_matrix(X, "test_file.hdf5", "data", "Y", force = T )
# Update diagonal
bdWriteOppsiteTriangularMatrix_hdf5(filename = "test_file.hdf5", group = "data", dataset = "Y", copytolower = F, elementsBlock = 10)

*/
