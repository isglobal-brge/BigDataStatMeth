#include "include/hdf5_getDiagonalMatrix.h"


//' Write diagonal matrix
//'
//' Write diagonal matrix to an existing dataset inside hdf5
//'
//' @param filename, character array with the name of an existin hdf5 data file containing the dataset to be modified
//' @param group, character array indicating the input group where the data set to be modified. 
//' @param datasets, character array indicating the input dataset to be modified
//' @return Numeric vector with all diagonal elements from hdf5 dataset
//' @examples
//' 
//' library(BigDataStatMeth)
//' library(rhdf5)
//' 
//' X <- matrix(rnorm(150), 10, 10)
//' diag(X) <- 0.5
//' # Create hdf5 data file with  data (Y)
//' bdCreate_hdf5_matrix_file("test_file2.hdf5", X, "data", "X", force = T)
//' # Update diagonal
//' diagonal <- bdgetDiagonal_hdf5("test_file.hdf5", "data", "X")
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdgetDiagonal_hdf5( std::string filename, std::string group, std::string dataset)
{
    
    H5File* file = nullptr;
    DataSet* pdataset = nullptr;
    Rcpp::NumericVector intNewDiagonal;

    try
    {
        
        Rcpp::IntegerVector offset = Rcpp::IntegerVector::create(0, 0);
        Rcpp::IntegerVector count = Rcpp::IntegerVector::create(1, 1);
        Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1, 1);
        Rcpp::IntegerVector block = Rcpp::IntegerVector::create(1, 1);
        
        std::string strDataset = group + "/" + dataset; 
        
        // Test file
        if( ResFileExist_filestream(filename) ) {
            file = new H5File( filename, H5F_ACC_RDWR ); 
        } else {
            Rcpp::Rcout<<"\nFile not exits, create file before bind matrices";
            return(Rcpp::wrap(0));
        }
        
        if( exists_HDF5_element_ptr(file, strDataset ) ) {
            pdataset = new DataSet(file->openDataSet(strDataset));
        } else {
            file -> close();
            ::Rf_error( "c++ exception in bdgetDiagonal_hdf5 (Dataset does not exist!)" );
            return(Rcpp::wrap(0));
        }
        
        Rcpp::IntegerVector dims_out = get_HDF5_dataset_size(*pdataset);
        Rcpp::NumericVector data(1);
        
        // H5Sselect_elements()
        for(int i=0; i < dims_out[0]; i++) {
            offset[0] = i; offset[1] = i;
            read_HDF5_matrix_subset(file, pdataset, offset, count, stride, block, REAL(data));
            intNewDiagonal.push_back(data[0]);
        }
        
    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception bdgetDiagonal_hdf5 (File IException)";
        return(Rcpp::wrap(0));
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdgetDiagonal_hdf5 (DataSet IException)";
        return(Rcpp::wrap(0));
    } catch(std::exception& ex) {
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdgetDiagonal_hdf5" << ex.what();
        return(Rcpp::wrap(0));
    }
    
    pdataset->close();
    file->close();
    return(intNewDiagonal);
    
}



/***R

library(BigDataStatMeth)
library(rhdf5)

setwd("C:/tmp_test/")

# devtools::reload(pkgload::inst("BigDataStatMeth"))

# Prepare data and functions
X <- matrix(rnorm(150), 10, 10)
diag(X) <- 0.5

# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("test_file2.hdf5", X, "data", "X", force = T)

# Update diagonal
diagonal <- bdgetDiagonal_hdf5("test_file.hdf5", "data", "X")

*/
