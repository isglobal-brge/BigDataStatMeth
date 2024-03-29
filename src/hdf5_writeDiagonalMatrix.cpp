#include "include/hdf5_writeDiagonalMatrix.h"


void Rcpp_setDiagonalMatrix( H5File* file, DataSet* pdataset, Rcpp::NumericVector intNewDiagonal)
{
    
    // Rcpp::NumericVector intNewDiagonal;
    Rcpp::IntegerVector count = Rcpp::IntegerVector::create(1, 1);
    Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1, 1);
    Rcpp::IntegerVector block = Rcpp::IntegerVector::create(1, 1);
    
    Rcpp::IntegerVector dims_out = get_HDF5_dataset_size(*pdataset);

    // H5Sselect_elements()
    for(int i=0; i < intNewDiagonal.size(); i++) {
        Rcpp::IntegerVector offset = Rcpp::IntegerVector::create(i, i);
        write_HDF5_matrix_subset_v2(file, pdataset, offset, count, stride, block, wrap(intNewDiagonal(i)) );
    }
    
    return void();
}





//' Write diagonal matrix
//'
//' Write diagonal matrix to an existing dataset inside hdf5
//'
//' @param diagonal, numeric vector with diagonal elements to be written in existing dataset. 
//' @param filename, character array with the name of an existin hdf5 data file containing the dataset to be modified
//' @param group, character array indicating the input group where the data set to be modified. 
//' @param dataset, character array indicating the input dataset to be modified
//' @return Original hdf5 dataset with new diagonal elements
//' @examples
//' library(BigDataStatMeth)
//' library(rhdf5)
//' 
//' # Prepare data and functions
//' X <- matrix(rnorm(100), 10, 10)
//' diagonal <- c(1,2,3,4,5,6,7, 8, 9, 10)
//' 
//' # Create hdf5 data file with  data (Y)
//' bdCreate_hdf5_matrix_file("test_file.hdf5", X, "data", "X", force = TRUE)
//' 
//' # Update diagonal
//' bdWriteDiagonal_hdf5(diagonal, "test_file.hdf5", "data", "X")
//' 
//' @export
// [[Rcpp::export]]
void bdWriteDiagonal_hdf5( Rcpp::RObject diagonal, std::string filename, std::string group, std::string dataset)
{
    
    H5File* file = nullptr;
    DataSet* pdataset = nullptr;
    Rcpp::NumericVector intNewDiagonal;

    try
    {
        
        // Rcpp::IntegerVector count = Rcpp::IntegerVector::create(1, 1);
        // Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1, 1);
        // Rcpp::IntegerVector block = Rcpp::IntegerVector::create(1, 1);
        
        std::string strDataset = group + "/" + dataset;
        
        if( !is<Rcpp::IntegerVector>(diagonal) || !is<Rcpp::NumericVector>(diagonal) ) {
            intNewDiagonal = Rcpp::as<Rcpp::NumericVector>(diagonal);
        } else {
            Rcpp::Rcout<<"\n Diagonal vector isn't a Numeric vector";
            return void();
        }
        
        // Test file
        if( ResFileExist_filestream(filename) ) {
            file = new H5File( filename, H5F_ACC_RDWR ); 
        } else {
            Rcpp::Rcout<<"\nFile not exits, create file before bind matrices";
            return void();
        }
        
        if( exists_HDF5_element_ptr(file, strDataset ) ) {
            pdataset = new DataSet(file->openDataSet(strDataset));
        } else {
            file -> close();
            ::Rf_error( "c++ exception in bdWriteDiagonal_hdf5 (Dataset does not exist!)" );
            return void();
        }
        
        Rcpp_setDiagonalMatrix( file, pdataset, intNewDiagonal);
            
        // for(int i=0; i < intNewDiagonal.size(); i++) {
        //     Rcpp::IntegerVector offset = Rcpp::IntegerVector::create(i, i);
        //     write_HDF5_matrix_subset_v2(file, pdataset, offset, count, stride, block, wrap(intNewDiagonal(i)) );
        // }
        
    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception bdWriteDiagonal_hdf5 (File IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdWriteDiagonal_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdWriteDiagonal_hdf5" << ex.what();
        return void();
    }
    
    pdataset->close();
    file->close();
    Rcpp::Rcout<<dataset<<" diagonal has been overwritten\n";
    return void();
}



/***R

library(BigDataStatMeth)
library(rhdf5)

# setwd("C:/tmp_test/")
setwd("/Volumes/XtraSpace/PhD_Test/BigDataStatMeth")

# devtools::reload(pkgload::inst("BigDataStatMeth"))

# Prepare data and functions
X <- matrix(rnorm(100), 10, 10)
Y <- matrix(rnorm(250), 50, 5)
diagonal <- c(1,2,3,4,5,6,7, 8, 9, 10)

# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("test_file.hdf5", X, "data", "X", force = TRUE)

# Update diagonal
bdWriteDiagonal_hdf5(diagonal, "test_file.hdf5", "data", "X")

*/
