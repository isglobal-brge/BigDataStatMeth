#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixTriangular.hpp"
#include "hdf5Utilities/hdf5Utilities.hpp"


//' Write Upper/Lower triangular matrix
//'
//' Write diagonal matrix to an existing dataset inside hdf5
//'
//' @param filename, character array with the name of an existin hdf5 data file containing the dataset to be modified
//' @param group, character array indicating the input group where the data set to be modified. 
//' @param dataset, character array indicating the input dataset to be modified
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
//' X <- matrix(rnorm(100), 10, 10)
//' X.1 <- X
//' X[lower.tri(X)] <- 0
//' # Create hdf5 data file with  data (Y)
//' bdCreate_hdf5_matrix_file("test_file.hdf5", X, "data", "X", force = TRUE)
//' # Update Lower triangular matrix in hdf5
//' bdWriteOppsiteTriangularMatrix_hdf5(filename = "test_file.hdf5", 
//'         group = "data", dataset = "X", copytolower = TRUE, elementsBlock = 10)
//' 
//' X <- X.1
//' X[upper.tri(X)] <- 0
//' # CAdd matrix data to a file
//' bdAdd_hdf5_matrix(X, "test_file.hdf5", "data", "Y", force = TRUE)
//' # Update Upper triangular matrix in hdf5
//' bdWriteOppsiteTriangularMatrix_hdf5(filename = "test_file.hdf5", 
//'         group = "data", dataset = "Y", copytolower = FALSE, elementsBlock = 10)
//' 
//' @export
// [[Rcpp::export]]
void bdWriteOppsiteTriangularMatrix_hdf5(std::string filename, 
                                      std::string group, std::string dataset, 
                                      Rcpp::Nullable<bool> copytolower = R_NilValue,
                                      Rcpp::Nullable<long> elementsBlock = 1000000)
 {
     
     Rcpp::NumericVector intNewDiagonal;
     bool blower;
     long dElementsBlock;
     
     try
     {
         
         // Get default values for Nullable variables
         if(copytolower.isNull()) { blower = false; } 
         else {  blower = Rcpp::as<bool>(copytolower); }
         
         if(elementsBlock.isNull()) { dElementsBlock = BigDataStatMeth::MAXELEMSINBLOCK; } 
         else { dElementsBlock = Rcpp::as<long>(elementsBlock); }
         
         std::string strDataset = group + "/" + dataset;
         
         
         BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
         dsA->openDataset();
         
         if(dsA->nrows() != dsA->ncols()) {
             Rcpp::Rcout<<"\nCan not write opposite triangular matrix - Non squuare matrix";
             delete dsA;
             return void();   
         }
         
         if( blower == false ) {
             setLowerTriangularMatrix( dsA, dElementsBlock);
         } else {
             setUpperTriangularMatrix( dsA, dElementsBlock);
         }
         
         delete dsA;
         
     }
     catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
         Rcpp::Rcout<<"c++ exception bdWriteOppsiteTriangularMatrix_hdf5 (File IException)";
         return void();
     } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
         Rcpp::Rcout << "c++ exception bdWriteOppsiteTriangularMatrix_hdf5 (DataSet IException)";
         return void();
     } catch(std::exception& ex) {
         Rcpp::Rcout << "c++ exception bdWriteOppsiteTriangularMatrix_hdf5" << ex.what();
         return void();
     }
     
     Rcpp::Rcout<<dataset<<" Triangular Matrix has been mirrored\n";
     return void();
 }
