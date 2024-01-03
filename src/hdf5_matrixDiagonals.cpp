#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixDiagonal.hpp"




//' Get diagonal matrix
//'
//' Gry diagonal matrix from an existing dataset inside hdf5
//'
//' @param filename, character array with the name of an existin hdf5 data file containing the dataset to be modified
//' @param group, character array indicating the input group where dataset is stored
//' @param dataset, character array indicating the input dataset name
//' @return Numeric vector with all diagonal elements from hdf5 dataset
//' @examples
//' 
//' library(BigDataStatMeth)
//' library(rhdf5)
//' 
//' X <- matrix(rnorm(100), 10, 10)
//' diag(X) <- 0.5
//' # Create hdf5 data file with  data (Y)
//' bdCreate_hdf5_matrix_file("test_file2.hdf5", X, "data", "X", force = TRUE)
//' # Update diagonal
//' diagonal <- bdgetDiagonal_hdf5("test_file.hdf5", "data", "X")
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdgetDiagonal_hdf5( std::string filename, std::string group, std::string dataset)
{
 
     Rcpp::NumericVector intNewDiagonal;
     
     try
     {
         BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
         dsA->openDataset();
         
         intNewDiagonal = getDiagonalfromMatrix(dsA);
         
         delete dsA;
         
     } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
         Rcpp::Rcout<<"c++ exception bdgetDiagonal_hdf5 (File IException)";
         return(Rcpp::wrap(0));
     } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
         Rcpp::Rcout << "c++ exception bdgetDiagonal_hdf5 (DataSet IException)";
         return(Rcpp::wrap(0));
     } catch(std::exception& ex) {
         Rcpp::Rcout << "c++ exception bdgetDiagonal_hdf5" << ex.what();
         return(Rcpp::wrap(0));
     }
     
     return(intNewDiagonal);
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
     
     Rcpp::NumericVector intNewDiagonal;
     
     try
     {
         
         std::string strDataset = group + "/" + dataset;
         
         if( !Rcpp::is<Rcpp::IntegerVector>(diagonal) || !Rcpp::is<Rcpp::NumericVector>(diagonal) ) {
             intNewDiagonal = Rcpp::as<Rcpp::NumericVector>(diagonal);
         } else {
             Rcpp::Rcout<<"\n Diagonal vector isn't a Numeric vector";
             return void();
         }
         
         BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
         dsA->openDataset();
         
         setDiagonalMatrix( dsA, intNewDiagonal);
         
         delete dsA;
         
     } catch(std::exception& ex) {
         Rcpp::Rcout << "c++ exception bdWriteDiagonal_hdf5" << ex.what();
         return void();
     }
     
     Rcpp::Rcout<<dataset<<" diagonal has been overwritten\n";
     return void();
 }
