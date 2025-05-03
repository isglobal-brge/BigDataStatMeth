#include <BigDataStatMeth.hpp>
// #include "hdf5Algebra/vectormatrix.hpp"



//' Apply vector calculus to a dataset in hdf5 file
//' 
//' This function applies a calculus with a vector to a matrix. Multiplies, 
//' sums, substract or divide each matrix row/column from a hdf5 dataset 
//' using a vector
//' 
//' @param filename string file name where dataset to apply weights is located
//' @param group string with the path inside the hdf5 data file where matrix 
//' is located
//' @param dataset string with the matrix name
//' @param vectorgroup string with the path inside the hdf5 data file where 
//' vector is located
//' @param vectordataset string with the vector name
//' @param outdataset character array with output dataset name where we want to 
//' store results
//' @param outgroup optional, character array with output group name where we 
//' want to store results if not provided then results are stored in the same 
//' group as original dataset
//' @param func, Character array, function to be applyed : 
//'"+" : to sum a vector to a matrix dataset by columns or rows
//'"-" : to substract a vector to a matrix dataset by columns or rows
//'"*" : to multiply a vector to a matrix dataset by columns or rows
//'"/" : to divide a vector to a matrix dataset by columns or rows
//' @param byrows logical (default = FALSE). By default weights are applied by 
//' columns but if byrows=TRUE then weights are applied by rows 
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel 
//' computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent 
//' threads in parallelization if threads is null then threads =  maximum 
//' number of threads available
//' @param overwrite, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
//' @return file with weighted dataset
//' @examples
//'library(BigDataStatMeth)
//'    
//'# Prepare data and functions
//'set.seed(123)
//'Y <- matrix(rnorm(100), 10, 10)
//'X <- matrix(rnorm(10), 10, 1)
//'        
//'# Create hdf5 data file with  data (Y)
//'bdCreate_hdf5_matrix("test.hdf5", Y, "data", "Y", overwriteFile = TRUE, 
//'                         overwriteDataset = FALSE, 
//'                         unlimited = FALSE)
//'bdCreate_hdf5_matrix("test.hdf5",  X, "data", "X", overwriteFile = FALSE, 
//'                         overwriteDataset = FALSE, 
//'                         unlimited = FALSE)
//'            
//'bdcomputeMatrixVector_hdf5("test.hdf5", 
//'                           group = "data", dataset = "Y",
//'                           vectorgroup = "data", vectordataset = "X", 
//'                           outdataset = "ProdComputed", 
//'                           func = "*",
//'                           byrows = TRUE, overwrite = TRUE)
//'    
//'bdcomputeMatrixVector_hdf5("test.hdf5", 
//'                           group = "data", dataset = "Y",
//'                           vectorgroup = "data", vectordataset = "X", 
//'                           outdataset = "SumComputed", 
//'                           func = "-",
//'                           byrows = TRUE, overwrite = TRUE)
//'    
//'bdcomputeMatrixVector_hdf5("test.hdf5", 
//'                           group = "data", dataset = "Y",
//'                           vectorgroup = "data", vectordataset = "X", 
//'                           outdataset = "SubsComputed", 
//'                           func = "-",
//'                           byrows = FALSE, overwrite = TRUE)
//'                           
//' # Remove file (used as example)
//' if (file.exists("test.hdf5")) {
//'   file.remove("test.hdf5")
//' }
//' 
//' @export
// [[Rcpp::export]]
 void bdcomputeMatrixVector_hdf5( std::string filename, std::string group, 
                                  std::string dataset,
                                  std::string vectorgroup, std::string vectordataset,
                                  std::string outdataset, 
                                  std::string func,
                                  Rcpp::Nullable<std::string> outgroup = R_NilValue,
                                  Rcpp::Nullable<bool> byrows = R_NilValue,
                                  Rcpp::Nullable<bool> paral = R_NilValue,
                                  Rcpp::Nullable<int> threads = R_NilValue,
                                  Rcpp::Nullable<int> overwrite  = false)
 {
     
     bool bbyrows, bparal, bforce;
     // int blocksize;
     std::string strgroupout;
     bool bError=false;
     
     Rcpp::NumericVector oper = {0, 1, 2, 3};
     oper.names() = Rcpp::CharacterVector({"+", "-", "*", "/"});
     BigDataStatMeth::hdf5Dataset* dsA = nullptr;
     BigDataStatMeth::hdf5Dataset* dsB = nullptr;
     BigDataStatMeth::hdf5Dataset* dsC = nullptr;
     
     try{
         
         if( byrows.isNull()) { bbyrows = false; } 
         else { bbyrows = Rcpp::as<bool> (byrows); }
         
         if( overwrite.isNull()) { bforce = true; } 
         else { bforce = Rcpp::as<bool> (overwrite); }
         
         if( outgroup.isNull()) { strgroupout = group; } 
         else { strgroupout = Rcpp::as<std::string> (outgroup); }
         
         if (paral.isNull()) { bparal = false; } 
         else { bparal = Rcpp::as<bool> (paral); }
         
         // Function exists?
         if( oper(oper.findName(func)) != 0 && oper(oper.findName(func)) != 1 &&
             oper(oper.findName(func)) != 2 && oper(oper.findName(func)) != 3)  {
             Rcpp::Rcout<<"Function does not exists, please use one of the following : '+', '-', '*', '/' ";
             return void();
         } 

         dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
         dsA->openDataset();
         
         dsB = new BigDataStatMeth::hdf5Dataset(filename, vectorgroup, vectordataset, false);
         dsB->openDataset();
         
         if( dsA->getDatasetptr() != nullptr && dsB->getDatasetptr() != nullptr)  {
             
             // Check Vector dataset and matrix dataset A should be matrix and B should be the vector
             if( dsA->nrows()==1 || dsA->ncols()==1) {
                 Rcpp::Rcout<<"\ndataset is not a matrix";
                 bError = true;
             }
             
             if( dsB->nrows()!=1 && dsB->ncols()!=1) {
                 Rcpp::Rcout<<"\nvectordataset is not a vector";
                 bError = true;
             }
             
             if( bError == false ) {
                 dsC = new BigDataStatMeth::hdf5Dataset(filename, strgroupout, outdataset, bforce);
                 
                 if( (bbyrows == false && dsB->ncols() != dsA->ncols()) || dsB->nrows() != 1 ) {
                     Rcpp::Rcout<<"\nNon-conformable dimensions or vector dataset is not a vector";
                     Rcpp::Rcout<<"\nCalculus has not been computed\n";
                 } else if( bbyrows == true && dsB->ncols() != dsA->nrows()) {
                     Rcpp::Rcout<<"\nNon-conformable dimensions - byrows = true";
                     Rcpp::Rcout<<"\nCalculus has not been computed\n";
                 } else {
                     
                     dsC = hdf5_matrixVector_calculus( dsA, dsB, dsC, oper(oper.findName(func)), bbyrows, bparal, threads);
                 }
                 
                 delete dsC; dsC = nullptr;
             }    
             
         }
         
         delete dsA; dsA = nullptr;
         delete dsB; dsB = nullptr;
         
     } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
         ::Rf_error( "c++ exception bdcomputeMatrixVector_hdf5 (File IException)" );
     } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
         ::Rf_error( "c++ exception bdcomputeMatrixVector_hdf5 (DataSet IException)" );
     } catch(std::exception &ex) {
         Rcpp::Rcout<< ex.what();
         return void();
     }
     
     return void();
 }
