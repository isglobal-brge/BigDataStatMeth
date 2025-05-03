#include <BigDataStatMeth.hpp>
// #include "hdf5Algebra/matrixSdMean.hpp"



//' Get sd and Mean by Rows or Columns
//' 
//' This functions gets Standard Deviation (sd) or Mean by Rows or Columns and
//' store results in hdf5 dataset inside the file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string Matrix
//' @param dataset string Matrix
//' @param sd logical (default = TRUE) if TRUE, standard deviation is computed
//' @param mean logical (default = TRUE) if TRUE, mean is computed 
//' @param byrows logical (default = FALSE) if TRUE, sd and mean are computed
//' by columns, if byrows=TRUE then sd and mean are computed by Rows.
//' @param wsize integer (default = 1000), file block size to read to 
//' perform calculus exitexit
//' @param overwrite, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
//' @return hdf5 data file containing two new datasets, one for sd (if sd is 
//' requested) and another for the mean (if mean is requested). Results are
//' stored inside a folder mean_sd inside hdf5 data file with names: 
//' sd.\<dataset\>, mean.\<dataset\> respectively
//' @examples
//' 
//' library(BigDataStatMeth)
//'     
//' # Prepare data and functions
//' set.seed(123)
//' Y <- matrix(rnorm(100), 10, 10)
//' X <- matrix(rnorm(10), 10, 1)
//'     
//' # Create hdf5 data file with  data (Y)
//' bdCreate_hdf5_matrix("test.hdf5", Y, "data", "Y", 
//'                         overwriteFile = TRUE, 
//'                         overwriteDataset = FALSE, 
//'                         unlimited = FALSE)
//' bdCreate_hdf5_matrix( "test.hdf5", X, "data", "X", 
//'                        overwriteFile = FALSE, 
//'                        overwriteDataset = FALSE, 
//'                        unlimited = FALSE)
//' 
//' # Get mean and sd        
//' bdgetSDandMean_hdf5(filename = "test.hdf5", group = "data", dataset = "Y",
//'                     sd = TRUE, mean = TRUE,byrows = TRUE)
//'                     
//' # Remove file (used as example)
//' if (file.exists("test.hdf5")) {
//'   file.remove("test.hdf5")
//' }
//'         
//' @export
// [[Rcpp::export]]
void bdgetSDandMean_hdf5( std::string filename, 
                          std::string group, std::string dataset,
                          Rcpp::Nullable<bool> sd = R_NilValue, 
                          Rcpp::Nullable<bool> mean  = R_NilValue,
                          Rcpp::Nullable<bool> byrows = R_NilValue,
                          Rcpp::Nullable<int> wsize  = R_NilValue, 
                          Rcpp::Nullable<bool> overwrite  = false)
{
    
    BigDataStatMeth::hdf5Dataset* dsA = nullptr;
    BigDataStatMeth::hdf5Dataset* dsmean = nullptr;
    BigDataStatMeth::hdf5Dataset* dssd = nullptr;
 
    try {
        
        bool bforce, bbyrows;
        hsize_t nrows, ncols;
        
        std::string strgroupout;
        Eigen::MatrixXd datanormal;
         
         
         if( byrows.isNull()) {   bbyrows = false;   } 
         else {   bbyrows = Rcpp::as<bool> (byrows);   }
         
         if( overwrite.isNull()) { bforce = false; } 
         else { bforce = Rcpp::as<bool> (overwrite);  }
         
         
         dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
         dsA->openDataset();
         
         if( dsA->getDatasetptr() != nullptr) {
             
             nrows = dsA->nrows();
             ncols = dsA->ncols();
             
             // Define blocksize atending number of elements in rows and cols
             if( bbyrows == false) {
                 datanormal = Eigen::MatrixXd::Zero(2,nrows);
                 get_HDF5_mean_sd_by_column( dsA, datanormal, wsize);
             } else {
                 datanormal = Eigen::MatrixXd::Zero(2,ncols);
                 get_HDF5_mean_sd_by_row( dsA, datanormal, wsize);
             }
             
             strgroupout = "mean_sd";
             std::string strdatasetmean = "mean." + dataset;
             std::string strdatasetsd = "sd." + dataset;
             
             BigDataStatMeth::hdf5Dataset* dsmean = new BigDataStatMeth::hdf5Dataset(filename, strgroupout, strdatasetmean, bforce);
             dsmean->createDataset( datanormal.cols(), 1, "real");
             dsmean->writeDataset( Rcpp::wrap(datanormal.row(0)) );
             
             BigDataStatMeth::hdf5Dataset* dssd = new BigDataStatMeth::hdf5Dataset(filename, strgroupout, strdatasetsd, bforce);
             dssd->createDataset( datanormal.cols(), 1, "real");
             dssd->writeDataset( Rcpp::wrap(datanormal.row(1)) );    
         }
         
         
         delete dsA; dsA = nullptr;
         delete dssd; dssd = nullptr;
         delete dsmean; dsmean = nullptr;
     
     } catch( H5::FileIException& error ) { 
         checkClose_file(dsA, dssd, dsmean);
         Rcpp::Rcerr<<"c++ exception bdgetSDandMean_hdf5 (File IException)";
         return void();
     } catch( H5::DataSetIException& error ) { 
         checkClose_file(dsA, dssd, dsmean);
         Rcpp::Rcerr << "c++ exception bdgetSDandMean_hdf5 (DataSet IException)";
         return void();
     } catch(std::exception& ex) {
         checkClose_file(dsA, dssd, dsmean);
         Rcpp::Rcerr << "c++ exception bdgetSDandMean_hdf5" << ex.what();
         return void();
     } catch (...) {
         checkClose_file(dsA, dssd, dsmean);
         Rcpp::Rcerr<<"\nC++ exception bdgetSDandMean_hdf5 (unknown reason)";
         return void();
     }
 
    return void();

}

