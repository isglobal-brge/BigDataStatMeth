#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixSdMean.hpp"



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
//' @param force, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
//' @return hdf5 data file containing two new datasets, one for sd (if sd is 
//' requested) and another for the mean (if mean is requested). Results are
//' stored inside a folder mean_sd inside hdf5 data file with names: 
//' sd.<dataset>, mean.<dataset> respectively
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
//' bdCreate_hdf5_matrix_file("test.hdf5", Y, "data", "Y", force = TRUE)
//' bdAdd_hdf5_matrix( X, "test.hdf5",  "data", "X", force = TRUE)
//' 
//' # Get mean and sd        
//' bdgetSDandMean_hdf5(filename = "test.hdf5", group = "data", dataset = "Y",
//'                     sd = TRUE, mean = TRUE,byrows = TRUE)
//'         
//' @export
// [[Rcpp::export]]
void bdgetSDandMean_hdf5( std::string filename, 
                          std::string group, std::string dataset,
                          Rcpp::Nullable<bool> sd = R_NilValue, 
                          Rcpp::Nullable<bool> mean  = R_NilValue,
                          Rcpp::Nullable<bool> byrows = R_NilValue,
                          Rcpp::Nullable<int> wsize  = R_NilValue, 
                          Rcpp::Nullable<bool> force  = false)
{
 
 
 bool bsd, bmean, bforce, bbyrows;
 hsize_t blocksize, nrows, ncols;
 
 std::string strgroupout;
 
 Eigen::MatrixXd datanormal;
 
 try {
     
     if( mean.isNull()) {   bmean = true;   } 
     else {  bmean = Rcpp::as<bool> (mean);   }
     
     if( sd.isNull()) {   bsd = true;   } 
     else {   bsd = Rcpp::as<bool> (sd);   }
     
     if( byrows.isNull()) {   bbyrows = false;   } 
     else {   bbyrows = Rcpp::as<bool> (byrows);   }
     
     if( force.isNull()) { bforce = false; } 
     else { bforce = Rcpp::as<bool> (force);  }
     
     BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
     dsA->openDataset();
     
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
     
     delete dsA;
     delete dssd;
     delete dsmean;
     
     } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
         Rcpp::Rcout<<"c++ exception bdgetSDandMean_hdf5 (File IException)";
         return void();
     } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdgetSDandMean_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdgetSDandMean_hdf5" << ex.what();
        return void();
    }
 
    return void();

}

