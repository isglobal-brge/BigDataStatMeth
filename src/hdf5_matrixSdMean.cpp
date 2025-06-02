#include <BigDataStatMeth.hpp>
// #include "hdf5Algebra/matrixSdMean.hpp"

/**
 * @file hdf5_matrixSdMean.cpp
 * @brief Implementation of matrix statistics (SD and mean) for HDF5-stored matrices
 * @details This file provides functionality for computing standard deviation
 * and mean statistics for matrices stored in HDF5 format. The implementation
 * supports:
 * - Row-wise and column-wise computations
 * - Block-based processing for large matrices
 * - Parallel computation capabilities
 * - Memory-efficient operations
 * 
 * Key features:
 * - Support for large matrices
 * - Configurable block size
 * - Flexible computation direction
 * - Memory-efficient implementation
 * - Comprehensive error handling
 */

/**
 * @brief Computes standard deviation and mean of HDF5 matrix
 * 
 * @details Implements efficient computation of standard deviation and mean
 * statistics for matrices stored in HDF5 format. The function supports both
 * row-wise and column-wise computations with block-based processing.
 * 
 * Implementation features:
 * - Block-based computation
 * - Configurable processing direction
 * - Memory-efficient operations
 * - Safe file operations
 * - Comprehensive error handling
 * 
 * @param filename Path to HDF5 file
 * @param group Group containing dataset
 * @param dataset Dataset name
 * @param sd Whether to compute standard deviation
 * @param mean Whether to compute mean
 * @param byrows Whether to compute by rows
 * @param wsize Block size for processing
 * @param overwrite Whether to overwrite existing results
 * 
 * @throws H5::FileIException for HDF5 file operation errors
 * @throws H5::DataSetIException for HDF5 dataset operation errors
 * @throws std::exception for other errors
 */

//' Compute Matrix Standard Deviation and Mean in HDF5
//'
//' @description
//' Computes standard deviation and/or mean statistics for a matrix stored in
//' HDF5 format, with support for row-wise or column-wise computations.
//'
//' @details
//' This function provides efficient statistical computation capabilities with:
//' 
//' * Computation options:
//'   - Standard deviation computation
//'   - Mean computation
//'   - Row-wise or column-wise processing
//' 
//' * Processing features:
//'   - Block-based computation
//'   - Memory-efficient processing
//'   - Configurable block size
//' 
//' * Implementation features:
//'   - Safe HDF5 file operations
//'   - Memory-efficient implementation
//'   - Comprehensive error handling
//'
//' Results are stored in a new group 'mean_sd' within the HDF5 file.
//'
//' @param filename Character string. Path to the HDF5 file.
//' @param group Character string. Path to the group containing the dataset.
//' @param dataset Character string. Name of the dataset to analyze.
//' @param sd Logical (optional). Whether to compute standard deviation.
//'   Default is TRUE.
//' @param mean Logical (optional). Whether to compute mean. Default is TRUE.
//' @param byrows Logical (optional). Whether to compute by rows (TRUE) or
//'   columns (FALSE). Default is FALSE.
//' @param wsize Integer (optional). Block size for processing. Default is 1000.
//' @param overwrite Logical (optional). Whether to overwrite existing results.
//'   Default is FALSE.
//'
//' @return No return value, called for side effects (statistics computation).
//'   Results are stored in the HDF5 file under the 'mean_sd' group with names:
//'   * 'sd.`dataset`' for standard deviation
//'   * 'mean.`dataset`' for mean
//'
//' @examples
//' \dontrun{
//' library(BigDataStatMeth)
//' 
//' # Create test matrices
//' set.seed(123)
//' Y <- matrix(rnorm(100), 10, 10)
//' X <- matrix(rnorm(10), 10, 1)
//' 
//' # Save to HDF5
//' bdCreate_hdf5_matrix("test.hdf5", Y, "data", "matrix1",
//'                      overwriteFile = TRUE)
//' bdCreate_hdf5_matrix("test.hdf5", X, "data", "vector1",
//'                      overwriteFile = FALSE)
//' 
//' # Compute statistics
//' bdgetSDandMean_hdf5(
//'   filename = "test.hdf5",
//'   group = "data",
//'   dataset = "matrix1",
//'   sd = TRUE,
//'   mean = TRUE,
//'   byrows = TRUE,
//'   wsize = 500
//' )
//' 
//' # Cleanup
//' if (file.exists("test.hdf5")) {
//'   file.remove("test.hdf5")
//' }
//' }
//'
//' @references
//' * The HDF Group. (2000-2010). HDF5 User's Guide.
//' * Welford, B. P. (1962). Note on a method for calculating corrected
//'   sums of squares and products. Technometrics, 4(3), 419-420.
//'
//' @seealso
//' * \code{\link{bdCreate_hdf5_matrix}} for creating HDF5 matrices
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
             if( dsmean->getDatasetptr() != nullptr ) {
                 dsmean->writeDataset( Rcpp::wrap(datanormal.row(0)) );
             } else {
                 checkClose_file(dsA, dsmean);
                 Rcpp::Rcerr << "c++ exception bdgetSDandMean_hdf5: " << "Error creating dataset";
                 return void();
             }
             
             BigDataStatMeth::hdf5Dataset* dssd = new BigDataStatMeth::hdf5Dataset(filename, strgroupout, strdatasetsd, bforce);
             dssd->createDataset( datanormal.cols(), 1, "real");
             if( dssd->getDatasetptr() != nullptr ) {
                 dssd->writeDataset( Rcpp::wrap(datanormal.row(1)) );
             } else {
                 checkClose_file(dsA, dssd, dsmean);
                 Rcpp::Rcerr << "c++ exception bdgetSDandMean_hdf5: " << "Error creating dataset";
                 return void();
             }
             
         } else {
             checkClose_file(dsA);
             Rcpp::Rcerr << "c++ exception bdgetSDandMean_hdf5: " << "Error opening dataset";
             return void();
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

