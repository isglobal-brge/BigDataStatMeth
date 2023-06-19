#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixSdMean.hpp"
#include "hdf5Algebra/matrixNormalization.hpp"


//' Normalize dataset in hdf5 file
//' 
//' This function normalize data scaling, centering or scaling and centering 
//' in a dataset stored in hdf5 file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string Matrix
//' @param dataset string Matrix
//' @param bcenter logical (default = TRUE) if TRUE, centering is done by 
//' subtracting the column means
//' @param bscale logical (default = TRUE) if TRUE, centering is done by 
//' subtracting the column means
//' @param byrows logical (default = FALSE) if TRUE, centering is done by 
//' subtracting the rows means, util when working with hdf5 datasets stored 
//' in Row Major format.
//' @param wsize integer (default = 1000), file block size to read to 
//' perform normalization
//' @param force, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
//' @return file with scaled, centered or scaled and centered dataset
//' @examples
//'   a = "See vignette"
//' @export
// [[Rcpp::export]]
void bdNormalize_hdf5( std::string filename, std::string group, std::string dataset,
                    Rcpp::Nullable<bool> bcenter = R_NilValue, Rcpp::Nullable<bool> bscale  = R_NilValue,
                    Rcpp::Nullable<bool> byrows = R_NilValue,
                    Rcpp::Nullable<int> wsize  = R_NilValue, Rcpp::Nullable<bool> force  = false)
{
 
     bool bc, bs, bforce, bbyrows, bgetTransposed;
     hsize_t blocksize, nrows, ncols, nRowsCols;
     std::string strgroupout;
     std::vector<hsize_t> stride = {1, 1},
                          block = {1, 1};
     
     Eigen::MatrixXd datanormal;
     
     try{
         
         if( bcenter.isNull()) {  bc = true;   }
         else {   bc = Rcpp::as<bool> (bcenter);  }
         
         if( bscale.isNull()) {   bs = true;   }
         else {   bs = Rcpp::as<bool> (bscale);   }
         
         if( byrows.isNull()) {  bbyrows = false;   }
         else {   bbyrows = Rcpp::as<bool> (byrows);   }
         
         if( force.isNull()) {   bforce = false;   }
         else {   bforce = Rcpp::as<bool> (force);   }
         
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
         
         strgroupout = "NORMALIZED/" + group;
         std::string strdatasetmean = "mean." + dataset;
         std::string strdatasetsd = "sd." + dataset;
         
         BigDataStatMeth::hdf5Dataset* dsmean = new BigDataStatMeth::hdf5Dataset(filename, strgroupout, strdatasetmean, bforce);
         dsmean->createDataset( datanormal.cols(), 1, "real");
         dsmean->writeDataset( Rcpp::wrap(datanormal.row(0)) );
         
         BigDataStatMeth::hdf5Dataset* dssd = new BigDataStatMeth::hdf5Dataset(filename, strgroupout, strdatasetsd, bforce);
         dssd->createDataset( datanormal.cols(), 1, "real");
         dssd->writeDataset( Rcpp::wrap(datanormal.row(1)) );
         
         delete dssd;
         delete dsmean;
         
         if(bbyrows == false ) {
             nRowsCols = nrows;
             bgetTransposed = true;
         } else {
             nRowsCols = ncols;
             bgetTransposed = false;
         }
         
         BigDataStatMeth::hdf5Dataset* dsNormal;
         blocksize = BigDataStatMeth::get_block_size(wsize, nrows, ncols);
         
         for(hsize_t i=0; i*blocksize <= nRowsCols ; i++)
         {
             std::vector<hsize_t> offset, 
                                  count;
             hsize_t sizetoread = 0;
             
             if( (i+1) * blocksize < nRowsCols ) {
                 sizetoread = blocksize;
             } else {
                 sizetoread = nRowsCols - ( i * blocksize );
             }
             
             // Prepare file and dataset
             if(i==0) {
                 dsNormal = new BigDataStatMeth::hdf5Dataset(filename, strgroupout, dataset, bforce);
                 dsNormal->createDataset( dsA, "real");
             }
             
             if(bbyrows == false) {
                 offset = { i*blocksize, 0 };
                 count = { sizetoread, ncols };    
             } else {
                 offset = {  0, i*blocksize };
                 count = { nrows, sizetoread };
             }
             
             // Normalize and write data
             std::vector<double> vdA( count[0] * count[1] ); 
             dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdA.data(), count[0], count[1] );
             
             
             if(bbyrows == false) {
                 X = BigDataStatMeth::RcppNormalize_Data_R_hdf5(X, bc, bs, bgetTransposed, datanormal.block(0,offset[0], 2, count[0]));    
             } else {
                 X = BigDataStatMeth::RcppNormalize_Data_R_hdf5(X, bc, bs, bgetTransposed, datanormal.block(0,offset[1], 2, count[1]));
             }

             dsNormal->writeRowMajorDatasetBlock( X, offset, count, stride, block);
             
         }
         
         delete dsNormal;
         delete dsA;

     } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
         Rcpp::Rcout<<"c++ exception bdNormalize_hdf5 (File IException)";
         return void();
     } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
         Rcpp::Rcout << "c++ exception bdNormalize_hdf5 (DataSet IException)";
         return void();
     } catch(std::exception& ex) {
         Rcpp::Rcout << "c++ exception bdNormalize_hdf5" << ex.what();
         return void();
     }

     return void();
}

