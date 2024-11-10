#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/tcrossprod.hpp"
#include "Utilities/Utilities.hpp"


//' tCrossprod with hdf5 matrix
//' 
//' This function performs the tcrossprod from a matrix inside and hdf5 data file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group, string, group name where dataset A is stored
//' @param A string name inside HDF5 file
//' @param groupB, string, group name where dataset b is stored
//' @param B string, dataset name for matrix B inside HDF5 file
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @param mixblock_size (optional) only for debug pourpose
//' @param outgroup (optional) group name to store results from tCrossprod inside hdf5 data file
//' @return no value
//' @examples
//'   
//'   library(BigDataStatMeth)
//'   library(rhdf5)
//'   
//'   N = 1000
//'   M = 1000
//'   
//'   set.seed(555)
//'   a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
//'   
//'   bdCreate_hdf5_matrix( filename = "test_temp.hdf5", 
//'                         object = a, group = "INPUT", 
//'                         dataset = "datasetA",
//'                         transp = FALSE,
//'                         overwriteFile = TRUE, 
//'                         overwriteDataset = FALSE, 
//'                         unlimited = FALSE)
//'                         
//'     file <- "test_temp.hdf5"
//'     dataset <- "results/res"
//'     
//'     bdtCrossprod_hdf5( filename = "test_temp.hdf5", group = "INPUT", 
//'                        A = "datasetA", outgroup = "results", 
//'                        outdataset = "res", force = TRUE ) # 
//'                        
//'     # Check results
//'     resr <- tcrossprod(a)
//'     res <-  h5read(file,dataset)
//'     all.equal( resr, res)
//'     
//'     bdtCrossprod_hdf5(filename = "test_temp.hdf5", group = "INPUT", 
//'                        A = "datasetA", outgroup = "results", 
//'                        outdataset = "res", block_size = 1024, 
//'                        force = TRUE ) # 
//'     
//'     # Check results
//'     resr <- tcrossprod(a)
//'     res <-  h5read(file,dataset)
//'     all.equal( resr, res)
//'   
//'     # Remove file (used as example)
//'     if (file.exists("test_temp.hdf5")) {
//'       file.remove("test_temp.hdf5")
//'     }
//'   
//' 
//' @export
// [[Rcpp::export]]
 void bdtCrossprod_hdf5( std::string filename, 
                        std::string group, 
                        std::string A, 
                        Rcpp::Nullable<std::string> B = R_NilValue, 
                        Rcpp::Nullable<std::string> groupB = R_NilValue, 
                        Rcpp::Nullable<int> block_size = R_NilValue,
                        Rcpp::Nullable<int> mixblock_size = R_NilValue,
                        Rcpp::Nullable<bool> paral = R_NilValue,
                        Rcpp::Nullable<int> threads = R_NilValue,
                        Rcpp::Nullable<std::string> outgroup = R_NilValue,
                        Rcpp::Nullable<std::string> outdataset = R_NilValue,
                        Rcpp::Nullable<bool> force = R_NilValue )                                
 {
     
     int iblock_size,
        iblockfactor = 2;
     bool bparal, bforce;
     
     std::string strsubgroupOut, 
     strdatasetOut, 
     strsubgroupIn,
     strsubgroupInB;
     std::string matB;
     
     try {
         
         H5::Exception::dontPrint();  
         
         strsubgroupIn = group;
         
         if( outgroup.isNull()) { strsubgroupOut = "OUTPUT";
         } else { strsubgroupOut = Rcpp::as<std::string> (outgroup); }
         
         if(B.isNotNull()){ matB =  Rcpp::as<std::string> (B) ; } 
         else { matB =  A; }
         
         if(groupB.isNotNull()){ strsubgroupInB =  Rcpp::as<std::string> (groupB) ; } 
         else { strsubgroupInB =  group; }
         
         if (paral.isNull()) { bparal = false; } 
         else { bparal = Rcpp::as<bool> (paral); }
         
         if (force.isNull()) { bforce = false; } 
         else { bforce = Rcpp::as<bool> (force); }
         
         if( outdataset.isNotNull()) { strdatasetOut =  Rcpp::as<std::string> (outdataset); } 
         else { strdatasetOut = "tCrossProd_" + A + "_x_" + matB; }
         
         
         BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupIn, A, false);
         dsA->openDataset();
         BigDataStatMeth::hdf5Dataset* dsB = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupInB, matB, false);
         dsB->openDataset();
         BigDataStatMeth::hdf5Dataset* dsC = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupOut, strdatasetOut, bforce);
         
         iblock_size = BigDataStatMeth::getMaxBlockSize( dsA->nrows(), dsA->ncols(), dsB->nrows(), dsB->ncols(), iblockfactor, block_size);
         
         if(bparal == true) { // parallel
             
             int memory_block; 
             if(mixblock_size.isNotNull()) {
                 memory_block = Rcpp::as<int> (mixblock_size);
             } else {
                 memory_block = iblock_size/2;
             }
             
             dsC = BigDataStatMeth::tcrossprod(dsA, dsB, dsC, iblock_size, memory_block, bparal, true, threads);
             
         } else if (bparal == false) { // Not parallel
             dsC = BigDataStatMeth::tcrossprod(dsA, dsB, dsC, iblock_size, 0, bparal, true, threads);
         }
         
         delete dsA;
         delete dsB;
         delete dsC;
         
     } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
         ::Rf_error( "c++ exception bdtCrossprod_hdf5 (File IException)" );
     } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
         ::Rf_error( "c++ exception bdtCrossprod_hdf5 (DataSet IException)" );
     } catch(std::exception &ex) {
     }
     
     // return List::create(Named("filename") = filename,
     //                     Named("dataset") = strsubgroupOut + "/" + strdatasetOut);
     return void();
     
     
 }