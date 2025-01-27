#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/crossprod.hpp"
#include "Utilities/Utilities.hpp"


//' Crossprod with hdf5 matrix
//' 
//' This function performs the cross product of one or two matrices inside 
//' and hdf5 data file
//' 
//' @inheritParams bdblockmult_hdf5
//' @param outdataset string (optional), An optional parameter specifying the 
//' dataset name for the output matrix. If NULL, the default name will be 
//' constructed as "CrossProd_" concatenated with the name of dataset A 
//' "_x_" and the name of dataset B.
//' @param mixblock_size only for debug pourpose
//' @details
//' For a single matrix \eqn{A}, the cross product is defined as \eqn{A^t A}, 
//' where \eqn{A^t} is the transpose of \eqn{A}. For two matrices \eqn{A} and 
//' \eqn{B}, the cross product is \eqn{A^5 B}. This operation is often used in 
//' linear algebra for projections and other computations.
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
//'     bdCrossprod_hdf5( filename = "test_temp.hdf5", group = "INPUT", 
//'                        A = "datasetA", outgroup = "results", 
//'                        outdataset = "res", overwrite = TRUE ) # 
//'                        
//'     # Check results
//'     resr <- tcrossprod(a)
//'     res <-  h5read(file,dataset)
//'     all.equal( resr, res)
//'     
//'     bdCrossprod_hdf5(filename = "test_temp.hdf5", group = "INPUT", 
//'                        A = "datasetA", outgroup = "results", 
//'                        outdataset = "res", block_size = 1024, 
//'                        overwrite = TRUE ) # 
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
void bdCrossprod_hdf5( std::string filename, 
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
                       Rcpp::Nullable<bool> overwrite = R_NilValue )                                
{
    
    int iblock_size,
        iblockfactor = 2;
    bool bparal, bforce;
    
    BigDataStatMeth::hdf5Dataset* dsA;
    BigDataStatMeth::hdf5Dataset* dsB;
    BigDataStatMeth::hdf5Dataset* dsC;
    
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
        
        if (overwrite.isNull()) { bforce = false; } 
        else { bforce = Rcpp::as<bool> (overwrite); }
        
        if( outdataset.isNotNull()) { strdatasetOut =  Rcpp::as<std::string> (outdataset); } 
        else { strdatasetOut = "CrossProd_" + A + "_x_" + matB; }
        
        
        dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupIn, A, false);
        dsA->openDataset();
        dsB = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupInB, matB, false);
        dsB->openDataset();
        
        if( dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            return void();
        
        dsC = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupOut, strdatasetOut, bforce);
        
        iblock_size = BigDataStatMeth::getMaxBlockSize( dsA->nrows(), dsA->ncols(), dsB->nrows(), dsB->ncols(), iblockfactor, block_size);

        if(bparal == true) { // parallel
            
            int memory_block; 
            if(mixblock_size.isNotNull()) {
                memory_block = Rcpp::as<int> (mixblock_size);
            } else {
                memory_block = iblock_size/2;
            }
            
            dsC = BigDataStatMeth::crossprod(dsA, dsB, dsC, iblock_size, memory_block, bparal, true, threads);
            
        } else if (bparal == false) { // Not parallel
            dsC = BigDataStatMeth::crossprod(dsA, dsB, dsC, iblock_size, 0, bparal, true, threads);
        }
        
        delete dsA;
        delete dsB;
        delete dsC;
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        if( dsA->isOpen()) dsA->close_file();
        Rcpp::Rcerr<<"\nc++ c++ exception bdCrossprod_hdf5 (File IException)\n";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        if( dsA->isOpen()) dsA->close_file();
        Rcpp::Rcerr<<"\nc++ exception bdCrossprod_hdf5 (DataSet IException)\n";
        return void();
    } catch(std::exception &ex) {
        if( dsA->isOpen()) dsA->close_file();
        Rcpp::Rcerr<<"\nc++ exception bdCrossprod_hdf5\n";
        return void();
    }
    
    // return List::create(Named("filename") = filename,
    //                     Named("dataset") = strsubgroupOut + "/" + strdatasetOut);
    return void();
    
}


/*** 
//' @param filename string file name where dataset to normalize is stored
//' @param group, string, group name where dataset A is stored
//' @param A string name inside HDF5 file
//' @param groupB, string, group name where dataset b is stored
//' @param B string, dataset name for matrix B inside HDF5 file
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @param outgroup (optional) group name to store Crossprod results inside hdf5 data file
//' @param outdataset (optional) dataset name to store Crossprod results inside hdf5 data file
//' @param overwrite, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
**/