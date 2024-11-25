#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixSum.hpp"
#include "Utilities/Utilities.hpp"


//' Hdf5 datasets sum
//'
//' Sum two existing datasets in hdf5 datafile and stores results i a new hdf5 dataset.
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string with the group name where matrix is stored inside HDF5 file
//' @param A string, datasetname with matrix to be multiplied
//' @param B string, datasetname with matrix to be multiplied
//' @param groupB, string, (optional) group name where dataset B is stored, if empty group folder is used
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, boolean (optional, default = FALSE) set paral = true to force parallel execution
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @param outgroup (optional) string with group name where we want to store the result matrix, by default out group = "OUTGROUP"
//' @param outdataset (optional) string with dataset name where we want to store the results
//' @param overwrite (optional) either a logical value indicating whether the results must be overwritten or not.
//' 
//' @return a dataset inside the hdf5 data file with A+B 
//' 
//' @examples
//' library("BigDataStatMeth")
//' 
//' N = 1500;  M = 1500
//' 
//' set.seed(555)
//' a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
//' b <- matrix( rnorm( N*M, mean=0, sd=1), M, N) 
//' 
//' bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
//'                      object = a, group = "groupA", 
//'                      dataset = "datasetA",
//'                      transp = FALSE,
//'                      overwriteFile = TRUE, 
//'                      overwriteDataset = FALSE, 
//'                      unlimited = FALSE)
//'                      
//' bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
//'                      object = t(b), 
//'                      group = "groupA", 
//'                      dataset = "datasetB",
//'                      transp = FALSE,
//'                      overwriteFile = FALSE, 
//'                      overwriteDataset = TRUE, 
//'                      unlimited = FALSE)
//'                      
//' # Multiply two matrix
//' bdblockSum_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", overwrite = TRUE ) 
//' bdblockSum_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 1024, overwrite = TRUE )  
//' 
//' @export
// [[Rcpp::export]]
void bdblockSum_hdf5(std::string filename, 
                   std::string group, 
                   std::string A, 
                   std::string B,
                   Rcpp::Nullable<std::string> groupB = R_NilValue, 
                   Rcpp::Nullable<int> block_size = R_NilValue, 
                   Rcpp::Nullable<bool> paral = R_NilValue,
                   Rcpp::Nullable<int> threads = R_NilValue,
                   Rcpp::Nullable<std::string> outgroup = R_NilValue,
                   Rcpp::Nullable<std::string> outdataset = R_NilValue,
                   Rcpp::Nullable<bool> overwrite = R_NilValue)
{
    
    int iblock_size;
    bool bparal, 
         bforce;
    
    std::string strsubgroupOut, 
                strdatasetOut, 
                strsubgroupIn,
                strsubgroupInB,
                strGroupB;
    
    try{
        
        H5::Exception::dontPrint();  
        
        if( outgroup.isNull()) { strsubgroupOut = "OUTPUT"; } 
        else { strsubgroupOut = Rcpp::as<std::string> (outgroup); }
        
        strsubgroupIn = group + "/";
        
        if(groupB.isNotNull()){
            strsubgroupInB =  Rcpp::as<std::string> (groupB) + "/";
            strGroupB = Rcpp::as<std::string> (groupB);
        } else {
            strsubgroupInB =  group + "/";
            strGroupB = group;
        }
        
        if (paral.isNull()) { bparal = false; } 
        else { bparal = Rcpp::as<bool> (paral); }
        
        if (overwrite.isNull()) { bforce = false; } 
        else { bforce = Rcpp::as<bool> (overwrite); }
        
        if( outdataset.isNotNull()) { strdatasetOut =  Rcpp::as<std::string> (outdataset); } 
        else { strdatasetOut =  A + "_+_" + B; }
        
        
        BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupIn, A, false);
        dsA->openDataset();
        BigDataStatMeth::hdf5Dataset* dsB = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupInB, B, false);
        dsB->openDataset();
        BigDataStatMeth::hdf5Dataset* dsC = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupOut, strdatasetOut, bforce);
        
        int irowsA = dsA->nrows(),
            icolsA = dsA->ncols(),
            irowsB = dsB->nrows(),
            icolsB = dsB->ncols();
        
        if (block_size.isNotNull()) {
            iblock_size = Rcpp::as<int> (block_size);
        } else {
            
            if( irowsA == 1 || icolsA == 1 || irowsB == 1 || icolsB == 1){
                iblock_size = BigDataStatMeth::getVectorBlockSize( irowsA*icolsA);
            } else{
                std::vector<hsize_t> blockSize = BigDataStatMeth::getMatrixBlockSize( irowsA, icolsA);
                if(irowsA < icolsA) {
                    iblock_size = blockSize.at(0);    
                } else {
                    iblock_size = blockSize.at(1);
                }
            }
        }
        
        if( irowsA != 1 && icolsA!= 1 && irowsB != 1 && icolsB!= 1) {
            Rcpp_block_matrix_sum_hdf5(dsA, dsB, dsC, iblock_size, bparal, threads);
        } else {

            if( irowsA==1 || icolsA==1 ) {
                Rcpp_block_matrix_vector_sum_hdf5(dsA, dsB, dsC, iblock_size, bparal, threads);
            } else {
                Rcpp_block_matrix_vector_sum_hdf5(dsB, dsA, dsC, iblock_size, bparal, threads);
            }
        }
        
        delete dsA;
        delete dsB;
        delete dsC;
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception bdblockSum_hdf5 (File IException)" );
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        error.printErrorStack();
        ::Rf_error( "c++ exception bdblockSum_hdf5 (DataSet IException)" );
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
    }
    
    // //..// return(C);
    // return List::create(Named("filename") = filename,
    //                     Named("dataset") = strsubgroupOut + "/" + strdatasetOut,
    //                     Named("result") = wrap(0));
    
    return void();
}
