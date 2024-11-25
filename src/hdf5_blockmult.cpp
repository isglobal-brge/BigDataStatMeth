#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/multiplication.hpp"

//' Hdf5 datasets multiplication
//'
//' Multiplies two existing datasets in hdf5 datafile and stores results i a new hdf5 dataset.
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string with the group name where matrix is stored inside HDF5 file
//' @param A string, datasetname with matrix to be multiplied
//' @param B string, datasetname with matrix to be multiplied
//' @param groupB, string, (optional) group name where dataset B is stored, if empty group folder is used
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, boolean (optional, default = FALSE) set paral = true to force parallel execution
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @param outgroup (optional) string with group name where we want to store the result matrix
//' @param outdataset (optional) string with dataset name where we want to store the results
//' @param overwrite (optional) either a logical value indicating whether the results must be overwritten or not.
//' @return a dataset inside the hdf5 data file with A*B 
//' @examples
//' library("BigDataStatMeth")
//' 
//' N = 1500
//' M = 1500
//' set.seed(555)
//' a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
//' b <- matrix( rnorm( N*M, mean=0, sd=1), M, N) 
//' bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
//'                      object = a, group = "groupA", 
//'                      dataset = "datasetA",
//'                      transp = FALSE,
//'                      overwriteFile = TRUE, 
//'                      overwriteDataset = FALSE, 
//'                      unlimited = FALSE)
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
//' bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", overwrite = TRUE ) 
//' bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 1024, overwrite = TRUE ) 
//' 
//' @export
// [[Rcpp::export]]
void bdblockmult_hdf5(std::string filename, 
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
    
    // int iblock_size;
        // iblockfactor = 2;
    bool bforce;
//.. 2024/03/28 ..//    bool bparal, bforce;
    
    std::string strsubgroupOut, 
                strdatasetOut, 
                strsubgroupIn,
                strsubgroupInB;
    
    BigDataStatMeth::hdf5Dataset* dsA;
    BigDataStatMeth::hdf5Dataset* dsB;
    BigDataStatMeth::hdf5Dataset* dsC;
    
    try {
        
        H5::Exception::dontPrint();  
        
        strsubgroupIn = group;
        
        if( outgroup.isNull()) { strsubgroupOut = "OUTPUT";
        } else { strsubgroupOut = Rcpp::as<std::string> (outgroup); }
        
        if(groupB.isNotNull()){ strsubgroupInB =  Rcpp::as<std::string> (groupB) ; } 
        else { strsubgroupInB =  group; }
        
        if( outdataset.isNotNull()) { strdatasetOut =  Rcpp::as<std::string> (outdataset); } 
        else { strdatasetOut =  A + "_x_" + B; }
        
        if (overwrite.isNull()) { bforce = false; } 
        else { bforce = Rcpp::as<bool> (overwrite); }
        
        dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupIn, A, false);
        dsA->openDataset();
        dsB = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupInB, B, false);
        dsB->openDataset();
        dsC = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupOut, strdatasetOut, bforce);
        
        BigDataStatMeth::multiplication(dsA, dsB, dsC, paral, block_size, threads);

        delete dsA;
        delete dsB;
        delete dsC;
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        delete dsA; delete dsB; delete dsC;
        ::Rf_error( "c++ exception blockmult_hdf5 (File IException)" );
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        delete dsA; delete dsB; delete dsC;
        ::Rf_error( "c++ exception blockmult_hdf5 (DataSet IException)" );
    } catch(std::exception &ex) {
        delete dsA; delete dsB; delete dsC;
    }
    
    // return List::create(Named("filename") = filename,
    //                     Named("dataset") = strsubgroupOut + "/" + strdatasetOut);
    return void();
}

