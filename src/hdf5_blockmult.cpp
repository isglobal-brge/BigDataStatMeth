#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/multiplication.hpp"

//' Hdf5 datasets multiplication
//'
//' Multiplies two existing datasets in hdf5 datafile and stores results i a new hdf5 dataset.
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
                    Rcpp::Nullable<bool> force = R_NilValue)
{
    
    int iblock_size, 
        iblockfactor = 2;
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
        
        // if (paral.isNull()) { bparal = false; } 
        // else { bparal = Rcpp::as<bool> (paral); }
        
        if (force.isNull()) { bforce = false; } 
        else { bforce = Rcpp::as<bool> (force); }
        
        dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupIn, A, false);
        dsA->openDataset();
        dsB = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupInB, B, false);
        dsB->openDataset();
        dsC = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupOut, strdatasetOut, bforce);
        
        //../ iblock_size = BigDataStatMeth::getMaxBlockSize( dsA->nrows(), dsA->ncols(), dsB->nrows(), dsB->ncols(), iblockfactor, block_size);

//.. 2024/03/28 ..//        if(bparal == true) { // parallel
            
            // int memory_block; 
            // if(mixblock_size.isNotNull()) {
            //     memory_block = Rcpp::as<int> (mixblock_size);
            // } else {
            //     memory_block = iblock_size/2;
            // }
            
            BigDataStatMeth::multiplication(dsA, dsB, dsC, paral, block_size, threads);

//.. 2024/03/28 ..//        } else if (bparal == false) { // Not parallel
//.. 2024/03/28 ..//            BigDataStatMeth::multiplication(dsA, dsB, dsC, paral, block_size, threads);
//.. 2024/03/28 ..//        }
        
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



/***R

library("BigDataStatMeth")
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

N = 200
M = 10

set.seed(555)
a <- matrix( rnorm(N*M,mean=0,sd=1), N, M) 
b <- matrix( rnorm(N*M,mean=0,sd=1), M, N) 

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = a, 
                     group = "pepet", 
                     dataset = "datasetpepet",
                     transp = FALSE,
                     overwriteFile = TRUE, 
                     overwriteDataset = FALSE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = b, 
                     group = "pepet", 
                     dataset = "tdatasetpepet",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = FALSE, 
                     unlimited = FALSE)
# blockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", force = TRUE )
blockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", block_size = 5, force = TRUE )


*/
