#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixSubstract.hpp"


/**
 *  // // [[Rcpp::export(.blockSum_hdf5)]]
 */

//' Hdf5 datasets substract
//'
//' Substracts two existing datasets in hdf5 datafile and stores results i a new hdf5 dataset.
//' 
//' @export
// [[Rcpp::export]]
void bdblockSubstract_hdf5(std::string filename, 
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
    
    int iblock_size;
    bool bparal, bforce;
    
    std::string strsubgroupOut, strdatasetOut,
                strsubgroupIn, 
                strsubgroupInB, strGroupB;
    
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
        
        if (force.isNull()) { bforce = false; } 
        else { bforce = Rcpp::as<bool> (force); }
        
        if( outdataset.isNotNull()) { strdatasetOut =  Rcpp::as<std::string> (outdataset); } 
        else { strdatasetOut =  A + "_-_" + B; }
        
        
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
            Rcpp_block_matrix_substract_hdf5(dsA, dsB, dsC, iblock_size, bparal, threads);
        } else {
            
            if( irowsA==1 || icolsA==1 ) {
                Rcpp_block_matrix_vector_substract_hdf5(dsA, dsB, dsC, iblock_size, bparal, threads);
            } else {
                Rcpp_block_matrix_vector_substract_hdf5(dsB, dsA, dsC, iblock_size, bparal, threads);
            }
        }
        
        delete dsA;
        delete dsB;
        delete dsC;
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception bdblockSubstract_hdf5 (File IException)" );
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        error.printErrorStack();
        ::Rf_error( "c++ exception bdblockSubstract_hdf5 (DataSet IException)" );
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
    }
    
    return void();
}
