#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/substract.hpp"


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
    
    int iblock_size,
        res;
    bool bparal, 
         bforce, 
         bexistgroup;
    
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
        
        if (force.isNull()) { bforce = false; } 
        else { bforce = Rcpp::as<bool> (force); }
        
        if( outdataset.isNotNull()) { strdatasetOut =  Rcpp::as<std::string> (outdataset); } 
        else { strdatasetOut =  A + "_-_" + B; }
        
        
        BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupIn, A, false);
        dsA->openDataset();
        BigDataStatMeth::hdf5Dataset* dsB = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupInB, B, false);
        dsB->openDataset();
        BigDataStatMeth::hdf5Dataset* dsC = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupOut, strdatasetOut, bforce);
        
        if (block_size.isNotNull()) {
            iblock_size = Rcpp::as<int> (block_size);
        } else {
            iblock_size = std::min(  std::min(dsA->nrows(),dsA->ncols()),  std::min(dsB->nrows(), dsB->ncols()));
            if (iblock_size>8192)
                iblock_size = 8192;
        }

        if( dsA->nrows() != 1 && dsA->ncols()!= 1 && dsB->nrows() != 1 && dsB->ncols()!= 1) {
            
            hdf5_block_matrix_substract_hdf5_indatasets_transposed(dsA, dsB, dsC, iblock_size, bparal, true, threads);
        } else {

            if(dsA->nrows()==1 || dsA->ncols()==1) {
                hdf5_block_matrix_vector_substract_hdf5_transposed(dsA, dsB, dsC, iblock_size, bparal, true, threads);
            } else {
                hdf5_block_matrix_vector_substract_hdf5_transposed(dsB, dsA, dsC, iblock_size, bparal, true, threads);
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
    
    // //..// return(C);
    // return List::create(Named("filename") = filename,
    //                     Named("dataset") = strsubgroupOut + "/" + strdatasetOut,
    //                     Named("result") = wrap(0));
    
    return void();
}
