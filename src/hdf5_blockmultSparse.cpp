#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/multiplicationSparse.hpp"


//' Block matrix multiplication
//' 
//' This function performs a block matrix-matrix multiplication with numeric matrix
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string path inside hdf5 data file where matrix A is stored
//' @param A, string with dataset name where matrix is stored
//' @param B, string with dataset name where matrix is stored
//' @param groupB string path inside hdf5 data file where matrix B is stored
//' @param block_size integer, block size used to perform calculus
//' @param mixblock_size integer
//' @param outgroup string with the group name under the matrix will be stored
//' @param outdataset string with the dataset name to store results
//' @param force, boolean
//' 
//' @return list with filename and the group and dataset name under the results are stored
//' @examples
//' 
//' library(Matrix)
//' library(BigDataStatMeth)
//' 
//' k <- 1e3
//' set.seed(1)
//' x_sparse <- sparseMatrix(
//'     i = sample(x = k, size = k),
//'     j = sample(x = k, size = k),
//'     x = rnorm(n = k)
//' )
//' set.seed(2)
//' y_sparse <- sparseMatrix(
//'     i = sample(x = k, size = k),
//'     j = sample(x = k, size = k),
//'     x = rnorm(n = k)
//' )
//' 
//' if( isTRUE(file.exists('BasicMatVect.hdf5'))) {
//'      file.remove('BasicMatVect.hdf5')
//' }
//' bdCreate_hdf5_matrix_file("BasicMatVect.hdf5", as.matrix(x_sparse), "SPARSE", "x_sparse")
//' bdAdd_hdf5_matrix(as.matrix(y_sparse), "BasicMatVect.hdf5", "SPARSE", "y_sparse")
//' 
//' d <- bdblockmult_sparse_hdf5("BasicMatVect.hdf5", "SPARSE", "x_sparse", "y_sparse")
//' 
//' # Remove file (used as example)
//' if (file.exists("BasicMatVect.hdf5")) {
//'   # Delete file if it exist
//'   file.remove("BasicMatVect.hdf5")
//' }
//' 
//' @export
// [[Rcpp::export]]
void bdblockmult_sparse_hdf5( std::string filename, std::string group, 
                          std::string A, std::string B,
                          Rcpp::Nullable<std::string> groupB = R_NilValue, 
                          Rcpp::Nullable<int> block_size = R_NilValue,
                          Rcpp::Nullable<int> mixblock_size = R_NilValue,
                          Rcpp::Nullable<bool> paral = R_NilValue,
                          Rcpp::Nullable<int> threads = R_NilValue,
                          Rcpp::Nullable<std::string> outgroup = R_NilValue,
                          Rcpp::Nullable<std::string> outdataset = R_NilValue,
                          Rcpp::Nullable<bool> force = R_NilValue )
{
     
    
    std::string strdataset;
    std::string spMatrix("dgCMatrix");
    
    int iblock_size;
    bool bparal, bforce;
    
    std::string strsubgroupOut,
                strdatasetOut, 
                strsubgroupIn,
                strsubgroupInB;
    
    bool bexistgroup;
    bool bsparseA, bsparseB;
    int iblockfactor = 2;

    try
    {
     
        H5::Exception::dontPrint();  
        
        strsubgroupIn = group;
        
        if( outgroup.isNull()) { strsubgroupOut = "OUTPUT";
        } else { strsubgroupOut = Rcpp::as<std::string> (outgroup); }
        
        if(groupB.isNotNull()){ strsubgroupInB =  Rcpp::as<std::string> (groupB) ; } 
        else { strsubgroupInB =  group; }
        
        if( outdataset.isNotNull()) { strdatasetOut =  Rcpp::as<std::string> (outdataset); } 
        else { strdatasetOut =  A + "_x_" + B; }
        
        if (paral.isNull()) { bparal = false; } 
        else { bparal = Rcpp::as<bool> (paral); }
        
        if (force.isNull()) { bforce = false; } 
        else { bforce = Rcpp::as<bool> (force); }
        
        
        BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupIn, A, false);
        dsA->openDataset();
        BigDataStatMeth::hdf5Dataset* dsB = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupInB, B, false);
        dsB->openDataset();
        BigDataStatMeth::hdf5Dataset* dsC = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupOut, strdatasetOut, bforce);
        
        iblock_size = BigDataStatMeth::getMaxBlockSize( dsA->nrows(), dsA->ncols(), dsB->nrows(), dsB->ncols(), iblockfactor, block_size);
     
     
     // if (block_size.isNotNull()) {
     //     iblock_size = Rcpp::as<int> (block_size);
     // } else {
     //     iblock_size = std::min(  std::min(dsA->nrows(),dsA->ncols()),  std::min(dsB->nrows(), dsB->ncols()));
     //     if (iblock_size>1024)
     //         iblock_size = 1024;
     // }
     
         if(bparal == true) { // parallel
             
             int memory_block; 
             if(mixblock_size.isNotNull()) {
                 memory_block = Rcpp::as<int> (mixblock_size);
             } else {
                 memory_block = iblock_size/2;
             }
             
             dsC = BigDataStatMeth::multiplicationSparse(dsA, dsB, dsC, iblock_size, memory_block, bparal, true, threads);
             
         } else if (bparal == false) { // Not parallel
             dsC = BigDataStatMeth::multiplicationSparse(dsA, dsB, dsC, iblock_size, 0, bparal, true, threads);
         }
     
         delete dsA; 
         delete dsB; 
         delete dsC; 
     
     
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception bdblockmult_sparse_hdf5 (File IException)" );
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception bdblockmult_sparse_hdf5 (DataSet IException)" );
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
    }
    
    return void();
     
     
}




