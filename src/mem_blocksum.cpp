#include <BigDataStatMeth.hpp>
#include "memAlgebra/memSum.hpp"

/**
 *  // // [[Rcpp::export(.blockSum_hdf5)]]
 */


//' @title Hdf5 datasets sum
//' @description Sum two existing datasets in hdf5 datafile and stores results i a new hdf5 dataset
//' @param A Matrix A
//' @param B Matrix B
//' @param block_size PARAM_DESCRIPTION, Default: NULL
//' @param paral if paral = TRUE performs parallel computation else performs seria computation, Default: FALSE
//' @param byBlocks If data matrix has more than 2.25e+08 (15000 x 15000) elements, by default the addition is done by blocks, but it can be forced not to be partitioned with parameter byblocks = FALSE, Default: TRUE
//' @param threads only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  (maximum number of threads available / 2), Default: NULL
//' @return new matrix with A + B
//' @examples 
//' \dontrun{
//' if(interactive()){
//'     N <- 2500
//'     M <- 400
//'     nc <-  4
//'     
//'     set.seed(555)
//'     mat <- matrix( rnorm( N*M, mean=0, sd=10), N, M) 
//'     
//'     sum_mem = bdblockSum(mat, mat, paral = TRUE, threads = nc)
//'  }
//' }
//' @rdname bdblockSum
//' @export 
// [[Rcpp::export]]
Rcpp::RObject bdblockSum(Rcpp::RObject A, Rcpp::RObject B,
                         Rcpp::Nullable<int> block_size = R_NilValue, 
                         Rcpp::Nullable<bool> paral = R_NilValue,
                         Rcpp::Nullable<bool> byBlocks = true,
                         Rcpp::Nullable<int> threads = R_NilValue)
{
    
    hsize_t iblock_size;
    bool bparal, bbyBlocks;
    
    Rcpp::NumericMatrix C;

    try{
        
        // if (paral.isNull()) { bparal = false; }
        // else { bparal = Rcpp::as<bool> (paral); }
        
        if (byBlocks.isNull()) { bbyBlocks = false; }
        else { bbyBlocks = Rcpp::as<bool> (byBlocks); }

        if( bparal==false || Rcpp::as<Rcpp::NumericVector>(A).size() < MAXELEMSINBLOCK || bbyBlocks == false) {
            if( Rcpp::is<Rcpp::NumericMatrix>(A) && Rcpp::is<Rcpp::NumericMatrix>(B) ) {
                return( BigDataStatMeth::Rcpp_matrix_sum(A, B) );
                
            } else if( Rcpp::is<Rcpp::NumericVector>(A) && Rcpp::is<Rcpp::NumericMatrix>(B)) {
                return( BigDataStatMeth::Rcpp_matrix_vect_sum( B, A) );
                
            } else if( Rcpp::is<Rcpp::NumericVector>(B) && Rcpp::is<Rcpp::NumericMatrix>(A)) {
                return( BigDataStatMeth::Rcpp_matrix_vect_sum( A, B) );
                
            } else if(Rcpp::is<Rcpp::NumericVector>(A) && Rcpp::is<Rcpp::NumericVector>(B)) {
                return( BigDataStatMeth::Rcpp_vector_sum(A, B));
                
            } else {
                Rcpp::Rcout<<"\nData type not allowed";
            }    
        } else {
            
            if( Rcpp::is<Rcpp::NumericMatrix>(A) && Rcpp::is<Rcpp::NumericMatrix>(B) ) {
                return( BigDataStatMeth::Rcpp_matrix_blockSum(A, B, threads) );
                
            } else if( Rcpp::is<Rcpp::NumericVector>(A) && Rcpp::is<Rcpp::NumericMatrix>(B)) {
                 return(BigDataStatMeth::Rcpp_matrix_vector_blockSum(B, A, paral, threads));
                
            } else if( Rcpp::is<Rcpp::NumericVector>(B) && Rcpp::is<Rcpp::NumericMatrix>(A)) {
                return(BigDataStatMeth::Rcpp_matrix_vector_blockSum(A, B, paral, threads));
                
            } else if(Rcpp::is<Rcpp::NumericVector>(A) && Rcpp::is<Rcpp::NumericVector>(B)) {
                
                // return( BigDataStatMeth::Rcpp_vector_sum(A, B));
                
            } else {
                Rcpp::Rcout<<"\nData type not allowed";
            }    
        }
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
    }

    // return(Rcpp::wrap(C));
    return(R_NilValue);
}

