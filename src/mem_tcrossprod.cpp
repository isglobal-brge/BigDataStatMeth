#include <BigDataStatMeth.hpp>
#include "memAlgebra/memOptimizedProducts.hpp"
#include "memAlgebra/memMultiplication.hpp"


//' Transpodsed Crossproduct 
//' 
//' This function performs a transposed crossproduct of a numerical matrix.
//' 
//' @export
//' 
//' @param A numerical matrix
//' @param B optional, numerical matrix
//' @param block_size (optional, defalut = NULL) block size to make matrix 
//' multiplication, if `block_size = 1` no block size is applied 
//' (size 1 = 1 element per block) if `block_size = NULL` (default) optimum 
//' block size is computed
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel 
//' computation else if paral = FALSE performs serial computation
//' @param threads (optional) only if bparal = true, number of concurrent threads
//' in parallelization if threads is null then threads =  maximum number of 
//' threads available
//' @return numerical matrix with transposed crossproduct
//' @examples
//' 
//' n <- 100
//' p <- 60
//' 
//' X <- matrix(rnorm(n*p), nrow=n, ncol=p)
//' res <- bdtCrossprod(X)
//' 
//' all.equal(tcrossprod(X), res)
//' 
//' n <- 100
//' p <- 100
//' 
//' Y <- matrix(rnorm(n*p), nrow=n)
//' 
//' # With two matrices
//' res <- bdtCrossprod(X,Y)
//' 
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd bdtCrossprod( Rcpp::RObject A, Rcpp::Nullable<Rcpp::RObject> B =  R_NilValue, 
                             Rcpp::Nullable<bool> transposed = R_NilValue,
                             Rcpp::Nullable<int> block_size = R_NilValue, 
                             Rcpp::Nullable<bool> paral = R_NilValue,
                             Rcpp::Nullable<int> threads = R_NilValue )
{
    
    
    Eigen::MatrixXd C;
    
    try {
        
        Eigen::MatrixXd mA;
        Eigen::MatrixXd mB;
        bool bparal;
        
        if( paral.isNull()) {
            bparal = false;
        } else {
            bparal = Rcpp::as<bool> (paral);
        }
        
        // Read DelayedmArray's A and b
        if ( Rcpp::is<Rcpp::NumericMatrix>(A) || Rcpp::is<Rcpp::IntegerMatrix>(A))    
        {
            try{  
                mA = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(A);
            } catch(std::exception &ex) { }
            
        } else {
            throw("Matrix A is not numeric - Only numeric matrix allowed");
            
        }
        
        if(B.isNull()) {
            C = BigDataStatMeth::bdtcrossproduct(mA);
        } else {
            
            
            if(Rcpp::is<Rcpp::NumericMatrix>(B) || Rcpp::is<Rcpp::IntegerMatrix>(B)) {
                try{  
                    mB = Rcpp::as<Eigen::MatrixXd>(B); 
                }
                catch(std::exception &ex) { }
            } else {
                throw("Matrix B is not numeric - Only numeric matrix allowed");
            }
            
            Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mTrans(mB.data(), mB.cols(), mB.rows());
            
            if(bparal == true) {
                C = BigDataStatMeth::Rcpp_block_matrix_mul_parallel(mA, mTrans, block_size, threads);
                
            } else if (bparal == false)  {
                C = BigDataStatMeth::Rcpp_block_matrix_mul(mA, mTrans, block_size);
            }
            
        }
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
    }
    
    return(C);
    
}
