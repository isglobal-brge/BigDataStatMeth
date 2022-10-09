#include "include/solveMatrixEquation.h"


//' Solve matrix equations
//' 
//' This function solve matrix equations 
//'  \code{A * X = B } 
//' where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//' 
//' @param A numerical matrix. 
//' @param B numerical matrix.
//' @return X numerical matrix. 
//' @examples
//' 
//' library(BigDataStatMeth)
//' 
//' n <- 500
//' m <- 500
//' 
//' # R Object
//' 
//' A <- matrix(runif(n*m), nrow = n, ncol = m)
//' B <- matrix(runif(n), nrow = n)
//' AS <- A%*%t(A)
//'       
//' X <- bdSolve(A, B)
//' XR <- solve(A,B)
//'       
//' all.equal(X, XR, check.attributes=FALSE)
//'   
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSolve(const Rcpp::RObject A, const Rcpp::RObject B) 
{
   
   ////// IMPORTANT !!!! :  op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   //                              where alpha is a scalar, X and B are m by n matrices
   
   auto dmtypeA = beachmat::find_sexp_type(A);
   auto dmtypeB = beachmat::find_sexp_type(B);
   
   Eigen::MatrixXd a; // = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(A);
   Eigen::MatrixXd b; // = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(B);

   char Uchar = 'U';
   int info = 0;
   
   
   try {
      
      if ( dmtypeA == INTSXP || dmtypeA==REALSXP ) {
         if ( A.isS4() == true){
            // a = read_DelayedArray(A);
            throw("Only numeric matrix allowd");
         }else {
            a = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(A);
         }
      } else {
         throw std::runtime_error("unacceptable matrix type - matrix A");
      }
      
      if ( dmtypeB == INTSXP || dmtypeB==REALSXP ) {
         if ( B.isS4() == true){
            // b = read_DelayedArray(B);
            throw("Only numeric matrix allowd");
         }else {
            b = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(B);
         }
      } else {
         throw std::runtime_error("unacceptable matrix type - matrix B");
      }
      
      // Declare matrix variables
      int n = a.rows();
      int nrhs = b.cols();
      
      std::vector<int> ipiv(n);
      //..// int ipiv[n];
      
      int lwork = std::max( 1, n );
      std::vector<double> work(lwork);
      int lda = std::max( 1, n );
      int ldb = std::max( 1, n );
      
      
      // Solve matrix equation
      if( a == a.transpose()  )
      {
         // dsysv_( char* UPLO, int* N , int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, double* WORK, int* LWORK, int* INFO);
         dsysv_( & Uchar, &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, work.data(), &lwork, &info);
      } else {
         
         // dgesv( int N, int NRHS, double A, int LDA, int IPIV, double B, int LDB, int INFO);
         dgesv_( &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, &info );
      }
      
   } catch(std::exception &ex) {
      Rcpp::Rcout<< ex.what();
      return Rcpp::wrap(-1);
   }
    
    return(Rcpp::wrap(b));

}





/***R

*/
