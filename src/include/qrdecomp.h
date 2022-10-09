#ifndef qrdecomp
#define qrdecomp

    #include <RcppEigen.h>
    //#include <omp.h>
    #include "pkg_omp.h"
    #include "pseudoinv.h"
    #include "beachmat/numeric_matrix.h"  // To access numeric matrix
    #include "beachmat/integer_matrix.h"  // To access numeric matrix
    #include "parallelBlockMult_hdf5.h"
  
  
   // Symbols in the LAPACK library files : 
   extern "C" {
   extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
   }
   
   // dgemv_ is a symbol in the LAPACK-BLAS library files 
   //    DGEMV  performs one of the matrix-vector operations : y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   extern "C" {
   extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
   }
   
   // dgemm_ is a symbol in the LAPACK-BLAS library files 
   //    DGEMM  performs one of the matrix-matrix operations : C := alpha*op( A )*op( B ) + beta*C,
   extern "C" {
   extern void dgemm_( char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int* );
   }
   
   // dtrsv_ is a symbol in the LAPACK-BLAS library files - 
   //    DTRSV  solves one of the systems of equations : A*x = b,   or   A**T*x = b,
   extern "C" {
   extern void dtrsv_(char*, char*, char*, int*, double*, int*, double*, int*);	
   }

   
   // dtrsm_ is a symbol in the LAPACK-BLAS library files - 
   //    DTRSM  solves one of the matrix equations : op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   extern "C" {
   extern void dtrsm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
   }
   
   // dtrtri_ is a symbol in the LAPACK-BLAS Level 3 
   //    DTRTRI computes the inverse of a real upper or lower triangular matrix A.
   extern "C" {
   extern void dtrtri_( char*, char*, int*, double*, int*, int*);
   }
   
   // dgeqrf_ is a symbol in the LAPACK-BLAS Level 3 
   //    DGEQRF computes a QR factorization of a real M-by-N matrix A
   extern "C" {
   extern void dgeqrf_( int*, int*, double*, int*, double*, double*, int*, int*);
   }
   
   
   
   // DGEQP3_ is computes a QR factorization with column pivoting of a matrix A:  A*P = Q*R  using Level 3 BLAS.
   extern "C" {
      extern void dgeqp3_( int*, int*, double*, int*, int*, double*, double*, int*, int*);
   }
   
   // dscal_ is a symbol in the LAPACK-BLAS Level 3 
   //    DSCAL scales a vector by a constant.
   extern "C" {
   extern void dscal_( int*, double*, double*, int*);
   }
   
   struct strQR {
      Eigen::MatrixXd Q;
      Eigen::MatrixXd R;
   };
   
   
   // C++ Functions
   strQR rcpp_bdQR( Eigen::MatrixXd & A, bool bthin);
   // strQR rcpp_bdQR_parallel( Eigen::MatrixXd A, Rcpp::Nullable<int> threads);
   
   // R Functions
   Rcpp::RObject bdQR( const Rcpp::RObject & X, Rcpp::Nullable<bool> thin);
   // Rcpp::RObject bdQR_compact( const Rcpp::RObject & A);
   Rcpp::RObject bddtrsm(Rcpp::RObject R, Rcpp::RObject Z, Rcpp::Nullable<int> threads) ; 


#endif
