#ifndef qrdecomp
#define qrdecomp

   #include <RcppEigen.h>
   #include "pseudoinv.h"
   #include "parallelBlockMult.h"
   #include "ReadDelayedData.h"
  
  
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
   
   
   // dgesvd_ is a symbol in the LAPACK-BLAS Level 3 
   //    DGESVD computes the singular value decomposition (SVD) of a real M-by-N matrix A, 
   //       optionally computing the left and/or right singular vectors
   extern "C" {
   extern void dgesvd_( char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
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
   
   strQR rcpp_bdQR( Eigen::MatrixXd & A, bool bthin);


#endif