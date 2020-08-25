#ifndef pseudoinv
#define pseudoinv

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"


  // dgemm_ is a symbol in the LAPACK-BLAS library files 
  //    DGEMM  performs one of the matrix-matrix operations : C := alpha*op( A )*op( B ) + beta*C,
  extern "C" {
    extern void dgemm_( char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int* );
  }
  
  
  // dgesvd_ is a symbol in the LAPACK-BLAS Level 3 
  //    DGESVD computes the singular value decomposition (SVD) of a real M-by-N matrix A, 
  //       optionally computing the left and/or right singular vectors
  extern "C" {
    extern void dgesvd( char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
  }
  
  // dscal_ is a symbol in the LAPACK-BLAS Level 3 
  //    DSCAL scales a vector by a constant.
  extern "C" {
    extern void dscal_( int*, double*, double*, int*);
  }
  
  
  Eigen::MatrixXd rcpp_bdpseudoinv(Eigen::MatrixXd A);



#endif