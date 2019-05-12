#ifndef olsgrid
#define olsgrid

  #include <RcppEigen.h>
  #include <iostream>
  #include <fstream>
  
  
  // dgemv_ symbol LAPACK-BLAS : performs one of the matrix-vector operations : y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
  extern "C" {
    extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
  }
  
  // dgemm_ symbol LAPACK-BLAS : performs one of the matrix-matrix operations : C := alpha*op( A )*op( B ) + beta*C,
  extern "C" {
    extern void dgemm_( char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int* );
  }
  
  // dtrsv_ symbol LAPACK-BLAS : solves one of the systems of equations : A*x = b,   or   A**T*x = b,
  extern "C" {
    extern void dtrsv_(char*, char*, char*, int*, double*, int*, double*, int*);	
  }
  
  // dtrsm_ symbol LAPACK-BLAS : solves one of the matrix equations : op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
  extern "C" {
    extern void dtrsm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
  }
  
  // dtrtri_ symbol LAPACK-BLAS Level 3 : computes the inverse of a real upper or lower triangular matrix A.
  extern "C" {
    extern void dtrtri_( char*, char*, int*, double*, int*, int*);
  }

#endif