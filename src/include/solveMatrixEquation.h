#ifndef solveMatrixEquation
#define solveMatrixEquation

   #include <RcppEigen.h>
   #include "ReadDelayedData.h"


   // Symbols in the LAPACK library files : 

   // DGESV computes the solution to a real system of linear equations : A * X = B, where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
   extern "C" {
      extern int dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
   }

   // DSYSV computes the solution to a real system of linear equations :   A * X = B, where A is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices.
   extern "C" {
      extern int dsysv_( char*, int*, int*, double*, int*, int*, double*, int*, double*, int*, int*);
   }
   
         
   // R functions
   Rcpp::RObject bdSolve(const Rcpp::RObject A, const Rcpp::RObject B);


#endif
