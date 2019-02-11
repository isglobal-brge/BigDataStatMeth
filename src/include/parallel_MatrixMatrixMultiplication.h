#ifndef prallel_MatrixMatrixMultiplication
#define prallel_MatrixMatrixMultiplication

  #include <RcppEigen.h>
  #include <RcppParallel.h>
  #include "matrix_utilities.h"

  Rcpp::NumericMatrix rcpp_parallel_XYProd(Rcpp::NumericMatrix matx, Rcpp::NumericMatrix maty);
  Rcpp::NumericMatrix rcpp_parallel_XYProdBlock(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY);
  Rcpp::NumericMatrix rcpp_parallel_XtYProd(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY);
  Rcpp::NumericMatrix rcpp_parallel_XYtProd(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY);
  


#endif
