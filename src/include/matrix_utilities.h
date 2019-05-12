#ifndef matrix_utilities
#define matrix_utilities

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"

  Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x);
  Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x);
  Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X );
  Rcpp::NumericMatrix RcppNormalize_Data_r ( Rcpp::NumericMatrix  x );
  
  Rcpp::RObject Normalize_Data ( Rcpp::RObject & x );
  

#endif


