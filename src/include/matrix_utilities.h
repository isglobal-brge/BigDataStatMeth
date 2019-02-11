#ifndef matrix_utilities
#define matrix_utilities

  #include <RcppEigen.h>

  Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x);
  Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x);
  Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X );
  
  Rcpp::RObject Normalize_Data ( Rcpp::RObject & x );

#endif


