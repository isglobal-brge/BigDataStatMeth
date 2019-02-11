#ifndef parallel_MatrixVectorMultiplication
#define parallel_MatrixVectorMultiplication

  #include <RcppEigen.h>
  #include <RcppParallel.h>
  #include <cmath>

  Rcpp::NumericMatrix rcpp_xwxt(Rcpp::NumericMatrix mat, Rcpp::NumericVector w);
  Rcpp::NumericMatrix rcpp_xtwx(Rcpp::NumericMatrix mat, Rcpp::NumericVector w);
  Rcpp::NumericMatrix rcpp_parallel_Xy(Rcpp::NumericMatrix mat, Rcpp::NumericVector y);
  Rcpp::RObject xwxt(Rcpp::RObject mat, Rcpp::RObject w);
  Rcpp::RObject xtwx(Rcpp::RObject mat, Rcpp::RObject w);                   

#endif