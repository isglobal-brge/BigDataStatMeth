#ifndef matrix_utilities
#define matrix_utilities

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"

  Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x);
  Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x);
  Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X, bool bc, bool bs );
  //..// Eigen::MatrixXd RcppNormalize_Data_hdf5 ( Eigen::MatrixXd  X, bool bc, bool bs );
  Eigen::MatrixXd RcppNormalize_Data_hdf5 ( Eigen::MatrixXd  X, bool bc, bool bs, bool btransp, Eigen::MatrixXd normdata );
  Rcpp::NumericMatrix RcppNormalize_Data_r ( Rcpp::NumericMatrix  x );
  
  Rcpp::RObject Normalize_Data ( Rcpp::RObject & x,
                                 Rcpp::Nullable<bool> bscale,
                                 Rcpp::Nullable<bool> bcenter);
  

#endif


