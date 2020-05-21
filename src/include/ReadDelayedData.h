#ifndef ReadDelayedData
#define ReadDelayedData

  #include <RcppEigen.h>
  #include "beachmat/numeric_matrix.h"  // To access numeric matrix
  #include "beachmat/integer_matrix.h"  // To access numeric matrix

  Eigen::MatrixXd read_DelayedArray_int( Rcpp::RObject A );
  Eigen::MatrixXd read_DelayedArray_real( Rcpp::RObject A );
  Eigen::MatrixXd read_DelayedArray( Rcpp::RObject A );
  
  Eigen::Vector2i get_DelayedArray_size(Rcpp::RObject A);
  
  Rcpp::NumericMatrix read_DelayedArray_int_r( Rcpp::RObject A );
  Rcpp::NumericMatrix read_DelayedArray_real_r( Rcpp::RObject A );
  Rcpp::NumericMatrix read_DelayedArray_rcpp( Rcpp::RObject A );

#endif