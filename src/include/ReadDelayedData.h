#ifndef ReadDelayedData
#define ReadDelayedData

  #include <RcppEigen.h>
  #include "beachmat/numeric_matrix.h"  // To access numeric matrix
  #include "beachmat/integer_matrix.h"  // To access numeric matrix

  Eigen::MatrixXd read_DelayedArray_int( Rcpp::RObject A );
  Eigen::MatrixXd read_DelayedArray_real( Rcpp::RObject A );
  Eigen::MatrixXd read_DelayedArray( Rcpp::RObject A );

#endif