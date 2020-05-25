#ifndef ReadDelayedData
#define ReadDelayedData

  #include <RcppEigen.h>
  #include "beachmat/numeric_matrix.h"  // To access numeric matrix
  #include "beachmat/integer_matrix.h"  // To access numeric matrix
  #include "H5Cpp.h"
  #include "rhdf5Utils.h"

  Eigen::MatrixXd read_DelayedArray_int( Rcpp::RObject A );
  Eigen::MatrixXd read_DelayedArray_real( Rcpp::RObject A );
  Eigen::MatrixXd read_DelayedArray( Rcpp::RObject A );
  
  Eigen::Vector2i get_DelayedArray_size(Rcpp::RObject A);
  
  int write_DelayedArray_to_hdf5(H5std_string filename, const std::string CDatasetName, Rcpp::RObject A);
  int write_DelayedArray_int_hdf5( H5std_string filename, const std::string CDatasetName, Rcpp::RObject A );
  int write_DelayedArray_real_hdf5( H5std_string filename, const std::string CDatasetName, Rcpp::RObject A );
  
  Rcpp::NumericMatrix read_DelayedArray_int_r( Rcpp::RObject A );
  Rcpp::NumericMatrix read_DelayedArray_real_r( Rcpp::RObject A );
  Rcpp::NumericMatrix read_DelayedArray_rcpp( Rcpp::RObject A );

#endif