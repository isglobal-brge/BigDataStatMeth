#ifndef matrix_utilities
#define matrix_utilities

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"
  #include "rhdf5Utils.h"
  #include "hdf5_to_Eigen.h"
  #include <regex>


  // functions called from c++
  Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x);
  Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x);
  Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X, bool bc, bool bs );
  //..// Eigen::MatrixXd RcppNormalize_Data_hdf5 ( Eigen::MatrixXd  X, bool bc, bool bs );
  /**Eigen::MatrixXd RcppNormalize_Data_hdf5 ( Eigen::MatrixXd  X, bool bc, bool bs, bool btransp, Eigen::MatrixXd normdata );**/
  Rcpp::NumericMatrix RcppNormalize_Data_r ( Rcpp::NumericMatrix  x );
  
  Eigen::MatrixXd GetCurrentBlock( Eigen::MatrixXd X, int startrow, int startcol, int nrows, int ncols );
  
  
  // functions called from R
  Rcpp::RObject bdNormalize_Data(Rcpp::RObject & X, Rcpp::Nullable<bool> bscale, Rcpp::Nullable<bool> bcenter);
  
  /**Rcpp::RObject Normalize_hdf5(std::string filename, std::string group, std::string dataset, 
                               Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale,
                               Rcpp::Nullable<int> wsize); **/
  

#endif


