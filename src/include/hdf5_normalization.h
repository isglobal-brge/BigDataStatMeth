#ifndef hdf5_normalization
#define hdf5_normalization

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  #include "hdf5_getSDandMean.h"
  
  // functions called from c++
  Eigen::MatrixXd RcppNormalize_Data_hdf5 ( Eigen::MatrixXd  X, bool bc, bool bs, bool btransp, Eigen::MatrixXd normdata );
  Eigen::MatrixXd RcppNormalize_Data_R_hdf5( Eigen::MatrixXd  X, bool bc, bool bs, bool btransp, Eigen::MatrixXd normdata);
  
  // functions called from R
  /**
  Rcpp::RObject Normalize_hdf5(std::string filename, std::string group, std::string dataset, 
                               Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale,
                               Rcpp::Nullable<int> wsize);
  **/
#endif
