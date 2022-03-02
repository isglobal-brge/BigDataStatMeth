#ifndef hdf5_blockmult
#define hdf5_blockmult

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  #include "parallelBlockMult_hdf5.h"
  
  // functions called from c++
  
  // functions called from R
  Rcpp::RObject blockmult_hdf5(std::string filename, const std::string group, 
                               std::string A, 
                               const std::string groupB, 
                               std::string B,
                               Rcpp::Nullable<int> block_size, 
                               Rcpp::Nullable<bool> paral,
                               Rcpp::Nullable<int> threads,
                               Rcpp::Nullable<double> mixblock_size,
                               Rcpp::Nullable<std::string> outgroup,
                               Rcpp::Nullable<std::string> outdataset);
#endif
