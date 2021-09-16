#ifndef hdf5_reduceMatrix
#define hdf5_reduceMatrix

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  #include "hdf5_checks.h"
  #include "H5Cpp.h"

  // cpp functions
  int RcppReduce_matrix_hdf5 ( H5File* file,  std::string strgroup, std::string stroutgroup, std::string stroutdataset, std::string strreducefunction, bool bremove );

  // R functions
  Rcpp::RObject bdReduce_matrix_hdf5( std::string filename, std::string group, 
                                      std::string reducefunction, 
                                      Rcpp::Nullable<std::string> outgroup, Rcpp::Nullable<std::string> outdataset,
                                      Rcpp::Nullable<bool> force , Rcpp::Nullable<bool> remove );

#endif