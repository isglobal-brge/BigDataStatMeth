#ifndef hdf5_blockmultSparse
#define hdf5_blockmultSparse

   #include <RcppEigen.h>
   #include <Eigen/Sparse>
   #include "H5Cpp.h"
   #include "parallelBlockMult_hdf5.h"
   #include "utils_sparse.h"
   #include "hdf5_to_Eigen.h"
   #include "rhdf5Utils.h"


   // functions called from c++

   // functions called from R
   Rcpp::RObject blockmult_sparse_hdf5(std::string filename, const std::string group, 
                                       std::string A, std::string B,
                                       Rcpp::Nullable<std::string> outgroup );
#endif