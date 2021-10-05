#ifndef hdf5_applyFunction
#define hdf5_applyFunction

#include <RcppEigen.h>
#include "rhdf5Utils.h"
#include "qrdecomp.h"
#include "optimizedproduct.h"
#include "svdDecomposition.h"
#include "hdf5_checks.h"
#include "H5Cpp.h"

// // cpp functions
// int RcppApply_Function_hdf5 ( H5File* file, DataSet* dataset, std::string stroutgroup, std::string stroutdataset);

// R functions
Rcpp::RObject bdapply_Function_hdf5( std::string filename, std::string group, Rcpp::RObject datasets, 
                                     std::string outgroup, std::string func,
                                     Rcpp::Nullable<bool> force );

#endif