#ifndef hdf5_sortMatrices
#define hdf5_sortMatrices

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    #include "hdf5_checks.h"
    #include "H5Cpp.h"
    #include "vector_utilities.h"
    
    // [[Rcpp::depends(RcppEigen)]]
    using namespace Rcpp;
    using namespace std;
    
    
    // R functions
    void bdSort_hdf5_dataset( std::string filename, std::string group, 
                              std::string datasets, std::string outdataset, 
                              Rcpp::List blockedSortlist, std::string func, 
                              Rcpp::Nullable<std::string> outgroup, 
                              Rcpp::Nullable<bool> force);

#endif
