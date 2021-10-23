#ifndef hdf5_bindMatrices
#define hdf5_bindMatrices

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    #include "hdf5_checks.h"
    #include "H5Cpp.h"
    
    // R functions
    void bdBind_hdf5( std::string filename, std::string group, Rcpp::StringVector datasets, 
                               std::string outgroup, std::string outdataset, std::string func,
                               Rcpp::Nullable<bool> force );

#endif