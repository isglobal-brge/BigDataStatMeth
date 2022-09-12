#ifndef hdf5_WeightedProduct
#define hdf5_WeightedProduct

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    
    // C++ functions
    
    
    // R functions
    void bdWeightedProduct_hdf5( std::string filename, std::string group, std::string dataset,
                                 std::string vectorgroup, std::string vectordataset,
                                 std::string outdataset, 
                                 Rcpp::Nullable<std::string> outgroup,
                                 Rcpp::Nullable<bool> byrows,
                                 Rcpp::Nullable<int> force);
#endif
