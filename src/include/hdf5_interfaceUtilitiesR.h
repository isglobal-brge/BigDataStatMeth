#ifndef hdf5_interfaceUtilitiesR
#define hdf5_interfaceUtilitiesR

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    #include "hdf5_checks.h"
    #include "H5Cpp.h"


    Rcpp::RObject bdgetDatasetsList_hdf5(std::string filename, std::string strgroup, Rcpp::Nullable<std::string> pref);


#endif