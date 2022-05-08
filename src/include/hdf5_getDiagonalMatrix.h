#ifndef hdf5_getDiagonalMatrix
#define hdf5_getDiagonalMatrix

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    #include "hdf5_checks.h"
    #include "H5Cpp.h"

    using namespace Rcpp;
    
    // R functions
    Rcpp::RObject bdgetDiagonal_hdf5( std::string filename, std::string group, std::string dataset);

#endif
