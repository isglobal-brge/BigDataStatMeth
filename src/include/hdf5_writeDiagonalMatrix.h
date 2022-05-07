#ifndef hdf5_writeDiagonalMatrix
#define hdf5_writeDiagonalMatrix

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    #include "hdf5_checks.h"
    #include "H5Cpp.h"

    using namespace Rcpp;
    
    // R functions
    void bdWriteDiagonal_hdf5( Rcpp::RObject diagonal, std::string filename, std::string group, Rcpp::StringVector dataset);

#endif
