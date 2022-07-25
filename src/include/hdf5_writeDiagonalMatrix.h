#ifndef hdf5_writeDiagonalMatrix
#define hdf5_writeDiagonalMatrix

    #include <RcppEigen.h>
    // #include "BigDataStatMeth.h"
    #include "rhdf5Utils.h"
    #include "hdf5_checks.h"
    #include "H5Cpp.h"

    // C++ functions
    void Rcpp_setDiagonalMatrix( H5File* file, DataSet* pdataset, Rcpp::NumericVector intNewDiagonal);
    
    // R functions
    void bdWriteDiagonal_hdf5( Rcpp::RObject diagonal, std::string filename, std::string group, Rcpp::StringVector dataset);

#endif
