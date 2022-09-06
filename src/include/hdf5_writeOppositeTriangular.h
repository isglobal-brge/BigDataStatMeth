#ifndef hdf5_writeOppositeTriangular
#define hdf5_writeOppositeTriangular
    
#include <RcppEigen.h>
#include "rhdf5Utils.h"
#include "hdf5_checks.h"
#include "H5Cpp.h"
    
    // C++ functions
    void Rcpp_setUpperTriangularMatrix( H5File* file, DataSet* pdataset, int dimensionSize, long dElementsBlock);
    void Rcpp_setLowerTriangularMatrix( H5File* file, DataSet* pdataset, int dimensionSize, long dElementsBlock);
    
    // R functions
    void bdWriteOppsiteTriangularMatrix_hdf5(std::string filename, std::string group, std::string dataset, 
                                             Rcpp::Nullable<bool> copytolower,
                                             Rcpp::Nullable<long> elementsBlock);
    
#endif
    