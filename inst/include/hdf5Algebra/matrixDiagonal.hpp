#ifndef BIGDATASTATMETH_HDF5_MATRIXDIAGONAL_HPP
#define BIGDATASTATMETH_HDF5_MATRIXDIAGONAL_HPP

// #include <RcppEigen.h>
// #include "H5Cpp.h"

namespace BigDataStatMeth {


extern inline Rcpp::NumericVector getDiagonalfromMatrix( BigDataStatMeth::hdf5Dataset* dsMat)
{
    
    Rcpp::NumericVector intNewDiagonal;
    
    try {
        
        std::vector<hsize_t> offset = {0, 0},
                             count = {1, 1},
                             stride = {1, 1},
                             block = {1, 1};
        
        Rcpp::NumericVector diagElem(1);
        
        for(hsize_t i=0; i < dsMat->nrows(); i++) {
            offset[0] = i; offset[1] = i;
            dsMat->readDatasetBlock(  offset, count, stride, block,  REAL(diagElem) );
            intNewDiagonal.push_back(diagElem[0]);
        }
        
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdgetDiagonal_hdf5" << ex.what();
        return(intNewDiagonal);
    }
    
    return(intNewDiagonal);
}



extern inline void setDiagonalMatrix( BigDataStatMeth::hdf5Dataset* dsMat, Rcpp::NumericVector intNewDiagonal)
{
    
    try{
        
        std::vector<hsize_t> count = {1, 1},
                             stride = {1, 1},
                             block = {1, 1};
        
        for( hsize_t i=0; i < intNewDiagonal.size(); i++) {
            std::vector<hsize_t> offset = {i, i};
            dsMat->writeDatasetBlock(Rcpp::wrap(intNewDiagonal(i)), offset, count, stride, block, true);
        }
        
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdgetDiagonal_hdf5" << ex.what();
        return void();
    }
    
    return void();
}



}

#endif // BIGDATASTATMETH_HDF5_MATRIXDIAGONAL_HPP

