#ifndef hdf5_invCholesky
#define hdf5_invCholesky

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    #include "BigDataStatMeth.h"
    #include "hdf5_to_Eigen.h"
    #include "hdf5_getDiagonalMatrix.h"
    #include "hdf5_writeDiagonalMatrix.h"
    #include "hdf5_writeOppositeTriangular.h"
    #include "pkg_omp.h"
    #include <iostream>
    #include <cstdlib>
    #include <chrono>
    #include <thread>
    #include <cmath>

    void Cholesky_decomposition_hdf5( H5File* file, DataSet* inDataset, DataSet* outDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads);
    void Inverse_of_Cholesky_decomposition_hdf5(  H5File* file, DataSet* InOutDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads);
    void Inverse_Matrix_Cholesky_parallel(  H5File* file, DataSet* InOutDataset, int idim0, int idim1, long dElementsBlock, Rcpp::Nullable<int> threads);
    // void Rcpp_bdInvCholesky_hdf5( H5File* file, DataSet* pdataset, std::string outgroup, std::string outdataset, bool bforce, Rcpp::Nullable<int> threads, Rcpp::Nullable<long> elementsBlock);
    void Rcpp_bdInvCholesky_hdf5( H5File* file, DataSet* pdataset, std::string outgroup, std::string outdataset,  bool bforce, bool fullMatrix,  Rcpp::Nullable<int> threads, Rcpp::Nullable<long> elementsBlock); 
    
    // R functions
    void bdInvCholesky_hdf5( std::string filename, std::string group, 
                             std::string dataset, std::string  outdataset, 
                             Rcpp::Nullable<std::string> outgroup, 
                             Rcpp::Nullable<bool> fullMatrix, 
                             Rcpp::Nullable<bool> force, 
                             Rcpp::Nullable<int> threads, 
                             Rcpp::Nullable<long> elementsBlock);
#endif
