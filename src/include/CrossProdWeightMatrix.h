#ifndef CrossProdWeightMatrix
#define CrossProdWeightMatrix
    
    #include <RcppEigen.h>
    #include "BigDataStatMeth.h"
    //#include <omp.h>
    #include "pkg_omp.h"
    #include <thread>
    #include "ReadDelayedData.h"
    
    // C++ functions :
    Eigen::MatrixXd Bblock_weighted_crossprod(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
    Eigen::MatrixXd Bblock_weighted_crossprod_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
                                                       int block_size, Rcpp::Nullable<int> threads  = R_NilValue);
    
    // R functions :
    Rcpp::RObject bdCrossprod_Weighted(Rcpp::RObject a, Rcpp::RObject w, 
                                       Rcpp::Nullable<int> block_size, 
                                       Rcpp::Nullable<bool> paral,
                                       Rcpp::Nullable<int> threads );

#endif
