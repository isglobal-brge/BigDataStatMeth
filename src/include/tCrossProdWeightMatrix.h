#ifndef tCrossProdWeightMatrix
#define tCrossProdWeightMatrix

    #include <RcppEigen.h>
    #include "BigDataStatMeth.h"
    #include "pkg_omp.h"
    #include <thread>
    
    // C++ functions :
    Eigen::MatrixXd Bblock_weighted_tcrossprod(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
    Eigen::MatrixXd Bblock_weighted_tcrossprod_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, 
                                                      int block_size, Rcpp::Nullable<int> threads  = R_NilValue);
    
    // R functions :
    Rcpp::RObject bdtCrossprod_Weighted(Rcpp::RObject a, Rcpp::RObject w, 
                                       Rcpp::Nullable<int> block_size, 
                                       Rcpp::Nullable<bool> paral,
                                       Rcpp::Nullable<int> threads );

#endif
