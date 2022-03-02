#ifndef mlr_mr
#define mlr_mr

    #include <RcppEigen.h>
    //#include <omp.h>
    #include "BigDataStatMeth.h"
    #include "pkg_omp.h"
    #include <thread>
    #include "rhdf5Utils.h"
    #include "matrix_utilities.h"
    #include "parallelBlockMult_hdf5.h"
    #include "svdDecomposition.h"
    #include "qrdecomp.h"
    
    // functions called from c++
    Eigen::MatrixXd Rcpp_mlr_mr(Eigen::MatrixXd x, Eigen::MatrixXd y, int iblocks, Rcpp::Nullable<int> threads );
    //..// Rcpp::RObject Rcpp_mlr_mr(Eigen::MatrixXd x, Eigen::MatrixXd y, int iblocks, Rcpp::Nullable<int> threads );
    
    // functions called from R
    Rcpp::RObject bdMLR_MR(Rcpp::RObject X, Rcpp::RObject y, int blocks, Rcpp::Nullable<int> threads);

#endif
