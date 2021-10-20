#ifndef hdf5_lm
#define hdf5_lm

  #include <RcppEigen.h>
  #include <omp.h>
  #include <thread>
  #include "rhdf5Utils.h"
  #include "matrix_utilities.h"
  #include "parallelBlockMult_hdf5.h"
  #include "svdDecomposition.h"
  #include "qrdecomp.h"
  #include "mlr_mr.h"

  // functions called from c++
  Eigen::MatrixXd Rcpp_mlr_mr_hdf5(Eigen::MatrixXd x, Eigen::MatrixXd y, int iblocks, Rcpp::Nullable<int> threads );
  //..// Rcpp::RObject Rcpp_mlr_mr(Eigen::MatrixXd x, Eigen::MatrixXd y, int iblocks, Rcpp::Nullable<int> threads );
  
  // functions called from R
  Rcpp::RObject bdMLR_MR_hdf5(std::string filename, 
                              const std::string group, 
                              std::string dataset,
                              const std::string betasgroup, 
                              std::string betasdataset,
                              int blocks, 
                              Rcpp::Nullable<std::string> outgroup,
                              Rcpp::Nullable<std::string> outdataset,
                              Rcpp::Nullable<int> threads) ;



#endif