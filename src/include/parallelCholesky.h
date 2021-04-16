#ifndef parallelCholesky
#define parallelCholesky

  #include <RcppEigen.h>
  #include <iostream>
  #include <cstdlib>
  #include <chrono>
  #include <omp.h>
  #include <thread>
  #include "ReadDelayedData.h"
  #include <cmath>

  Eigen::MatrixXd Cholesky_decomposition_parallel( Eigen::MatrixXd& A, Rcpp::Nullable<int> threads );
  Eigen::VectorXd Forward_Substituion_parallel(Eigen::MatrixXd L, Eigen::VectorXd y, Rcpp::Nullable<int> threads);
  Eigen::MatrixXd Inverse_Matrix_Cholesky_parallel( Eigen::MatrixXd L, Rcpp::Nullable<int> threads );
  Eigen::MatrixXd Inverse_of_Cholesky_decomposition_parallel( Eigen::MatrixXd& A, Eigen::MatrixXd& L, Rcpp::Nullable<int> threads );
  
#endif
