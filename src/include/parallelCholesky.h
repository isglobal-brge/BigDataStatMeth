#ifndef parallelCholesky
#define parallelCholesky

  #include <RcppEigen.h>
  #include <iostream>
  #include <cstdlib>
  #include <chrono>
  #include <omp.h>
  #include "ReadDelayedData.h"
  // #include </usr/lib/gcc/x86_64-linux-gnu/7/include/omp.h>
  #include <cmath>

  // Eigen::MatrixXd Cholesky_decomposition_parallel();
  Eigen::MatrixXd Cholesky_decomposition_parallel( Eigen::MatrixXd& A );
  Eigen::VectorXd Forward_Substituion_parallel(Eigen::MatrixXd L, Eigen::VectorXd y);
  Eigen::MatrixXd Inverse_Matrix_Cholesky_parallel( Eigen::MatrixXd L );
  Eigen::MatrixXd Inverse_of_Cholesky_decomposition_parallel( Eigen::MatrixXd& A, Eigen::MatrixXd& L );
  
#endif
