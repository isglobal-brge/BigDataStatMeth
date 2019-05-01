#ifndef parallelBlockMult
#define parallelBlockMult

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"
  #include <cstdlib>
  #include <cmath>
  #include <omp.h>

  Eigen::MatrixXd block_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
  Eigen::MatrixXd block_matrix_mul(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
  Eigen::MatrixXd blockmult(Rcpp::RObject a, Rcpp::RObject b, Rcpp::Nullable<int> block_size, Rcpp::Nullable<bool> paral);

#endif