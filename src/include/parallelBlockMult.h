#ifndef parallelBlockMult
#define parallelBlockMult

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"
  #include <mpi.h>
  #include <cstdlib>
  #include <cmath>
  #include <omp.h>

  Eigen::MatrixXd block_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
  Eigen::MatrixXd block_matrix_mul(Eigen::MatrixXd& A, Eigen::MatrixXd& B, int block_size);

#endif