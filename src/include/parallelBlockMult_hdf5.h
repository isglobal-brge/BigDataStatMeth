#ifndef parallelBlockMult_hdf5
#define parallelBlockMult_hdf5

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"
  #include "H5Cpp.h"
  #include "rhdf5Utils.h"
  #include <cstdlib>
  #include <cmath>
  #include <omp.h>
  #include <thread>
  

  //..// Eigen::MatrixXd block_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
  Eigen::MatrixXd RowMajorVector_to_ColMajorMatrix(double* datablock, int countx, int county);
  Eigen::MatrixXd hdf5_block_matrix_mul( IntegerVector sizeA, IntegerVector sizeB, 
                                         int hdf5_block, std::string filename, std::string strsubgroup );
  Eigen::MatrixXd Bblock_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, 
                                            int block_size, Rcpp::Nullable<int> threads );
  Eigen::MatrixXd Bblock_matrix_mul(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
  Eigen::MatrixXd Bblockmult(Rcpp::RObject a, Rcpp::RObject b, Rcpp::Nullable<int> block_size, 
                            Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads );

#endif