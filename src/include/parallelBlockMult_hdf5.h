#ifndef parallelBlockMult_hdf5
#define parallelBlockMult_hdf5

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"
  #include "H5Cpp.h"
  #include "rhdf5Utils.h"
  #include "hdf5_to_Eigen.h"
  #include <cstdlib>
  #include <cmath>
  #include <omp.h>
  #include <thread>
  

  
  /*** TODO : 
      Simplify all this functions 
  ***/  
   
  //..// Eigen::MatrixXd block_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
  //....// Eigen::MatrixXd RowMajorVector_to_ColMajorMatrix(double* datablock, int countx, int county);
  Eigen::MatrixXd hdf5_block_matrix_mul( IntegerVector sizeA, IntegerVector sizeB, 
                                         int hdf5_block, std::string filename, std::string strsubgroup );
  Eigen::MatrixXd Bblock_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, 
                                            int block_size, Rcpp::Nullable<int> threads );
  Eigen::MatrixXd Bblock_matrix_mul(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size);
  
  int hdf5_block_matrix_mul_hdf5( IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                  std::string filename, std::string strsubgroupIN, std::string strsubgroupOUT,
                                  int mem_block_size, bool bparal, Rcpp::Nullable<int> threads);

  int hdf5_block_matrix_mul_hdf5_transposed( IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                             std::string filename, std::string strsubgroupIN, std::string strsubgroupOUT, 
                                             int mem_block_size, bool bparal, bool browmajor, 
                                             Rcpp::Nullable<int> threads);
  
  // Functions : Matrix multiplication directly from file
  int hdf5_block_matrix_mul_hdf5_indatasets_transposed( std::string matA, std::string matB, 
                                                      IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                                      std::string filename, std::string strsubgroupIN, std::string strsubgroupOUT, 
                                                      int mem_block_size, bool bparal, bool browmajor, 
                                                      Rcpp::Nullable<int> threads );
  
  int hdf5_block_matrix_mul_hdf5_indatasets( std::string matA, std::string matB,
                                             IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                             std::string filename, std::string strsubgroupIN, std::string strsubgroupOUT, 
                                             int mem_block_size, bool bparal, bool browmajor, 
                                             Rcpp::Nullable<int> threads);
    
  /*Rcpp::List Bblockmult(Rcpp::RObject a, Rcpp::RObject b, Rcpp::Nullable<int> block_size, 
                             Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads,
                             Rcpp::Nullable<double> bigmatrix, Rcpp::Nullable<std::string> outfile);*/
  
  Rcpp::List blockmult(Rcpp::RObject a, Rcpp::RObject b, 
                        Rcpp::Nullable<int> block_size,
                        Rcpp::Nullable<bool> paral,
                        Rcpp::Nullable<int> threads,
                        Rcpp::Nullable<double> bigmatrix,
                        Rcpp::Nullable<std::string> outfile,
                        Rcpp::Nullable<double> mixblock_size,
                        Rcpp::Nullable<bool> onmemory);
  
#endif