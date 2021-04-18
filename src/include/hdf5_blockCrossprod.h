#ifndef hdf5_blocCrossprod
#define hdf5_blocCrossprod

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  #include "parallelBlockMult_hdf5.h"
  
  // functions called from c++
  // int hdf5_block_matrix_crossprod_hdf5( std::string matA, IntegerVector sizeA, int hdf5_block,
  //                                       std::string filename, std::string strsubgroupIN, std::string strsubgroupOUT, 
  //                                       int mem_block_size, bool bparal, bool browmajor, 
  //                                       Rcpp::Nullable<int> threads);
  // int hdf5_block_matrix_crossprod_hdf5( std::string matA, IntegerVector sizeA, 
  //                                       std::string matB, IntegerVector sizeB, 
  //                                       int hdf5_block,
  //                                       std::string filename, std::string strsubgroupIN, std::string strsubgroupOUT, 
  //                                       int mem_block_size, bool bparal, bool browmajor, 
  //                                       Rcpp::Nullable<int> threads);
  
  int hdf5_block_matrix_crossprod_hdf5( std::string matA, IntegerVector sizeA, 
                                        std::string matB, IntegerVector sizeB, 
                                        int hdf5_block,
                                        std::string filename, 
                                        std::string strsubgroupIN, std::string strsubgroupINB, 
                                        std::string strsubgroupOUT, 
                                        int mem_block_size, bool bparal, bool browmajor, 
                                        Rcpp::Nullable<int> threads);
    
  // functions called from R
  // Rcpp::RObject blockCrossprod_hdf5(std::string filename, const std::string group, 
  //                                   std::string A, std::string B,
  //                                   Rcpp::Nullable<int> block_size, 
  //                                   Rcpp::Nullable<bool> paral,
  //                                   Rcpp::Nullable<int> threads,
  //                                   Rcpp::Nullable<double> mixblock_size,
  //                                   Rcpp::Nullable<std::string> outgroup);
  
  Rcpp::RObject Crossprod_hdf5(std::string filename, const std::string group,
                               std::string A, 
                               Rcpp::Nullable<std::string> groupB, 
                               Rcpp::Nullable<std::string>  B, 
                               Rcpp::Nullable<int> block_size, 
                               Rcpp::Nullable<bool> paral,
                               Rcpp::Nullable<int> threads,
                               Rcpp::Nullable<double> mixblock_size,
                               Rcpp::Nullable<std::string> outgroup);
    
#endif