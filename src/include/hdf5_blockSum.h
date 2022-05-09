#ifndef hdf5_blockSum
#define hdf5_blockSum

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  
  
  
  
  
  // functions called from c++
  
  void hdf5_block_matrix_mul_hdf5_indatasets_transposed( std::string matA, std::string matB, 
                                                         IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                                         std::string filename, std::string strsubgroupIN, 
                                                         std::string strsubgroupINB, std::string strsubgroupOUT, 
                                                         std::string strdatasetOUT, 
                                                         int mem_block_size, bool bparal, bool browmajor, 
                                                         Rcpp::Nullable<int> threads );
  
  
  Rcpp::RObject blockSum_hdf5(std::string filename, const std::string group,
                                std::string A, 
                                const std::string groupB, 
                                std::string B,
                                Rcpp::Nullable<int> block_size, 
                                Rcpp::Nullable<bool> paral,
                                Rcpp::Nullable<int> threads,
                                Rcpp::Nullable<std::string> outgroup,
                                Rcpp::Nullable<std::string> outdataset);
#endif
