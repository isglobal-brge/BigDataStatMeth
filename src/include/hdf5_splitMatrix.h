#ifndef hdf5_splitMatrix
#define hdf5_splitMatrix

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  #include "hdf5_checks.h"
  #include "H5Cpp.h"

  // cpp functions
  int RcppSplit_matrix_hdf5 ( H5File* file, DataSet* dataset, bool bycols, std::string stroutgroup, std::string stroutdataset, int blocksize, int irows, int icols );

  // R functions
  void bdSplit_matrix_hdf5( std::string filename, std::string group, std::string dataset, 
                            Rcpp::Nullable<std::string> outgroup, Rcpp::Nullable<std::string> outdataset, 
                            Rcpp::Nullable<int> nblocks,  Rcpp::Nullable<int> blocksize,
                            Rcpp::Nullable<bool> bycols, Rcpp::Nullable<bool> force  );

#endif