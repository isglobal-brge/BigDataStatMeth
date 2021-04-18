#ifndef svdBlockDecomposition_hdf5
#define svdBlockDecomposition_hdf5

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"
  #include "svdDecomposition.h"
  #include "svdutils.h"
  #include "rhdf5Utils.h"
  #include "matrix_utilities.h"
  #include "hdf5_normalization.h"

/*
  struct svdeig {
    Eigen::VectorXd d;
    Eigen::MatrixXd u;
    Eigen::MatrixXd v;
    bool bokuv;
    bool bokd;
  };
*/  
  
  //..// svdeig bdSVD_hdf5(Rcpp::RObject X);
  int First_level_SvdBlock_decomposition_hdf5(H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
                                              int irows, int icols, Rcpp::Nullable<int> threads);
  
  int Next_level_SvdBlock_decomposition_hdf5(H5File* file, std::string strGroupName, int k, int q, 
                                             bool bcenter, bool bscale, Rcpp::Nullable<int> threads);

#endif