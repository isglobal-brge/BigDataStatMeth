#ifndef hdf5_pca
#define hdf5_pca

  #include <RcppEigen.h>
  #include "H5Cpp.h"
  #include "matrix_utilities.h"
  #include "svdDecomposition.h"
  #include "optimizedproduct.h"
  #include "hdf5_to_Eigen.h"
  #include "rhdf5Utils.h"
  #include "parallelBlockMult_hdf5.h"
  

  struct pcaeig {
    Eigen::MatrixXd varCoord;
    Eigen::MatrixXd Y;
    Eigen::MatrixXd var;
    Eigen::MatrixXd percvar;
    Eigen::MatrixXd components;
    std::string hdf5file = "";
  };
  
  
  int get_HDF5_PCA_variance_ptr(  H5File* file, std::string strdataset);
  // Rcpp::RObject bdPCA_hdf5(std::string filename, std::string strsubgroup, std::string strdataset);
  //..// Rcpp::RObject bdPCA_hdf5(std::string filename, std::string strsubgroup, std::string strdataset, Rcpp::Nullable<int> threads = R_NilValue);
  
  //..//Rcpp::RObject bdPCA_hdf5(std::string filename, std::string group, std::string dataset, Rcpp::Nullable<int> threads = R_NilValue);
  Rcpp::RObject bdPCA_hdf5(std::string filename, std::string group, std::string dataset, 
                           Rcpp::Nullable<int> bcenter = false, Rcpp::Nullable<int> bscale = false, 
                           Rcpp::Nullable<int> threads = R_NilValue);


#endif