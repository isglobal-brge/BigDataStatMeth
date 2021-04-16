#ifndef hdf5_to_Eigen
#define hdf5_to_Eigen


  #include <RcppEigen.h>
  #include "H5Cpp.h"
  #include "rhdf5Utils.h"

  // [[Rcpp::depends(RcppEigen)]]

  using namespace H5;
  using namespace Rcpp;
  
  Eigen::MatrixXd GetCurrentBlock_hdf5( H5File* file, DataSet* dataset,
                                      hsize_t offsetx, hsize_t offsety, 
                                      hsize_t countx, hsize_t county);
  
  Eigen::MatrixXd GetCurrentBlock_hdf5_Original( H5File* file, DataSet* dataset,
                                                 hsize_t offsetx, hsize_t offsety, 
                                                 hsize_t countx, hsize_t county);

  Eigen::MatrixXd RowMajorVector_to_ColMajorMatrix(double* datablock, int countx, int county);

#endif