#ifndef hdf5_QC
#define hdf5_QC

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  #include "hdf5_imputation.h"

  // functions called from c++
  void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
  void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
  int Remove_snp_low_data_HDF5( H5File* file, DataSet* dataset, bool bycols, std::string stroutdata, double pcent);

  // functions called from R
  Rcpp::RObject bdRemovelowdata( std::string filename, std::string group, std::string dataset, std::string outgroup, std::string outdataset, 
                                 Rcpp::Nullable<double> pcent = 0.5, Rcpp::Nullable<bool> bycols = true ); 

#endif
