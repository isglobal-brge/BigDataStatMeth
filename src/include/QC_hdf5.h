#ifndef QC_hdf5
#define QC_hdf5

  #include <RcppEigen.h>
  #include "H5Cpp.h"
  #include "rhdf5Utils.h"
  #include "hdf5_imputation.h"

  // functions called from c++
  int Remove_snp_low_data_HDF5( H5File* file, DataSet* dataset, bool bycols, std::string stroutdata, double pcent);

  // functions called from R
  void bdRemovelowdata( std::string filename, std::string group, std::string dataset, std::string outgroup, std::string outdataset, 
                           Rcpp::Nullable<double> pcent = 0.5, Rcpp::Nullable<bool> SNPincols = true );

#endif