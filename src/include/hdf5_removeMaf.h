#ifndef hdf5_removeMaf
#define hdf5_removeMaf

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  #include "hdf5_getAlleleFreq.h"
  #include "hdf5_QC.h"



  // functions called from c++
  //..// int Remove_MAF_HDF5( H5File* file, DataSet* dataset, bool bycols, std::string stroutdata, double pcent);
  int Remove_MAF_HDF5( H5File* file, DataSet* dataset, bool bycols, std::string stroutdata, double pcent, int blocksize);


  // functions called from R
  Rcpp::RObject bdremove_maf_hdf5( std::string filename, std::string group, std::string dataset, std::string outgroup, std::string outdataset, 
                                   Rcpp::Nullable<double> maf = 0.05, Rcpp::Nullable<bool> bycols = false, Rcpp::Nullable<int> blocksize = 100 );
#endif