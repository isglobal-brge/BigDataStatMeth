#ifndef hdf5_removeMaf
#define hdf5_removeMaf

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
#include "hdf5_getAlleleFreq.h"
#include "hdf5_QC.h"


  // functions called from c++
  int Remove_MAF_HDF5( H5File* file, DataSet* dataset, bool bycols, std::string stroutdata, double pcent);


  // functions called from R

#endif