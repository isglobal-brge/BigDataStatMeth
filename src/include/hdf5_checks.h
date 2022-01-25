#ifndef hdf5_checks
#define hdf5_checks

  #include <RcppEigen.h>
  #include "rhdf5Utils.h"
  #include "H5Cpp.h"

  bool exist_FileGroupDataset(std::string filename, std::string group, std::string dataset);
  double prepare_outDataset(H5File* file, std::string outDataset, bool bforce);
  
  double prepare_outGroup(H5File* file, std::string outGroup, bool bforce);
    
  

#endif