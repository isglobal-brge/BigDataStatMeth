#ifndef utils_sparse
#define utils_sparse

   #include <RcppEigen.h>
   #include "hdf5_to_Eigen.h"
   #include "rhdf5Utils.h"

   bool is_hdf5_readed_block_sparse(H5File* file, DataSet* dataset,
                                    hsize_t offsetx, hsize_t offsety, 
                                    hsize_t countx, hsize_t county); // Read data from hdf5 file and Test if loaded block is sparse

   bool is_sparse( Eigen::MatrixXd mat) ; // Test if eigen matrix is sparse


#endif