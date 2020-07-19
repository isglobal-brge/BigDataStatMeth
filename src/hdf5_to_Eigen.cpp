#include "include/hdf5_to_Eigen.h"


// Convert array with readed data in rowmajor to colmajor matrix : 
//    Eigen and R works with ColMajor and hdf5 in RowMajor
Eigen::MatrixXd RowMajorVector_to_ColMajorMatrix(double* datablock, int countx, int county)
{
  Eigen::MatrixXd mdata(countx, county);
  
  for (size_t i=0; i<countx; i++)
    for(size_t j=0;j<county;j++)
      mdata(i,j) = datablock[i*county+j];
  
  return(mdata);
}


// Read block from hdf5 matrix
Eigen::MatrixXd GetCurrentBlock_hdf5( H5File* file, DataSet* dataset,
                                      hsize_t offsetx, hsize_t offsety, 
                                      hsize_t countx, hsize_t county)
{
  
  IntegerVector offset = IntegerVector::create(offsetx, offsety) ;
  IntegerVector count = IntegerVector::create(countx, county) ;
  IntegerVector stride = IntegerVector::create(1, 1) ;
  IntegerVector block = IntegerVector::create(1, 1) ;
  
  NumericMatrix data(countx, county);
  
  // read_HDF5_matrix_subset(filename, dataset, offset, count, stride, block, REAL(data));
  read_HDF5_matrix_subset(file, dataset, offset, count, stride, block, REAL(data));
  
  Eigen::MatrixXd mat = RowMajorVector_to_ColMajorMatrix(REAL(data), countx, county);
  
  return( mat );
  
}




// Read block from hdf5 matrix
Eigen::MatrixXd GetCurrentBlock_hdf5_Original( H5File* file, DataSet* dataset,
                                      hsize_t offsetx, hsize_t offsety, 
                                      hsize_t countx, hsize_t county)
{
  
  IntegerVector offset = IntegerVector::create(offsetx, offsety) ;
  IntegerVector count = IntegerVector::create(countx, county) ;
  IntegerVector stride = IntegerVector::create(1, 1) ;
  IntegerVector block = IntegerVector::create(1, 1) ;
  
  NumericMatrix data(countx, county);
  
  // read_HDF5_matrix_subset(filename, dataset, offset, count, stride, block, REAL(data));
  read_HDF5_matrix_subset(file, dataset, offset, count, stride, block, REAL(data));
  
  Eigen::MatrixXd mat = as<Eigen::MatrixXd>(data);
  
  return( mat );
  
}