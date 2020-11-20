#include "include/ReadDelayedData.h"


// [[Rcpp::depends(RcppEigen)]]


Eigen::MatrixXd read_DelayedArray_int( Rcpp::RObject A )
{

  auto dmat = beachmat::create_integer_matrix(A);
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  Eigen::MatrixXd eX(nrows, ncols);
  
  if ( nrows < ncols )
  {
    Rcpp::IntegerVector output(dmat->get_nrow());
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      dmat->get_col(ncol, output.begin());
      eX.col(ncol) = Rcpp::as<Eigen::VectorXd >(output);
    }
  }else {
    Rcpp::IntegerVector output(dmat->get_ncol());
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      eX.row(nrow) = Rcpp::as<Eigen::VectorXd >(output);
    }
  }
  return eX;
}


Eigen::MatrixXd read_DelayedArray_real( Rcpp::RObject A )
{
  // size_t ncols = 0, nrows=0;
  
  auto dmat = beachmat::create_numeric_matrix(A);

  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  Eigen::MatrixXd eX(nrows, ncols);
  // eX.resize(nrows, ncols);
  
  
  if ( ncols < nrows )
  {
    Rcpp::NumericVector output(nrows);
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      dmat->get_col(ncol, output.begin());
      eX.col(ncol) = Rcpp::as<Eigen::VectorXd >(output);
    }
  }else {
    Rcpp::NumericVector output(ncols);
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      eX.row(nrow) = Rcpp::as<Eigen::VectorXd >(output);
    }
  }
  return eX;
}


// Get matrix data in memory from DelayedArray
Eigen::MatrixXd read_DelayedArray( Rcpp::RObject A  )
{
  auto dmtypex = beachmat::find_sexp_type(A);

  if ( dmtypex == INTSXP )
  {
    return (read_DelayedArray_int(A));
  } else if (dmtypex==REALSXP)
  {
    return (read_DelayedArray_real(A));
  }else
  {
    stop("unacceptable matrix type");
  }
  
}


// Get number of rows and columns from delayed array data
Eigen::Vector2i get_DelayedArray_size(Rcpp::RObject A)
{
  auto dmtypex = beachmat::find_sexp_type(A);
  Eigen::Vector2i size;
  
  if ( dmtypex == INTSXP )
  {
    
    auto dmat = beachmat::create_integer_matrix(A);
    size[0] = dmat->get_nrow(); 
    size[1] = dmat->get_ncol();
    
  } else if (dmtypex==REALSXP) {
    
    auto dmat = beachmat::create_numeric_matrix(A);
    size[0] = dmat->get_nrow(); 
    size[1] = dmat->get_ncol();
    
  }else {
    stop("unacceptable matrix type");
  }
  return(size);
}


// Write DelayedArray to hdf5 data file 
int write_DelayedArray_to_hdf5(H5std_string filename, const std::string CDatasetName, Rcpp::RObject A)
{
  
  int res = 0;
  auto dmtypex = beachmat::find_sexp_type(A);
  
  if ( dmtypex == INTSXP )  {
    res = write_DelayedArray_int_hdf5(filename, CDatasetName, A);
    
  } else if (dmtypex==REALSXP)  {
    res = write_DelayedArray_real_hdf5_transposed(filename, CDatasetName, A);
    
  }else  {
    stop("unacceptable matrix type");
  }
  
  return(res);
}


// Write DelayedArray to hdf5 data file 
int write_DelayedArray_to_hdf5_ptr(H5File* file, const std::string CDatasetName, Rcpp::RObject A, bool transp)
{
  
  int res = 0;
  auto dmtypex = beachmat::find_sexp_type(A);
  
  if ( dmtypex == INTSXP )  {
    if ( transp == true)
      res = write_DelayedArray_int_hdf5_ptr(file, CDatasetName, A);
    else
      res = write_DelayedArray_int_hdf5_transposed_ptr(file, CDatasetName, A);
    
  } else if (dmtypex==REALSXP)  {
    if(transp==true)
      res= write_DelayedArray_real_hdf5_ptr(file,CDatasetName, A);
    else
      res = write_DelayedArray_real_hdf5_transposed_ptr(file, CDatasetName, A);
    
  }else  {
    stop("unacceptable matrix type");
  }
  
  return(res);
}




// Write DealyedArray with integer data type to hdf5 file as double
int write_DelayedArray_int_hdf5( H5std_string filename, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_integer_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();
    
    // Create empty dataset in hdf5 file
    create_HDF5_dataset( filename, CDatasetName, nrows, ncols, "real");
    
    // Open file and dataset
    H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file.openDataSet(CDatasetName);
    
    if ( ncols < nrows ) 
    {
      offset[0] = 0;
      dims[0] = count[0] = nrows;
      dims[1] = count[1] = 1;
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[1] = ncol;
        dmat->get_col(ncol, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    else 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = ncols;
      
      Rcpp::NumericVector output(ncols);
      for (size_t nrow=0; nrow<nrows; ++nrow)   // Write by rows
      {
        offset[0] = nrow;
        dmat->get_row(nrow, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    dataset.close();
    file.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}


/***
// Write DealyedArray with integer data type to hdf5 file as double
int write_DelayedArray_int_hdf5_ptr( H5File* file, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_integer_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();

    // Create empty dataset in hdf5 file
    create_HDF5_dataset_ptr( file, CDatasetName, nrows, ncols, "real");
    
    // Open file and dataset
    //..// H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file->openDataSet(CDatasetName);
    
    if ( ncols < nrows ) 
    {
      offset[0] = 0;
      dims[0] = count[0] = nrows;
      dims[1] = count[1] = 1;
      
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[1] = ncol;
        dmat->get_col(ncol, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
      
    }
    else 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = ncols;
      
      Rcpp::NumericVector output(ncols);
      for (size_t nrow=0; nrow<nrows; ++nrow)   // Write by rows
      {
        offset[0] = nrow;
        dmat->get_row(nrow, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    dataset.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}
***/



// Write DealyedArray with integer data type to hdf5 file as double
int write_DelayedArray_int_hdf5_ptr( H5File* file, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_integer_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();
    
    // Create empty dataset in hdf5 file
    create_HDF5_dataset_ptr( file, CDatasetName, ncols, nrows, "real");
    
    // Open file and dataset
    //..// H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file->openDataSet(CDatasetName);
    
    if ( ncols < nrows ) 
    {
      offset[0] = 0;
      dims[0] = count[0] = nrows;
      dims[1] = count[1] = 1;
      
      DataSpace dataspace(RANK2, dims);
      DataSpace memspace(RANK2, dims, NULL);
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[1] = ncol;
        dmat->get_col(ncol, output.begin());
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        
      }
      
      dataspace.close();
      
    }
    else 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = ncols;
      
      DataSpace dataspace(RANK2, dims);
      DataSpace memspace(RANK2, dims, NULL);
      
      Rcpp::NumericVector output(ncols);
      for (size_t nrow=0; nrow<nrows; ++nrow)   // Write by rows
      {
        offset[0] = nrow;
        dmat->get_row(nrow, output.begin());

        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        
      }
      dataspace.close();
    }
    dataset.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}






// Write DealyedArray with integer data type to hdf5 file as double
// writes data as row-major (hdf5 file format), equivalent to write transposed data.
int write_DelayedArray_int_hdf5_transposed( H5std_string filename, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_integer_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();
    
    // Create empty dataset in hdf5 file (transposed data nrows = ncols)
    create_HDF5_dataset( filename, CDatasetName, nrows, ncols,"real");
    
    // Open file and dataset
    H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file.openDataSet(CDatasetName);
    
    Rcpp::Rcout<<"\nInteger\n";
    if ( ncols < nrows ) 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = nrows;
      
      DataSpace dataspace(RANK2, dims);
      DataSpace memspace(RANK2, dims, NULL);
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[0] = ncol;
        dmat->get_col(ncol, output.begin());
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        
      }
      
      dataspace.close();
    }
    else 
    {
      offset[1] = 0;
      dims[0] = count[0] = ncols;
      dims[1] = count[1] = 1;

      DataSpace dataspace(RANK2, dims);
      DataSpace memspace(RANK2, dims, NULL);
      
      Rcpp::NumericVector output(ncols);
      for (size_t nrow=0; nrow<nrows; ++nrow)   // Write by rows
      {
        offset[0] = nrow;
        dmat->get_row(nrow, output.begin());

        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        
      }
      dataspace.close();
    }
    dataset.close();
    file.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}


// Write DealyedArray with integer data type to hdf5 file as double
// writes data as row-major (hdf5 file format), equivalent to write transposed data.
int write_DelayedArray_int_hdf5_transposed_ptr( H5File* file, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_integer_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();
    
    // Create empty dataset in hdf5 file (transposed data nrows = ncols)
    create_HDF5_dataset_ptr( file, CDatasetName, ncols, nrows, "real");
    
    // Open file and dataset
    //..// H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file->openDataSet(CDatasetName);
    
    Rcpp::Rcout<<"\nInteger\n";
    if ( ncols <= nrows ) 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = nrows;
      
      DataSpace dataspace(RANK2, dims);
      DataSpace memspace(RANK2, dims, NULL);
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[0] = ncol;
        dmat->get_col(ncol, output.begin());
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        
      }
      dataspace.close();
    }
    else 
    {
      offset[1] = 0;
      dims[0] = count[0] = ncols;
      dims[1] = count[1] = 1;

      DataSpace dataspace(RANK2, dims);
      DataSpace memspace(RANK2, dims, NULL);
      
            
      Rcpp::NumericVector output(ncols,1);
      for (size_t nrow=0; nrow<nrows; ++nrow)   // Write by rows
      {
        offset[0] = 0;
        offset[1] = nrow;
        dmat->get_row(nrow, output.begin());
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        
      }
      dataspace.close();
    }
    dataset.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}





// Write DealyedArray with double data type to hdf5 file
int write_DelayedArray_real_hdf5( H5std_string filename, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_numeric_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();
    
    // Create empty dataset in hdf5 file
    create_HDF5_dataset( filename, CDatasetName, nrows, ncols,"real");
    
    // Open file and dataset
    H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file.openDataSet(CDatasetName);
    
    if ( ncols < nrows ) 
    {
      offset[0] = 0;
      dims[0] = count[0] = nrows;
      dims[1] = count[1] = 1;
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[1] = ncol;
        dmat->get_col(ncol, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    else 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = ncols;
      
      Rcpp::NumericVector output(ncols);
      // Write by rows
      for (size_t nrow=0; nrow<nrows; ++nrow) 
      {
        offset[0] = nrow;
        dmat->get_row(nrow, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    
    dataset.close();
    file.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}




// Write DealyedArray with double data type to hdf5 file
int write_DelayedArray_real_hdf5_ptr( H5File* file, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_numeric_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();
    
    // Create empty dataset in hdf5 file
    create_HDF5_dataset_ptr( file, CDatasetName, nrows, ncols,"real");
    
    // Open file and dataset
    //..// H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file->openDataSet(CDatasetName);
    
    if ( ncols < nrows ) 
    {
      offset[0] = 0;
      dims[0] = count[0] = nrows;
      dims[1] = count[1] = 1;
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[1] = ncol;
        dmat->get_col(ncol, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    else 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = ncols;
      
      Rcpp::NumericVector output(ncols);
      // Write by rows
      for (size_t nrow=0; nrow<nrows; ++nrow) 
      {
        offset[0] = nrow;
        dmat->get_row(nrow, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    
    dataset.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}






// Write DealyedArray with integer data type to hdf5 file as double
// writes data as row-major (hdf5 file format), equivalent to write transposed data.
int write_DelayedArray_real_hdf5_transposed( H5std_string filename, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_numeric_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();
    
    // Create empty dataset in hdf5 file (transposed data nrows = ncols)
    create_HDF5_dataset( filename, CDatasetName, ncols, nrows, "real");
    
    // Open file and dataset
    H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file.openDataSet(CDatasetName);

    if ( ncols < nrows ) 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = nrows;
      /*dims[0] = 1;
      dims[1] = nrows;
      count[0] = 1;
      count[1] = nrows;*/
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[0] = ncol;
        dmat->get_col(ncol, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    else 
    {
      offset[0] = 0;
      dims[0] = count[0] = ncols;
      dims[1] = count[1] = 1;
      
      Rcpp::NumericVector output(ncols);
      for (size_t nrow=0; nrow<nrows; ++nrow)   // Write by rows
      {
        offset[1] = nrow;
        dmat->get_row(nrow, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    dataset.close();
    file.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}



// Write DealyedArray with integer data type to hdf5 file as double
// writes data as row-major (hdf5 file format), equivalent to write transposed data.
int write_DelayedArray_real_hdf5_transposed_ptr( H5File* file, const std::string CDatasetName, Rcpp::RObject A )
{
  
  hsize_t offset[2], count[2], dims[2];
  hsize_t stride[2] = {1,1};
  hsize_t block[2] = {1,1};
  
  auto dmat = beachmat::create_numeric_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  try
  {
    Exception::dontPrint();
    
    // Create empty dataset in hdf5 file (transposed data nrows = ncols)
    //..// int res = create_HDF5_dataset( filename, CDatasetName, ncols, nrows, "real");
    create_HDF5_dataset_ptr( file, CDatasetName, ncols, nrows, "real");
    
    // Open file and dataset
    //..// H5File file(filename, H5F_ACC_RDWR);
    DataSet dataset = file->openDataSet(CDatasetName);
    
    if ( ncols < nrows ) 
    {
      offset[1] = 0;
      dims[0] = count[0] = 1;
      dims[1] = count[1] = nrows;
      /*dims[0] = 1;
       dims[1] = nrows;
       count[0] = 1;
       count[1] = nrows;*/
      
      Rcpp::NumericVector output(nrows);
      for (size_t ncol=0; ncol<ncols; ++ncol) 
      {
        offset[0] = ncol;
        dmat->get_col(ncol, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    else 
    {
      offset[0] = 0;
      dims[0] = count[0] = ncols;
      dims[1] = count[1] = 1;
      
      Rcpp::NumericVector output(ncols);
      for (size_t nrow=0; nrow<nrows; ++nrow)   // Write by rows
      {
        offset[1] = nrow;
        dmat->get_row(nrow, output.begin());
        
        DataSpace dataspace(RANK2, dims);
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(&output[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
      }
    }
    dataset.close();
    
  } // end of try block
  catch(FileIException error) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    error.printErrorStack();
    return -1;
  }
  
  return 0;
}




/*** ELIMINAR  QUAN TOT OK

int write_DelayedArray_real_hdf5_v2( H5std_string filename, const std::string CDatasetName, Rcpp::RObject A )
{
  // size_t ncols = 0, nrows=0;
  
  Rcpp::IntegerVector ivoffset(2);
  Rcpp::IntegerVector ivcount(2);
  Rcpp::IntegerVector ivstride = {1,1};
  Rcpp::IntegerVector ivblock = {1,1};
  
  auto dmat = beachmat::create_numeric_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  int res = create_HDF5_dataset( filename, CDatasetName, nrows, ncols,"real");
  
  if ( ncols < nrows )
  {
    ivoffset[0] = 0;
    ivcount[0] = nrows;
    ivcount[1] = 1;
    
    Rcpp::NumericVector output(nrows);
    // Write by columns
    Rcpp::Rcout<<"Write by columns..";
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      Rcpp::Rcout<<"Writign column ..."<<ncol<<"\n";
      ivoffset[1] = ncol;
      dmat->get_col(ncol, output.begin());
      write_HDF5_matrix_subset( filename, CDatasetName, ivoffset, ivcount,
                               ivstride, ivblock,  Rcpp::as<Rcpp::NumericVector >(output) );
    }
  }else {
    
    ivoffset[1] = 0;
    ivcount[0] = 1;
    ivcount[1] = ncols;
    
    Rcpp::NumericVector output(ncols);
    // Write by rows
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      ivoffset[0] = nrow;
      
      write_HDF5_matrix_subset( filename, CDatasetName, ivoffset, ivcount,
                                ivstride, ivblock,  Rcpp::as<Rcpp::NumericVector >(output) );
    }
  }
  return 0;
}

***/



Rcpp::NumericMatrix read_DelayedArray_int_r( Rcpp::RObject A )
{

  auto dmat = beachmat::create_integer_matrix(A);
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  Rcpp::NumericMatrix eX(nrows, ncols);

  if ( nrows < ncols )
  {
    Rcpp::IntegerVector output(dmat->get_nrow());
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      dmat->get_col(ncol, output.begin());
      eX.column(ncol) = output;
    }
  }else {
    Rcpp::IntegerVector output(dmat->get_ncol());
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      eX.row(nrow) = output;
    }
  }
  return eX;
}


Rcpp::NumericMatrix read_DelayedArray_real_r( Rcpp::RObject A )
{
  auto dmat = beachmat::create_numeric_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  Rcpp::NumericMatrix eX(nrows, ncols);
  
  if ( ncols < nrows )
  {
    Rcpp::NumericVector output(nrows);
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      dmat->get_col(ncol, output.begin());
      eX.column(ncol) = output;
    }
  }else {
    Rcpp::NumericVector output(ncols);
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      eX.row(nrow) = output;
    }
  }
  return eX;
}




Rcpp::NumericMatrix read_DelayedArray_rcpp( Rcpp::RObject A  )
{
  auto dmtypex = beachmat::find_sexp_type(A);
  // size_t ncols = 0, nrows=0;
  
  if ( dmtypex == INTSXP )
  {
    return (read_DelayedArray_int_r(A));
  } else if (dmtypex==REALSXP)
  {
    return (read_DelayedArray_real_r(A));
  }else
  {
    stop("unacceptable matrix type");
  }
  
}



