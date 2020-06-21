#include "include/rhdf5Utils.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;


// Check if file exists
bool ResFileExist(const std::string& name) {
  
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 

}


// Remove file if exists
bool RemoveFile(std::string filename) 
{
  try {
    
    if( remove( filename.c_str() ) != 0 )
      return 0;
    return 0;
    
  }catch(FileIException error){
    error.printErrorStack();
    return -1;
  }
   
}



extern "C" {
  
  
  
  
  // Create or open a file
  H5FilePtr Open_hdf5_file(const std::string& fname)
  {
    H5::Exception::dontPrint();
    H5::H5File* file = 0;
    
    try {
      file = new H5::H5File(fname.c_str(), H5F_ACC_RDWR);
    } catch(const H5::FileIException&) {
      file = new H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
    }
    
    return H5FilePtr(file);
  }
  
  

  /* Create a group in hdf5 file */
  int create_HDF5_group(H5std_string filename, const H5std_string mGroup)
  {
    try
    {
      Exception::dontPrint();
      
      // Open file and write group
      H5File file(filename, H5F_ACC_RDWR);
      Group group(file.createGroup(mGroup));
      
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
  
  
  
  /* Create empty dataset in hdf5 file */
  int create_HDF5_dataset(H5std_string filename, const std::string CDatasetName, 
                        const size_t rows, const size_t cols, std::string strdatatype)
  {
    
    try
    {
      Exception::dontPrint();
      
      hsize_t     dimsf[2];              // dataset dimensions
      dimsf[0] = rows;
      dimsf[1] = cols;
      
      H5File file(filename, H5F_ACC_RDWR);
      
      DataSpace dataspace( RANK2, dimsf );
      
      if( strdatatype == "int") {
        IntType datatype( PredType::NATIVE_INT );
        DataSet dataset = file.createDataSet( CDatasetName, datatype, dataspace );
        dataset.close();
      } else {
        IntType datatype( PredType::NATIVE_DOUBLE ); 
        DataSet dataset = file.createDataSet( CDatasetName, datatype, dataspace );
        dataset.close();
      }
      
      dataspace.close();
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
  
  

  int write_HDF5_matrix(H5std_string filename, const std::string CDatasetName, RObject DatasetValues)
  {
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      // Open file
      H5File file(filename, H5F_ACC_RDWR);
      
      // Create the data space for the dataset.
      std::vector<int> dims;
      if(is<NumericMatrix>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        dims[0] = as<NumericMatrix>(DatasetValues).rows();
        dims[1] = as<NumericMatrix>(DatasetValues).cols();
        DataSpace dataspace(RANK2, dims);
        
        std::vector<double> matHiCValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));
        
        DataSet dataset = file.createDataSet(CDatasetName,PredType::NATIVE_DOUBLE, dataspace);
        dataset = file.openDataSet(CDatasetName);

        dataset.write( &matHiCValues[0] , PredType::NATIVE_DOUBLE);
        
        dataset.close();
        dataspace.close();
          

      } 
      else if( Rcpp::is<IntegerMatrix>(DatasetValues)) 
      {
        hsize_t dims[2];
        dims[0] = as<IntegerMatrix>(DatasetValues).rows();
        dims[1] = as<IntegerMatrix>(DatasetValues).cols();
        DataSpace dataspace(RANK2, dims);
        
        std::vector<double> matHiCValues = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));
  
        DataSet dataset = file.createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
        dataset = file.openDataSet(CDatasetName);
        dataset.write( &matHiCValues[0], PredType::NATIVE_DOUBLE);
        
        dataset.close();
        dataspace.close();
      } 
      
      file.close();
    
    } 
    catch(FileIException error) { // catch failure caused by the H5File operations
      error.printErrorStack();
      return -1;
    } catch(DataSetIException error) { // catch failure caused by the DataSet operations
      error.printErrorStack();
      return -1;
    } catch(DataSpaceIException error) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    }
    
    return 0;
  }
  

  
  
  
  int read_HDF5_matrix_subset(H5File* file, DataSet* dataset, 
                              IntegerVector ivoffset, IntegerVector ivcount,
                              IntegerVector ivstride, IntegerVector ivblock,
                              double* rdatablock);
  
  

  int write_HDF5_matrix_subset_v2( H5File* file, DataSet* dataset,
                                   IntegerVector ivoffset, IntegerVector ivcount,
                                   IntegerVector ivstride, IntegerVector ivblock,
                                   RObject DatasetValues)
  {
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      // Open file
      //..// H5File file(filename, H5F_ACC_RDWR);
      //..// DataSet dataset = file.openDataSet(CDatasetName);
      
      hsize_t offset[2], count[2], stride[2], block[2];
      
      // Specify size and shape of subset to write
      offset[0] = ivoffset(0); offset[1] = ivoffset(1);
      count[0]  = ivcount(0); count[1]  = ivcount(1);
      stride[0] = ivstride(0); stride[1] = ivstride(1); // default 1
      block[0] = ivblock(0); block[1] = ivblock(1); // default 1
      
      
      // Create the data space for the dataset.
      std::vector<int> dims;
      if(is<NumericMatrix>(DatasetValues) || Rcpp::is<IntegerMatrix>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        dims[0] = ivcount[0];
        dims[1] = ivcount[1];
        DataSpace dataspace(RANK2, dims);
        
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset->getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
        
        // Write a subset of data to the dataset
        std::vector<double> matdata = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));
        
        dataset->write(&matdata[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();
        
      } 
      else if(is<NumericVector>(DatasetValues) || is<IntegerVector>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        dims[0] = ivcount[0];
        dims[1] = ivcount[1];
        DataSpace dataspace(RANK2, dims);
        
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset->getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
        
        // Write a subset of data to the dataset
        std::vector<double> matdata = as<std::vector<double> >(as<NumericVector>(DatasetValues));
        
        dataset->write(&matdata[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        dataspace.close();  
      } 
      
      
      

    } 
    catch(FileIException error) { // catch failure caused by the H5File operations
      error.printErrorStack();
      return -1;
    } catch(DataSetIException error) { // catch failure caused by the DataSet operations
      error.printErrorStack();
      return -1;
    } catch(DataSpaceIException error) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    }
    
    return 0;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
    
  

  int write_HDF5_matrix_subset(H5std_string filename, const std::string CDatasetName,
                               IntegerVector ivoffset, IntegerVector ivcount,
                               IntegerVector ivstride, IntegerVector ivblock,
                               RObject DatasetValues)
  {
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      // Open file
      H5File file(filename, H5F_ACC_RDWR);
      DataSet dataset = file.openDataSet(CDatasetName);
      
      hsize_t offset[2], count[2], stride[2], block[2];
      
      // Specify size and shape of subset to write
      offset[0] = ivoffset(0); offset[1] = ivoffset(1);
      count[0]  = ivcount(0); count[1]  = ivcount(1);
      stride[0] = ivstride(0); stride[1] = ivstride(1); // default 1
      block[0] = ivblock(0); block[1] = ivblock(1); // default 1
      

      // Create the data space for the dataset.
      std::vector<int> dims;
      if(is<NumericMatrix>(DatasetValues) || Rcpp::is<IntegerMatrix>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        dims[0] = ivcount[0];
        dims[1] = ivcount[1];
        DataSpace dataspace(RANK2, dims);
        
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
        
        // Write a subset of data to the dataset
        std::vector<double> matdata = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));
        
        dataset.write(&matdata[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        
        
      } 
      else if(is<NumericVector>(DatasetValues) || is<IntegerVector>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        dims[0] = ivcount[0];
        dims[1] = ivcount[1];
        DataSpace dataspace(RANK2, dims);
        
        DataSpace memspace(RANK2, dims, NULL);
        
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 

        // Write a subset of data to the dataset
        std::vector<double> matdata = as<std::vector<double> >(as<NumericVector>(DatasetValues));
        
        dataset.write(&matdata[0], PredType::NATIVE_DOUBLE, memspace, dataspace);
        
        dataset.close();
        dataspace.close();
        
      } 
      
      file.close();
      
    } 
    catch(FileIException error) { // catch failure caused by the H5File operations
      error.printErrorStack();
      return -1;
    } catch(DataSetIException error) { // catch failure caused by the DataSet operations
      error.printErrorStack();
      return -1;
    } catch(DataSpaceIException error) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    }
    
    return 0;
  }
  
  
  
  // Read rhdf5 data matrix subset, readed data in rdatablock
  int read_HDF5_matrix_subset(H5File* file, DataSet* dataset, 
                                  IntegerVector ivoffset, IntegerVector ivcount,
                                  IntegerVector ivstride, IntegerVector ivblock,
                                  double* rdatablock)
  {
    
    try
    {
      
      hsize_t offset[2], count[2], stride[2], block[2];
      offset[0] = ivoffset[0]; offset[1] = ivoffset[1];
      count[0] = ivcount[0]; count[1] = ivcount[1];
      stride[0] = ivstride[0]; stride[1] = ivstride[1];
      block[0] = ivblock[0]; block[1] = ivblock[1];

      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      // Define Memory Dataspace. Get file dataspace and select a subset from the file dataspace.
      hsize_t dimsm[2];
      dimsm[0] = count[0]; 
      dimsm[1] = count[1];

      DataSpace memspace(RANK2, dimsm, NULL);
      
      //  Get dataspace of the dataset.
      DataSpace dataspace = dataset->getSpace();
      dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
      
      H5T_class_t type_class = dataset->getTypeClass();
      
      // Get class of datatype and print message if it's an integer.
      if( type_class == H5T_INTEGER )
        dataset->read( rdatablock, PredType::NATIVE_INT, memspace, dataspace );
      else if (type_class == H5T_FLOAT)
        dataset->read( rdatablock, PredType::NATIVE_DOUBLE, memspace, dataspace );
      
      dataspace.close();
    }  // end of try block
    
    catch( FileIException error ){ // catch failure caused by the H5File operations
      error.printErrorStack();
      return -1;
    } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
      error.printErrorStack();
      return -1;
    } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    }
    
    return 0;  // successfully terminated
    
  }
  
  
  
  
  
  
  
  
  
  
 /***  A ELIMINAR ....

  // Read rhdf5 data matrix subset, readed data in rdatablock
  int read_HDF5_matrix_subset (H5std_string filename, const std::string CDatasetName,
                         IntegerVector ivoffset, IntegerVector ivcount,
                         IntegerVector ivstride, IntegerVector ivblock,
                         double* rdatablock)
  {
    
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      hsize_t offset[2], count[2], stride[2], block[2];
      hsize_t dimsm[2];
      
      // Specify size and shape of subset to read
      offset[0] = ivoffset(0); offset[1] = ivoffset(1);
      count[0]  = ivcount(0); count[1]  = ivcount(1);
      stride[0] = ivstride(0); stride[1] = ivstride(1); // default 1
      block[0] = ivblock(0); block[1] = ivblock(1); // default 1
      
      // Define Memory Dataspace. Get file dataspace and select a subset from the file dataspace.
      dimsm[0] = ivcount(0); dimsm[1] = ivcount(1);

      Exception::dontPrint();
      
      //  Open the specified file and the specified dataset in the file.
      H5File file( filename, H5F_ACC_RDONLY );
      DataSet dataset = file.openDataSet( CDatasetName );
      
      DataSpace memspace(RANK2, dimsm, NULL);
      
      //  Get dataspace of the dataset.
      DataSpace dataspace = dataset.getSpace();
      dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
      
      H5T_class_t type_class = dataset.getTypeClass();
      
      // Get class of datatype and print message if it's an integer.
      if( type_class == H5T_INTEGER )
        dataset.read( rdatablock, PredType::NATIVE_INT, memspace, dataspace );
      else if (type_class == H5T_FLOAT)
        dataset.read( rdatablock, PredType::NATIVE_DOUBLE, memspace, dataspace );
  
      dataset.close();
      dataspace.close();
      file.close();
    }  // end of try block
    
    catch( FileIException error ){ // catch failure caused by the H5File operations
      error.printErrorStack();
      return -1;
    } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
      error.printErrorStack();
      return -1;
    } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    }
    
    return 0;  // successfully terminated
    
  }

  ***/
 
 
  
  /*** FER LA VERSIÓ PARAL·LELA

  // Read rhdf5 data matrix subset, readed data in rdatablock
  int read_HDF5_matrix_subset_parallel (H5std_string filename, const std::string CDatasetName,
                                        IntegerVector ivoffset, IntegerVector ivcount,
                                        IntegerVector ivstride, IntegerVector ivblock,
                                        double* rdatablock)
  {
    
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      const int RANK = 2; // 2D-matrix
      
      hsize_t offset[2], count[2], stride[2], block[2];
      hsize_t dimsm[2];
      
      // Specify size and shape of subset to read
      offset[0] = ivoffset(0); offset[1] = ivoffset(1);
      count[0]  = ivcount(0); count[1]  = ivcount(1);
      stride[0] = ivstride(0); stride[1] = ivstride(1); // default 1
      block[0] = ivblock(0); block[1] = ivblock(1); // default 1
      
      // Define Memory Dataspace. Get file dataspace and select a subset from the file dataspace.
      dimsm[0] = ivcount(0); dimsm[1] = ivcount(1);
      
      Exception::dontPrint();
      
      //  Open the specified file and the specified dataset in the file.
      H5File file( filename, H5F_ACC_RDONLY | H5F_ACC_SWMR_READ, H5P_DEFAULT );
      DataSet dataset = file.openDataSet( CDatasetName );
      
      DataSpace memspace(RANK, dimsm, NULL);
      
      //  Get dataspace of the dataset.
      DataSpace dataspace = dataset.getSpace();
      dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
      
      H5T_class_t type_class = dataset.getTypeClass();
      
      // Get class of datatype and print message if it's an integer.
      if( type_class == H5T_INTEGER )
        dataset.read( rdatablock, PredType::NATIVE_INT, memspace, dataspace );
      else if (type_class == H5T_FLOAT)
        dataset.read( rdatablock, PredType::NATIVE_DOUBLE, memspace, dataspace );
      
      dataset.close();
      dataspace.close();
      file.close();
    }  // end of try block
    
    catch( FileIException error ){ // catch failure caused by the H5File operations
      error.printErrorStack();
      return -1;
    } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
      error.printErrorStack();
      return -1;
    } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
      error.printErrorStack();
      return -1;
    }
    
    return 0;  // successfully terminated
    
  }
  ***/
  
  
}


// Create hdf5 data file with matrix data
int Create_hdf5_file(std::string filename)
{
  
  
  try { 
    
    // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
    Exception::dontPrint();
    
    H5File file(filename, H5F_ACC_RDWR); // Open HDF5 file
    file.close();
    
  } catch(const H5::FileIException&) {
    
    // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
    Exception::dontPrint();
    
    H5File file(filename, H5F_ACC_TRUNC); // Create HDF5 file
    file.close();
  }
    
  return(0);
}





//' Write data to hdf5 file
//'
//' Creates a hdf5 file with numerical data matrix,
//' 
//' @param filename, character array indicating the name of the file to create
//' @param folder, character array indicating folder name to put the matrix in hdf5 file
//' @param mat numerical data matrix
//' @return none
//' @export
// [[Rcpp::export]]
int Create_HDF5_matrix_file(std::string filename, RObject mat, Rcpp::Nullable<std::string> folder = R_NilValue)
{
  
  try
  {
    int res;
    std::string strsubgroup;
    
    if ( mat.sexp_type()==0   )
      throw std::range_error("Data matrix must exsits and mustn't be null");
    
    Eigen::MatrixXd matdat = as<Eigen::MatrixXd>(mat);
    
    // Create HDF5 file
    H5File file(filename, H5F_ACC_TRUNC);
    
    if( folder.isNull()) {
      strsubgroup = "Base.matrices";
    } else {
      strsubgroup = Rcpp::as<std::string> (folder);
    }

    res = create_HDF5_group(filename, strsubgroup );

    res = write_HDF5_matrix(filename, strsubgroup + "/matrix", as<Rcpp::NumericMatrix>(mat) );
    
    IntegerVector offset = IntegerVector::create(1,2) ;
    IntegerVector count = IntegerVector::create(3,4) ;
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1); 
    
    NumericMatrix data(count(1), count(0));
    
    //s'hauria de canviar la funció, aquesta està obsoleta....// read_HDF5_matrix_subset(filename, strsubgroup + "/Base.matrices", offset, count, stride, block, REAL(data));
    Rcpp::Rcout<<data;

    file.close();
    
  }
  catch( FileIException error ) { // catch failure caused by the H5File operations
    error.printErrorStack();
    return -1;
  }
  
  return(0);
}

/***R

library(BigDataStatMeth)
library(DelayedArray)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/BigDataStatMeth/tmp")

N <- 6
M <- 2
A <- matrix(rnorm(N*M,mean=0,sd=1), N, M)
Ad <- DelayedArray(A)
B <- matrix(sample.int(15, 9*100, TRUE), 20, 20)

A[1:10,1:10]


A <- matrix(c(1,2,3,4,5,6,7,8,9,10),2,5, byrow = TRUE)
Create_HDF5_matrix_file("prova.hdf5" ,A)
# Create_HDF5_matrix_file("prova.hdf5" ,B)



h5f

*/