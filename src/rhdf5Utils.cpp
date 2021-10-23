#include "include/rhdf5Utils.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;


// Check if file exists
bool ResFileExist(const std::string& name) {
  
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
  
}


// // Check if file exists
bool ResFileExist_filestream(std::string name) {

  bool exists = true;

  fstream fileStream;
  fileStream.open(name);

  if (fileStream.fail()) {
    exists = false;
  }

  return(exists);

}






// Remove file if exists
bool RemoveFile(std::string filename) 
{
  try {
    
    if( remove( filename.c_str() ) != 0 )
      return 0;
    return 0;
    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception RemoveFile (File IException)" );
    return -1;
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception RemoveFile (DataSet IException)" );
    return -1;
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception RemoveFile (Group IException)" );
    return -1;
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception RemoveFile (DataSpace IException)" );
    return -1;
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception RemoveFile (Data TypeIException)" );
    return -1;
  }
   
}



extern "C" {
  
  
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
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception create_HDF5_group (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception create_HDF5_group (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception create_HDF5_group (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_group (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_group (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  
  
  bool pathExists(hid_t id, const std::string& path)
  {
    return H5Lexists( id, path.c_str(), H5P_DEFAULT ) > 0;
  }
  
  
  
  
  // Join multiple datasets in one dataset in the same group
  int join_datasets(H5File* file, std::string strsubgroup, StringVector strinput, std::string strasout)
  {
    try{
      
      Exception::dontPrint();
      
      IntegerVector stride = IntegerVector::create(1, 1);
      IntegerVector block = IntegerVector::create(1, 1);
      IntegerVector offset = IntegerVector::create(0, 0);
      IntegerVector count = IntegerVector::create(0, 0);
      string stroutdataset = strsubgroup + "/" + strasout;
      
      DataSet dataset = file->openDataSet(strsubgroup + "/" + strinput[0]);
      IntegerVector dims_out = get_HDF5_dataset_size(dataset);

      // Create unlimited dataset and add first dataset
      if( exists_HDF5_element_ptr(file,stroutdataset))
        remove_HDF5_element_ptr(file,stroutdataset);
      
      create_HDF5_unlimited_matrix_dataset_ptr( file, stroutdataset, (unsigned long long)dims_out[0], (unsigned long long)dims_out[1], "real");
      NumericMatrix readeddata(dims_out[0], dims_out[1]);
      // read first dataset
      read_HDF5_matrix_subset( file, &dataset, offset, dims_out, stride, block, REAL(readeddata) );
      dataset.close();  
      
      DataSet* unlimDataset = new DataSet(file->openDataSet(stroutdataset));
      
      // write first dataset to new dataset
      write_HDF5_matrix_subset_v2(file, unlimDataset, offset, dims_out, stride, block, transpose(readeddata));  

      // Update offset to new position
      offset[1] = offset[1] + dims_out[1];
      
      for( int i=1; i<strinput.size(); i++)
      {
        dataset = file->openDataSet(strsubgroup + "/" + strinput[i]);
        dims_out = get_HDF5_dataset_size(dataset);
        
        // Update columns size
        count = dims_out;
        
        // Extend dataset before put data
        extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, dims_out[1]);

        // Read data
        NumericMatrix newdata(dims_out[0], dims_out[1]);
        IntegerVector offsetRead = IntegerVector::create(0, 0);

        // read dataset and close it.
        read_HDF5_matrix_subset( file, &dataset, offsetRead, dims_out, stride, block, REAL(newdata) );
        dataset.close();
        
        write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, transpose(newdata)  );  

        // Update offset
        offset[1] = offset[1] + dims_out[1];
      }
      
      unlimDataset->close();

    
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception join_datasets (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception join_datasets (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception join_datasets (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception join_datasets (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception join_datasets (Data TypeIException)" );
      return -1;
    }
    return(0);
  }
  
    
  /* Create a group in hdf5 file */
  int create_HDF5_group_ptr( H5File* file, const H5std_string mGroup)
  {
    try
    {
      Exception::dontPrint();
      
      std::string strgroup = mGroup;
      
      // if group no exists -> Create group to file
      if(!pathExists( file->getId(), strgroup)) 
        file->createGroup("/"+ mGroup);

    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_group_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  
  bool exists_HDF5_element_ptr(H5File* file, const H5std_string element)
  {
    bool bexists = false;
    try
    {
      Exception::dontPrint();

      // Search dataset
      if(pathExists( file->getId(), element)) 
        bexists = true;
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception exists_HDF5_element_ptr (Data TypeIException)" );
      return -1;
    }
    return bexists;
  }
  
  
  // Create multiple group in hdf5 data file, groups must be separated by "/"
  int create_HDF5_groups_ptr( H5File* file, const H5std_string mGroup)
  {
    try
    {
      Exception::dontPrint();
      
      std::string strgroup = mGroup;
      std::string results = "";
      vector<string> result; 
      
      boost::split(result, mGroup, boost::is_any_of("/")); 
      
      for (int i = 0; i < result.size(); i++) {
        if(!pathExists( file->getId(), results + result[i])) 
          file->createGroup(results + result[i]);
        
        results = result[i] + "/";
      }
        
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_groups_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  
  
  
  
  bool remove_HDF5_element_ptr(H5File* file, const H5std_string element)
  {
    
    bool bremok = true;
    
    try
    {
      Exception::dontPrint();

      int result = H5Ldelete(file->getId(), element.data(), H5P_DEFAULT);
      if(result<0)
        bremok = false;
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception remove_HDF5_element_ptr (Data TypeIException)" );
      return -1;
    }
    
    return(bremok);
    
  }
  
  
  bool remove_HDF5_multiple_elements_ptr(H5File* file, std::string strgroup, StringVector elements)
  {
    
    bool bremok = true;
    
    try
    {
      Exception::dontPrint();
      
      for (int i=0; i<elements.size(); i++)
      {
        H5std_string element = strgroup + "/" + elements[i];
        
        int result = H5Ldelete(file->getId(), element.data(), H5P_DEFAULT);  
        if(result<0)
          bremok = false;
        // else
        //  Rcpp::Rcout<<"\n Element : "<<element<<" removed\n";
      }
      
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (Data TypeIException)" );
      return -1;
    }
    
    return(bremok);
    
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
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception create_HDF5_dataset (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception create_HDF5_dataset (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception create_HDF5_dataset (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_dataset (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_dataset (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  
  
  
  /* Create empty dataset in hdf5 file (pointer to file) */
  int create_HDF5_dataset_ptr(H5File* file, const std::string CDatasetName, 
                          const size_t rows, const size_t cols, std::string strdatatype)
  {
    
    try
    {
      Exception::dontPrint();
      
      hsize_t     dimsf[2];              // dataset dimensions
      dimsf[0] = rows;
      dimsf[1] = cols;
      
      DataSpace dataspace( RANK2, dimsf );
      
      if( strdatatype == "int") {
        IntType datatype( PredType::NATIVE_INT );
        DataSet dataset = file->createDataSet( CDatasetName, datatype, dataspace );
        dataset.close();
      } else {
        IntType datatype( PredType::NATIVE_DOUBLE ); 
        DataSet dataset = file->createDataSet( CDatasetName, datatype, dataspace );
        dataset.close();
      }
      
      dataspace.close();

    } catch(FileIException& error) { // catch failure caused by the H5File operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_dataset_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  
  /* Create empty matrix dataset in hdf5 file */
  int create_HDF5_unlimited_matrix_dataset_ptr(H5File* file, const std::string CDatasetName, 
                          const size_t rows, const size_t cols, std::string strdatatype)
  {
    
    try
    {
      Exception::dontPrint();
      
      herr_t status;
      
      hsize_t     dimsf[2];              // dataset dimensions
      dimsf[0] = rows;
      dimsf[1] = cols;
      hid_t cparms; 

      // Declare unlimited dimensions
      hsize_t  maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
      DataSpace dataspace ( RANK2, dimsf, maxdims );
      
      // Enabling chunking
      hsize_t chunk_dims[2];
      chunk_dims[0] = rows;
      chunk_dims[1] = cols;
      
      cparms = H5Pcreate(H5P_DATASET_CREATE);
      status = H5Pset_chunk( cparms, RANK2, chunk_dims);
      
      // Create dataset
      if( strdatatype == "int") {
        IntType datatype( PredType::NATIVE_INT );
        DataSet dataset = file->createDataSet( CDatasetName, datatype, dataspace, cparms);
        dataset.close();
      } else {
        IntType datatype( PredType::NATIVE_DOUBLE ); 
        DataSet dataset = file->createDataSet( CDatasetName, datatype, dataspace, cparms);
        dataset.close();
      }

      dataspace.close();

    } catch(FileIException& error) { // catch failure caused by the H5File operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_unlimited_matrix_dataset_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_unlimited_matrix_dataset_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_unlimited_matrix_dataset_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_unlimited_matrix_dataset_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      file->close();
      ::Rf_error( "c++ exception create_HDF5_unlimited_matrix_dataset_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  
  
  
  /* Create empty vector dataset in hdf5 file */
  int create_HDF5_unlimited_vector_dataset_ptr(H5File* file, const std::string CDatasetName, 
                                               const size_t length, std::string strdatatype)
  {
    
    try
    {
      Exception::dontPrint();
      
      herr_t status;
      
      hsize_t     dimsf[1];              // dataset dimensions
      dimsf[0] = length;
      hid_t cparms; 
      
      // Declare unlimited dimensions
      hsize_t  maxdims[1] = {H5S_UNLIMITED};
      DataSpace dataspace ( RANK1, dimsf, maxdims );
      
      // Enabling chunking
      hsize_t chunk_dims[1];
      chunk_dims[0] = length;
      
      cparms = H5Pcreate(H5P_DATASET_CREATE);
      status = H5Pset_chunk( cparms, RANK1, chunk_dims);
      
      // Create dataset
      if( strdatatype == "int") {
        IntType datatype( PredType::NATIVE_INT );
        DataSet dataset = file->createDataSet( CDatasetName, datatype, dataspace, cparms);
        dataset.close();
      } else if( strdatatype == "char" || strdatatype == "character") {
        
        
        
        // define name struct
        typedef struct name {
          char chr[MAXSTRING];
        } name;
        
        DataSpace dataspace(RANK1, chunk_dims);
        
        cparms = H5Pcreate(H5P_DATASET_CREATE);
        status = H5Pset_chunk( cparms, RANK1, chunk_dims);
        
        // Create the memory datatype.
        H5::CompType mtype(sizeof(name));
        mtype.insertMember("chr", HOFFSET(name, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
        
        DataSet dataset = file->createDataSet( CDatasetName, mtype, dataspace, cparms);
        dataset.close();
        
      } else {
        IntType datatype( PredType::NATIVE_DOUBLE ); 
        DataSet dataset = file->createDataSet( CDatasetName, datatype, dataspace, cparms);
        dataset.close();
      }
      
      dataspace.close();
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception create_HDF5_unlimited_vector_dataset_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception create_HDF5_unlimited_vector_dataset_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception create_HDF5_unlimited_vector_dataset_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_unlimited_vector_dataset_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception create_HDF5_unlimited_vector_dataset_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  

  int extend_HDF5_matrix_subset_ptr(H5File* file, DataSet* dataset, const size_t rows, const size_t cols)
  {
    try
    {
      Exception::dontPrint();
      
      int rank, ndims;
      
      // Get dataspace from dataset
      DataSpace dataspace = dataset->getSpace();
      
      // Get the number of dimensions in the dataspace.
      rank = dataspace.getSimpleExtentNdims();
      
      // Get the dimension size of each dimension in the dataspace and
      hsize_t dims_out[2];
      ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
      
      // Create new dataset size from new dims and old dims
      hsize_t   newdims[2];
      newdims[0] = rows;
      newdims[1] = cols;
      
      hsize_t      size[2];
      size[0]   = dims_out[0] + newdims[0];
      size[1]   = dims_out[1] + newdims[1];

      dataset->extend( size );
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (Data TypeIException)" );
      return -1;
    }
    return 0;
    
  }
  
  
  
  
  
  // This function does'nt work --> Problems wiht !!! =>  dataset->extend( size ); ==> What happens??? 
  int extend_HDF5_vector_subset_ptr(H5File* file, DataSet* dataset, const size_t length)
  {
    try
    {
      Exception::dontPrint();
      
      int rank, ndims;
      
      // Get dataspace from dataset
      DataSpace dataspace = dataset->getSpace();
      
      Rcpp::Rcout<<"\n\tExtenem - 2 \n";
      
      // Get the number of dimensions in the dataspace.
      rank = dataspace.getSimpleExtentNdims();
      
      // Get the dimension size of each dimension in the dataspace and
      hsize_t dims_out[1];
      ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
      
      // Create new dataset size from new dims and old dims
      hsize_t   newdims[1];
      newdims[0] = length;
      hsize_t      size[1];
      size[0]   = dims_out[0] + newdims[0];

      dataset->extend( size );

    } catch(FileIException& error) { // catch failure caused by the H5File operations
      dataset->close();
      file->close();
      ::Rf_error( "c++ exception extend_HDF5_vector_subset_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      dataset->close();
      file->close();
      ::Rf_error( "c++ exception extend_HDF5_vector_subset_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      dataset->close();
      file->close();
      ::Rf_error( "c++ exception extend_HDF5_vector_subset_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      dataset->close();
      file->close();
      ::Rf_error( "c++ exception extend_HDF5_vector_subset_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      dataset->close();
      file->close();
      ::Rf_error( "c++ exception extend_HDF5_vector_subset_ptr (Data TypeIException)" );
      return -1;
    }
    return 0;
    
  }
  
/*
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
    
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_HDF5_matrix (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
*/  
  
  int write_HDF5_matrix_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues)
  {
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();

      // Create the data space for the dataset.
      std::vector<int> dims;
      
      if(is<NumericMatrix>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        dims[0] = as<NumericMatrix>(DatasetValues).rows();
        dims[1] = as<NumericMatrix>(DatasetValues).cols();
        DataSpace dataspace(RANK2, dims);
        
        std::vector<double> matHiCValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));
        
        DataSet dataset = file->createDataSet(CDatasetName,PredType::NATIVE_DOUBLE, dataspace);
        dataset = file->openDataSet(CDatasetName);
        
        dataset.write( &matHiCValues[0] , PredType::NATIVE_DOUBLE);
        
        dataset.close();
        dataspace.close();

      } 
      else if( Rcpp::is<IntegerMatrix>(DatasetValues) ) 
      {
        hsize_t dims[2];
        dims[0] = as<IntegerMatrix>(DatasetValues).rows();
        dims[1] = as<IntegerMatrix>(DatasetValues).cols();
        DataSpace dataspace(RANK2, dims);
        
        std::vector<double> matHiCValues = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));
        
        DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
        dataset = file->openDataSet(CDatasetName);
        dataset.write( &matHiCValues[0], PredType::NATIVE_DOUBLE);
        
        dataset.close();
        dataspace.close();
      } else if(is<NumericVector>(DatasetValues) || is<IntegerVector>(DatasetValues)) 
      {
        
        hsize_t vectorsize;
        
        if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues)) 
          vectorsize = as<IntegerVector>(DatasetValues).length();
        else if (is<NumericVector>(DatasetValues))
          vectorsize = as<NumericVector>(DatasetValues).length();
        else 
          vectorsize = 1;
        
        hsize_t dims[] = {vectorsize};
        DataSpace dataspace(RANK1, dims);
        
        
        if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues) ) 
        {
          std::vector<int> vectHiCValues(dims[0]);
          //..// int vectHiCValues[dims[0]];
          for(int i=0;i<dims[0]; i++)
            vectHiCValues[i] = as<IntegerVector>(DatasetValues)(i);
          
          DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_INT, dataspace);
          dataset = file->openDataSet(CDatasetName);
          dataset.write( vectHiCValues.data(), PredType::NATIVE_INT);
          dataspace.close();
          dataset.close();
        } 
        else if(is<NumericVector>(DatasetValues) ) 
        {
          std::vector<double> vectValues(dims[0]);
          //..// double vectValues[dims[0]];
          for(int i=0;i<dims[0]; i++)
            vectValues[i] = as<NumericVector>(DatasetValues)(i);
          
          DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
          dataset = file->openDataSet(CDatasetName);
          //..// dataset.write(vectValues, PredType::NATIVE_DOUBLE);
          dataset.write( vectValues.data() , PredType::NATIVE_DOUBLE);
          dataspace.close();
          dataset.close();
        } 

      } 
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }

  
  int write_HDF5_matrix_transposed_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues)
  {
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      // Create the data space for the dataset.
      std::vector<int> dims;
      
      if(is<NumericMatrix>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        dims[1] = as<NumericMatrix>(DatasetValues).rows();
        dims[0] = as<NumericMatrix>(DatasetValues).cols();
        DataSpace dataspace(RANK2, dims);
        
        std::vector<double> matValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));
        
        DataSet dataset = file->createDataSet(CDatasetName,PredType::NATIVE_DOUBLE, dataspace);
        dataset = file->openDataSet(CDatasetName);
        
        dataset.write( &matValues[0] , PredType::NATIVE_DOUBLE);
        
        dataset.close();
        dataspace.close();
        
      } 
      else if( Rcpp::is<IntegerMatrix>(DatasetValues) ) 
      {
        hsize_t dims[2];
        dims[1] = as<IntegerMatrix>(DatasetValues).rows();
        dims[0] = as<IntegerMatrix>(DatasetValues).cols();
        DataSpace dataspace(RANK2, dims);
        
        std::vector<double> matValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));
        
        DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
        dataset = file->openDataSet(CDatasetName);
        dataset.write( &matValues[0], PredType::NATIVE_DOUBLE);
        
        dataset.close();
        dataspace.close();
      } else if(is<NumericVector>(DatasetValues) || is<IntegerVector>(DatasetValues)) 
      {
        
        hsize_t vectorsize;
        
        if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues)) 
          vectorsize = as<IntegerVector>(DatasetValues).length();
        else if (is<NumericVector>(DatasetValues))
          vectorsize = as<NumericVector>(DatasetValues).length();
        else 
          vectorsize = 1;
        
        hsize_t dims[] = {vectorsize};
        DataSpace dataspace(RANK1, dims);
        
        
        if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues) ) 
        {
          //..// int vectHiCValues[dims[0]];
          std::vector<int> vectHiCValues(dims[0]);
          
          for(int i=0;i<dims[0]; i++)
            vectHiCValues[i] = as<IntegerVector>(DatasetValues)(i);
          
          DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_INT, dataspace);
          dataset = file->openDataSet(CDatasetName);
          //..// dataset.write( vectHiCValues, PredType::NATIVE_INT);
          dataset.write( vectHiCValues.data(), PredType::NATIVE_INT);
          dataspace.close();
          dataset.close();
        } 
        else if(is<NumericVector>(DatasetValues) ) 
        {
          //..// double vectValues[dims[0]];
          std::vector<double> vectValues(dims[0]);
          for(int i=0;i<dims[0]; i++)
            vectValues[i] = as<NumericVector>(DatasetValues)(i);
          
          DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
          dataset = file->openDataSet(CDatasetName);
          //..// dataset.write(vectValues, PredType::NATIVE_DOUBLE);
          dataset.write(vectValues.data(), PredType::NATIVE_DOUBLE);
          dataspace.close();
          dataset.close();
        } 
        
      } 
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix_transposed_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix_transposed_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_HDF5_matrix_transposed_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_transposed_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_transposed_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  
  
  
  
  int write_HDF5_matrix_from_R_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues, bool transposed)
  {
    
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      // Create the data space for the dataset.
      std::vector<int> dims;
      
      if(is<NumericMatrix>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        if(transposed == true){
          dims[0] = as<NumericMatrix>(DatasetValues).rows();
          dims[1] = as<NumericMatrix>(DatasetValues).cols();
        } else{
          dims[0] = as<NumericMatrix>(DatasetValues).cols();
          dims[1] = as<NumericMatrix>(DatasetValues).rows();
        }
        
        DataSpace dataspace(RANK2, dims);
        
        DataSet dataset = file->createDataSet(CDatasetName,PredType::NATIVE_DOUBLE, dataspace);
        dataset = file->openDataSet(CDatasetName);
        
        if(transposed == true){
          std::vector<double> matHiCValues = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));
          dataset.write( &matHiCValues[0] , PredType::NATIVE_DOUBLE);
        }else{
          std::vector<double> matHiCValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));
          dataset.write( &matHiCValues[0] , PredType::NATIVE_DOUBLE);
        }
        
        dataset.close();
        dataspace.close();
        
      } 
      else if( Rcpp::is<IntegerMatrix>(DatasetValues) ) 
      {
        
        hsize_t dims[2];
        if(transposed == true){
          dims[0] = as<IntegerMatrix>(DatasetValues).rows();
          dims[1] = as<IntegerMatrix>(DatasetValues).cols();
        }else {
          dims[0] = as<IntegerMatrix>(DatasetValues).cols();
          dims[1] = as<IntegerMatrix>(DatasetValues).rows();
        }
        
        DataSpace dataspace(RANK2, dims);
        
        DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
        dataset = file->openDataSet(CDatasetName);
        
        if(transposed == true){
          std::vector<double> matHiCValues = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));
          dataset.write( &matHiCValues[0], PredType::NATIVE_DOUBLE);
        } else {
          std::vector<double> matHiCValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));
          dataset.write( &matHiCValues[0], PredType::NATIVE_DOUBLE);
        }
        
        dataset.close();
        dataspace.close();
      } else if(is<NumericVector>(DatasetValues) || is<IntegerVector>(DatasetValues)) 
      {
        
        hsize_t vectorsize;
        
        if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues)) 
          vectorsize = as<IntegerVector>(DatasetValues).length();
        else if (is<NumericVector>(DatasetValues))
          vectorsize = as<NumericVector>(DatasetValues).length();
        else 
          vectorsize = 1;
        
        hsize_t dims[] = {vectorsize};
        DataSpace dataspace(RANK1, dims);
        
        
        if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues) ) 
        {
          //..// int vectHiCValues[dims[0]];
          std::vector<int> vectHiCValues(dims[0]);
          for(int i=0;i<dims[0]; i++)
            vectHiCValues[i] = as<IntegerVector>(DatasetValues)(i);
          
          DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_INT, dataspace);
          dataset = file->openDataSet(CDatasetName);
          //..// dataset.write( vectHiCValues, PredType::NATIVE_INT);
          dataset.write( vectHiCValues.data(), PredType::NATIVE_INT);
          dataspace.close();
          dataset.close();
        } 
        else if(is<NumericVector>(DatasetValues) ) 
        {
          //..// double vectValues[dims[0]];
          std::vector<double> vectValues(dims[0]);
          for(int i=0;i<dims[0]; i++)
            vectValues[i] = as<NumericVector>(DatasetValues)(i);
          
          DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
          dataset = file->openDataSet(CDatasetName);
          //..// dataset.write(vectValues, PredType::NATIVE_DOUBLE);
          dataset.write(vectValues.data(), PredType::NATIVE_DOUBLE);
          dataspace.close();
          dataset.close();
        } 
        
      } else {
        Rcpp::Rcout<<"\n\n Unknown data type - Dataset not created !!\n";
      }
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  

  
/*  
  // Write dat from R to a HDF5 file by blocks ( faster? )
  int write_HDF5_matrix_by_blocks_from_R_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues, bool transposed)
  {
    
    IntegerVector stride = {1,1};
    IntegerVector block = {1,1};
    
    DataSet* dataset;
    
    
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      // Create the data space for the dataset.
      std::vector<int> dims;
      
      if(is<NumericMatrix>(DatasetValues)) 
      {
        
        hsize_t dims[2];
        
        int icols = as<NumericMatrix>(DatasetValues).cols(), 
            irows = as<NumericMatrix>(DatasetValues).rows();
        
        // Reserve space
        if(transposed == true){
          // dims[0] = as<NumericMatrix>(DatasetValues).rows();
          // dims[1] = as<NumericMatrix>(DatasetValues).cols();
          dims[0] = irows;
          dims[1] = icols;
        } else{
          // dims[0] = as<NumericMatrix>(DatasetValues).cols();
          // dims[1] = as<NumericMatrix>(DatasetValues).rows();
          dims[0] = icols;
          dims[1] = irows;
        }
        
        DataSpace dataspace(RANK2, dims);
        
        dataset = new DataSet(file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace));
        // dataset = file->openDataSet(CDatasetName);

        // Get block dimensions and number of blocks to read        
        int rowsInBlocks =  MAXELEMSINBLOCK / icols ;
        int nBlocks = irows / rowsInBlocks;

        
        for (int i=0; i<nBlocks; i++ )
        {
          int maxoffset;
          
          if( ((i+1) * rowsInBlocks ) < irows) { maxoffset = rowsInBlocks;
          } else { maxoffset = irows - ((i+1) * rowsInBlocks); }
          
          NumericMatrix BlockMatrix =  as<NumericMatrix>(DatasetValues)(Range(i*rowsInBlocks,(i+1)*rowsInBlocks), Range(0, icols ) );
          
          if(transposed == true){

            IntegerVector count = {icols, irows };
            IntegerVector offset = {icols, maxoffset};
            

            
            //.good.// std::vector<double> matHiCValues = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)(Range(i*rowsInBlocks, (i+1)*rowsInBlocks), 
            //.good.//                                                                                                        Range(0, icols ) )));
            //.good.// dataset.write( &matHiCValues[0] , PredType::NATIVE_DOUBLE);
            write_HDF5_matrix_subset_v2(file, dataset, offset, count, stride, block, transpose( BlockMatrix) );
              
          }else{
            //.good.// std::vector<double> matHiCValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));
            //.good.// dataset.write( &matHiCValues[0] , PredType::NATIVE_DOUBLE);
            
            IntegerVector count = {irows, icols };
            IntegerVector offset = {maxoffset, icols };
            
            write_HDF5_matrix_subset_v2(file, dataset, offset, count, stride, block, BlockMatrix);
          }
        }
        
        dataset->close();
        dataspace.close();
        
      } 
      else if( Rcpp::is<IntegerMatrix>(DatasetValues) ) 
      {
        
        int icols = as<IntegerMatrix>(DatasetValues).cols(), 
          irows = as<IntegerMatrix>(DatasetValues).rows();
        
        // Get block dimensions and number of blocks to read        
        int rowsInBlocks =  MAXELEMSINBLOCK / icols ;
        int nBlocks = irows / rowsInBlocks;
        
        hsize_t dims[2];
        if(transposed == true){
          // dims[0] = as<IntegerMatrix>(DatasetValues).rows();
          // dims[1] = as<IntegerMatrix>(DatasetValues).cols();
          dims[0] = irows;
          dims[1] = icols;
        }else {
          // dims[0] = as<IntegerMatrix>(DatasetValues).cols();
          // dims[1] = as<IntegerMatrix>(DatasetValues).rows();
          dims[0] = icols;
          dims[1] = irows;
        }
        
        DataSpace dataspace(RANK2, dims);
        
        dataset = new DataSet(file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace));
        // dataset = file->openDataSet(CDatasetName);
        
        for (int i=0; i<nBlocks; i++ )
        {
          int maxoffset;
          
          if( ((i+1) * rowsInBlocks ) < irows) { maxoffset = rowsInBlocks;
          } else { maxoffset = irows - ((i+1) * rowsInBlocks); }
          
          NumericMatrix BlockMatrix =  as<NumericMatrix>(DatasetValues)(Range(i*rowsInBlocks,(i+1)*rowsInBlocks), Range(0, icols ) );
        
          if(transposed == true){
            //.good.// std::vector<double> matHiCValues = as<std::vector<double> >(transpose(as<NumericMatrix>(DatasetValues)));
            //.good.// dataset.write( &matHiCValues[0], PredType::NATIVE_DOUBLE);
            
            IntegerVector count = {icols, irows };
            IntegerVector offset = {icols, maxoffset};
            
            write_HDF5_matrix_subset_v2(file, dataset, offset, count, stride, block, transpose(BlockMatrix));
            
          } else {
            //.good.// std::vector<double> matHiCValues = as<std::vector<double> >(as<NumericMatrix>(DatasetValues));
            //.good.// dataset.write( &matHiCValues[0], PredType::NATIVE_DOUBLE);
            
            IntegerVector count = {irows, icols };
            IntegerVector offset = {maxoffset, icols };
            
            write_HDF5_matrix_subset_v2(file, dataset, offset, count, stride, block, BlockMatrix);
          }
        }
        
        dataset->close();
        dataspace.close();
      } else if(is<NumericVector>(DatasetValues) || is<IntegerVector>(DatasetValues)) 
      {
        
        hsize_t vectorsize;
        
        if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues)) 
          vectorsize = as<IntegerVector>(DatasetValues).length();
        else if (is<NumericVector>(DatasetValues))
          vectorsize = as<NumericVector>(DatasetValues).length();
        else 
          vectorsize = 1;
        
        hsize_t dims[] = {vectorsize};
        DataSpace dataspace(RANK1, dims);
        
        
        if(is<IntegerVector>(DatasetValues) || is<LogicalVector>(DatasetValues) ) 
        {
          int vectHiCValues[dims[0]];
          for(int i=0;i<dims[0]; i++)
            vectHiCValues[i] = as<IntegerVector>(DatasetValues)(i);
          
          DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_INT, dataspace);
          dataset = file->openDataSet(CDatasetName);
          dataset.write( vectHiCValues, PredType::NATIVE_INT);
          dataspace.close();
          dataset.close();
        } 
        else if(is<NumericVector>(DatasetValues) ) 
        {
          double vectValues[dims[0]];
          for(int i=0;i<dims[0]; i++)
            vectValues[i] = as<NumericVector>(DatasetValues)(i);
          
          DataSet dataset = file->createDataSet(CDatasetName, PredType::NATIVE_DOUBLE, dataspace);
          dataset = file->openDataSet(CDatasetName);
          dataset.write(vectValues, PredType::NATIVE_DOUBLE);
          dataspace.close();
          dataset.close();
        } 
        
      } else {
        Rcpp::Rcout<<"\n\n Unknown data type - Dataset not created !!\n";
      }
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_from_R_ptr (Data TypeIException)" );
      return -1;
    }
    
    return 0;
  }
  
*/  
  
  
  

  int write_HDF5_matrix_subset_v2( H5File* file, DataSet* dataset,
                                   IntegerVector ivoffset, IntegerVector ivcount,
                                   IntegerVector ivstride, IntegerVector ivblock,
                                   RObject DatasetValues)
  {
    try
    {
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
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
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset_v2 (Data TypeIException)" );
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
    catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset (DataSet IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_HDF5_matrix_subset (DataSpace IException)" );
      return -1;
    }
    
    return 0;
  }
  
  
  // Create dataset from stringVector
  int write_hdf5_string_vector(H5File* file, std::string datasetname, StringVector DatasetValues)
  {
    
    //..// DataSet* unlimDataset;
    
    try
    {
      
      typedef struct name {
        char chr[MAXSTRING];
      } name;
      
      
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      // Open file
      //..// DataSet dataset = file->openDataSet(datasetname);
      
      // Create the data space for the dataset.
      hsize_t vectorsize;
      
      if (is<StringVector>(DatasetValues))
      {
        vectorsize = DatasetValues.length();
        
        // Define hdf5 dataspace size
        hsize_t dims[] = {vectorsize};
        DataSpace dataspace(RANK1, dims);
        
        // Create the memory datatype.
        H5::CompType mtype(sizeof(name));
        mtype.insertMember("chr", HOFFSET(name, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
        
        // Create the dataset.
        DataSet* dataset;
        dataset = new DataSet(file->createDataSet(datasetname, mtype, dataspace));
        
        // Get dataspace of the dataset.
        dataspace = dataset->getSpace();
        
        if(vectorsize > MAXSTRBLOCK) {

          // Number of blocks to process
          int iblocsks = vectorsize/MAXSTRBLOCK;
           
          for(int i=0; i<=iblocsks; i++)
          {
            
            // Gets block size to read
            hsize_t ilength = MAXSTRBLOCK;
            if(i == iblocsks){
              ilength = vectorsize - (i * MAXSTRBLOCK);
            }
            
            // Convert Dataframe to range list
            name *names_list = new name[ilength];
            
            for(int row=0; row< ilength; row++ )
            {
              String wchrom = as<StringVector>(DatasetValues)((i*MAXSTRBLOCK) + row);
              std::string word = wchrom.get_cstring();
              
              int j=0;
              for( j=0; j < word.size() && j < (MAXSTRING-1); j++ ){
                names_list[row].chr[j] = word[j]; }
              
              names_list[row].chr[j] = '\0'; // insert hdf5 end of string
            }
              
            // HyperSlab position and length
            hsize_t start[1];
            start[0] = (i*MAXSTRBLOCK);
            hsize_t count[] = {ilength};

            DataSpace memspace(RANK1, count, NULL);
            
            // Get position and write data in dataset
            dataspace.selectHyperslab(H5S_SELECT_SET, count, start); 
            dataset->write(names_list, mtype, memspace, dataspace);
            
            
            // Release resources
            delete[] names_list;
            memspace.close();

            
          }
          
        } else {

          int datarows = as<StringVector>(DatasetValues).size();
          
          // Convert Dataframe to range list
          name *names_list = new name[datarows];
          
          for(int i=0; i< datarows; i++ )
          {
            //..// name n;
            String wchrom = as<StringVector>(DatasetValues)(i);
            std::string word = wchrom.get_cstring();
            
            int j=0;
            for( j=0; j < word.size() && j < (MAXSTRING-1); j++ )
              names_list[i].chr[j] = word[j];
            
            names_list[i].chr[j] = '\0'; // insert hdf5 end of string
            
          }

          dataset->write(names_list, mtype);
          delete[] names_list;
        }
        
        // Release resources
        dataspace.close();
        dataset->close();
        delete dataset;
        
        
      }
      
    } 
    catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_hdf5_string_vector (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_hdf5_string_vector (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_hdf5_string_vector (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_hdf5_string_vector (DataSpace IException)" );
      return -1;
    }
    
    
    
    
    return 0;
    
  }
  
  
  
  
  

  
  // Write dimnames from Matrix RObjects to hdf5 data file 
  int write_hdf5_matrix_dimnames(H5File* file, std::string groupname, std::string datasetname, StringVector rownames, StringVector colnames )
  {
    
    try{
      
      
      Exception::dontPrint();
      
      std::string strGroup = groupname + "/." + datasetname + "_dimnames";
      
      // Add rownames
      create_HDF5_groups_ptr(file, strGroup);
      
      //..// create_HDF5_group_ptr(file, strGroup);
      if( rownames.length()>1 )
        write_hdf5_string_vector(file, strGroup + "/1" , rownames);
      // else
      //   Rcpp::Rcout<<"Info : no rownames to save";
      
      // Add colnames
      if( colnames.length()>1 )
        write_hdf5_string_vector(file, strGroup + "/2", colnames);
      // else
      //   Rcpp::Rcout<<"Info : no colnames to save";
      

    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception write_hdf5_matrix_dimnames (Data TypeIException)" );
      return -1;
    }
    
    return(0);
  }
  
  
  
  
  
  // Read rhdf5 data matrix subset, 
  // input : 
  //      ivoffset : start position
  //      ivcount : block size
  //      ivstride :(1,1) by default.
  //                Allows to sample elements along a dimension. A stride of one (or NULL) will select every element along a dimension,
  //                a stride of two will select every other element, and a stride of three will select an element after every two elements. 
  //      ivblock : (1,1) by default.
  //                The block array determines the size of the element block selected from a dataspace. If the block size is one or NULL 
  //                then the block size is a single element in that dimension.
  // output : 
  //    rdatablock : matrix block
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
    } catch(FileIException& error) { // catch failure caused by the H5File operations
        dataset->close();
        file->close();
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
        dataset->close();
        file->close();
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
        dataset->close();
        file->close();
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
        dataset->close();
        file->close();
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
        dataset->close();
        file->close();
      ::Rf_error( "c++ exception read_HDF5_matrix_subset (Data TypeIException)" );
      return -1;
    }
    
    return 0;  // successfully terminated
    
  }
  
  

  // Get mean and corrected sd from each column in dataset in the case of n<<m, this information is used
  // to normalize data, center or scale.
  int get_HDF5_mean_sd_by_column_ptr(H5File* file, DataSet* dataset, Eigen::MatrixXd& normalize )
  {
    
    IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
    // Eigen::MatrixXd meansd (2,dims_out[0]);
    
    try
    {
      
      // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
      Exception::dontPrint();
      
      int block_size = 500;
      
      IntegerVector stride = IntegerVector::create(1, 1);
      IntegerVector block = IntegerVector::create(1, 1);
      IntegerVector offset = IntegerVector::create(0, 0);
      IntegerVector count = IntegerVector::create(0, 0);
      
      
      // IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
      // Eigen::MatrixXd meansd (2,dims_out[0]);
      
      
      count[1] = dims_out[1];
      if( block_size < dims_out[0] )
        count[0] = block_size;
      else
        count[0] = dims_out[0];
      
      
      // Read data in blocks of 500 columns
      for(int i=0; (i < floor(dims_out[0]/block_size)) || i==0; i++)
      {

        if(i>0){
          offset[0] = offset[0] + block_size;
        
          if( offset[0] + block_size <= dims_out[0] ) {
            count[0] = block_size;
          }else {
            count[0] = dims_out[0] - offset[1]+block_size; 
          }
        }
        
        Eigen::MatrixXd X = GetCurrentBlock_hdf5(file, dataset, offset[0], offset[1], count[0], count[1]);

        Eigen::VectorXd mean = X.rowwise().mean();
        Eigen::VectorXd std = ((X.colwise() - mean).array().square().rowwise().sum() / (X.cols() - 1)).sqrt();

        normalize.block( 0, offset[0] , 1, mean.size()) = mean.transpose();
        normalize.block( 1, offset[0], 1, std.size()) = std.transpose();
        
      }
      
      // normalize = meansd;
      
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (Data TypeIException)" );
      return -1;
    }
    
    return(0);  // successfully terminated
    
  }
  

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


// Read dimnames from Hdf5 vector dataset 
StringVector get_hdf5_matrix_dimnames(H5File* file, std::string groupname, std::string datasetname, int idim )
{
  StringVector dimnames;
  
  try{
    
    Exception::dontPrint();
    
    std::string strGroup = "." + datasetname + "_dimnames";
    
    // Get dataset
    std::string strDataset = strGroup + "/" + std::to_string(idim);
    DataSet* dataset = new DataSet(file->openDataSet(strDataset));
    
    // Define data types
    H5::DataSpace dataspace  = dataset->getSpace();
    H5::StrType   datatype   = dataset->getStrType();
    
    // Get array dimension !!!
    IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
    
    for( size_t i = 0; i<dims_out[0]; i++ )
    {
      std::string field_value;
      dataset->read(field_value, datatype, dataspace);
      dimnames.push_back(field_value);
      
    }
    
    datatype.close();
    dataspace.close();
    
    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (File IException)" );
    return -1;
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (DataSet IException)" );
    return -1;
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (Group IException)" );
    return -1;
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (DataSpace IException)" );
    return -1;
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception get_hdf5_matrix_dimnames (Data TypeIException)" );
    return -1;
  }
  return(dimnames);
}




int Create_hdf5_matrix_unlimited_ptr( H5File* file, const H5std_string dataset , RObject mat)
{
  try
  {
    Exception::dontPrint();
    
    if ( mat.sexp_type()==0   )
      throw std::range_error("Data matrix must exsits and mustn't be null");
    
    Eigen::MatrixXd matdat = as<Eigen::MatrixXd>(mat);
    
    /* Create the data space with unlimited dimensions.   */
    hsize_t      dims[2]  = { (unsigned long long)as<NumericMatrix>(mat).rows(), 
                              (unsigned long long)as<NumericMatrix>(mat).cols()};  // dataset dimensions at creation
    hsize_t      maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
    DataSpace mspace1( RANK2, dims, maxdims);
    
    
    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_unlimited_ptr (File IException)" );
    return -1;
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_unlimited_ptr (DataSet IException)" );
    return -1;
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_unlimited_ptr (Group IException)" );
    return -1;
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_unlimited_ptr (DataSpace IException)" );
    return -1;
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_unlimited_ptr (Data TypeIException)" );
    return -1;
  }
  
  return 0;  // successfully terminated
  
  
}


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





//' Create hdf5 data file and write data to it
//'
//' Creates a hdf5 file with numerical data matrix,
//' 
//' @param filename, character array indicating the name of the file to create
//' @param object numerical data matrix
//' @param group, character array indicating folder name to put the matrix in hdf5 file
//' @param dataset, character array indicating the dataset name to store the matix data
//' @param transp boolean, if trans=true matrix is stored transposed in hdf5 file
//' @param force, optional boolean if true and file exists, removes old file and creates a new file with de dataset data.
//' @return none
//' @export
// [[Rcpp::export]]
void bdCreate_hdf5_matrix_file(std::string filename, RObject object, 
                                        Rcpp::Nullable<std::string> group = R_NilValue, Rcpp::Nullable<std::string> dataset = R_NilValue,
                                        Rcpp::Nullable<bool> transp = R_NilValue, Rcpp::Nullable<bool> force = R_NilValue )
{
  
  H5File* file;
  
  try
  {


    std::string strsubgroup, strdataset;
    bool transposed, bforce;
    
    if(group.isNull())  strsubgroup = "INPUT" ;
    else    strsubgroup = Rcpp::as<std::string>(group);
    

    if(dataset.isNull())  strdataset = "A" ;
    else    strdataset = Rcpp::as<std::string>(dataset);
    
    if(transp.isNull())  transposed = false ;
    else    transposed = Rcpp::as<bool>(transp);
    
    if(force.isNull())  bforce = false ;
    else    bforce = Rcpp::as<bool>(force);
    
    if ( object.sexp_type()==0   )
      throw std::range_error("Data matrix must exsits and mustn't be null");
    
    
    CharacterVector svrows, svrcols;
    List dimnames;
    
    if( ResFileExist_filestream(filename)) {
        if(bforce == true) {
            RemoveFile(filename);
        } else {
            Rcpp::Rcout<<"\n File already exits, please set force = TRUE if you want to overwrite the old file";
        }
    } 
    
    // Create HDF5 file
    file = new H5File( filename, H5F_ACC_TRUNC );
    create_HDF5_group_ptr(file, strsubgroup);

    if ( object.isS4() == true)    
    {
      //..// write_DelayedArray_to_hdf5(filename, strsubgroup + "/" + strdataset, mat);
      write_DelayedArray_to_hdf5_ptr(file, strsubgroup + "/" + strdataset, object, transposed);
      // Get attribs from DelayedArray Object
      RObject S4content = object.slot("seed");
      // Get dimnames attribs 
      dimnames = S4content.attr("dimnames");
      
    } else {
      
      try{  
        
        //DEBUG : Rcpp::Rcout<<"\n typeof : "<<TYPEOF(object)<<"\n";
        
        if ( TYPEOF(object) == INTSXP ) {
          //..// write_HDF5_matrix_ptr(file, strsubgroup + "/" + strdataset, Rcpp::as<IntegerMatrix>(object));
          write_HDF5_matrix_from_R_ptr(file, strsubgroup + "/" + strdataset, Rcpp::as<IntegerMatrix>(object), transposed);
        } else{
          //..// write_HDF5_matrix_ptr(file, strsubgroup + "/" + strdataset, Rcpp::as<NumericMatrix>(mat));
          write_HDF5_matrix_from_R_ptr(file, strsubgroup + "/" + strdataset, Rcpp::as<NumericMatrix>(object), transposed);
        }
      
        // Get dimnames from R matrix object
        dimnames = object.attr( "dimnames" );
        
      }
      catch(std::exception &ex) { }
    }
    
    // Write dimnaes from RObject to hdf5 data file
    if( (dimnames.size()>0) && (!Rf_isNull(dimnames)) ) {
      
      //..// Rcpp::Rcout<<"\n Dimnames [1] : "<<Rf_isNull(dimnames[1])<<"\n";
      //..// Rcpp::Rcout<<"\n Dimnames [0] : "<<Rf_isNull(dimnames[0])<<"\n";
      
      if(!Rf_isNull(dimnames[1]) ){  svrcols = dimnames[1]; }
      if(!Rf_isNull(dimnames[0]) ){ svrows = dimnames[0]; }


      if(transposed == false){
        write_hdf5_matrix_dimnames(file, strsubgroup, strdataset, svrows, svrcols );
      }else{
        write_hdf5_matrix_dimnames(file, strsubgroup, strdataset, svrcols, svrows );
      }
      
    }

    // Read dimnames and colnames from hdf5 data file (working ok)
    //..// StringVector rownames,colnames;
    //..// rownames = get_hdf5_matrix_dimnames(&file, strdataset, 1);
    //..// colnames = get_hdf5_matrix_dimnames(&file, strsubgroup, strdataset, 2);
    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_file (File IException)" );
    //..// return(wrap(-1));
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_file (DataSet IException)" );
    //..// return(wrap(-1));
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_file (Group IException)" );
    //..// return(wrap(-1));
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_file (DataSpace IException)" );
    //..// return(wrap(-1));
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception Create_hdf5_matrix_file (Data TypeIException)" );
    //..// return(wrap(-1));
  }
  file->close();
  Rcpp::Rcout<<"\nFile and dataset has been created\n";
  //..// return(wrap(0));

}


//' Write matrix to existing hdf5 file
//'
//' Creates a hdf5 file with numerical data matrix,
//' 
//' @param filename, character array indicating the name of the file to create
//' @param object numerical data matrix
//' @param group, character array indicating folder or group name to put the matrix in hdf5 file
//' @param dataset, character array indicating the dataset name that contains the matix data
//' @param transp, boolean if true, data is manipulated in transposed form
//' @return none
//' @export
// [[Rcpp::export]]
void bdAdd_hdf5_matrix(RObject object, std::string filename, std::string group, std::string dataset, Rcpp::Nullable<bool> transp = R_NilValue, Rcpp::Nullable<bool> force = false )
{
    
  H5File* file;
  try
  {

    if( object.sexp_type()==0   )
      throw std::range_error("Data matrix must exsits and mustn't be null");

    if(!ResFileExist(filename))
      throw std::range_error("File not exits, create file before add new dataset");

    CharacterVector svrows, svrcols;
    bool transposed,
         bforce;
    
    if(transp.isNull())  transposed = false ;
    else    transposed = Rcpp::as<bool>(transp);
    
    if(force.isNull())  bforce = false ;
    else    bforce = Rcpp::as<bool>(force);
    
    
    file = new H5File( filename, H5F_ACC_RDWR );
    
    if(!exists_HDF5_element_ptr(file, group)) 
      create_HDF5_group_ptr(file, group);
    
    if( exists_HDF5_element_ptr(file, group + "/" + dataset) ) {
        if(bforce == true) {
            remove_HDF5_element_ptr( file, group + "/" + dataset);
            
            // Remove rownames and colnames if exsits
            std::string strdimnames = group + "/." + dataset + "_dimnames";
            if(exists_HDF5_element_ptr(file, strdimnames)) {
                remove_HDF5_element_ptr( file, strdimnames);
            }
            
        } else {
            Rcpp::Rcout<<"\n ERROR - dataset exists, please set force = TRUE to force overwrite";
          //..// return(wrap(-1));
        }
    }

    //..04-07-2020..// file->close();

    if ( object.isS4() == true) {
      write_DelayedArray_to_hdf5_ptr(file, group + "/" + dataset, object, transposed);
    } else {
      
      try{  
        
        if ( TYPEOF(object) == INTSXP ) {
          write_HDF5_matrix_from_R_ptr(file, group + "/" + dataset, Rcpp::as<IntegerMatrix>(object), transposed);
        } else{
          write_HDF5_matrix_from_R_ptr(file, group + "/" + dataset, Rcpp::as<NumericMatrix>(object), transposed);
        }
      }
      catch(std::exception &ex) { }
    } 
    
    
    // Write dimnames to dimnames datasets in hdf5 file
    List dimnames = object.attr( "dimnames" );

    // Write dimnaes from RObject to hdf5 data file
    if( (dimnames.size()>0) && (!Rf_isNull(dimnames)) ) {
      
      if(!Rf_isNull(dimnames[1]) ){  svrcols = dimnames[1]; }
      if(!Rf_isNull(dimnames[0]) ){ svrows = dimnames[0]; }
      
      
      if(transposed == false){
        write_hdf5_matrix_dimnames(file, group, dataset, svrows, svrcols );
      }else{
        write_hdf5_matrix_dimnames(file, group, dataset, svrcols, svrows );
      }
    }

    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception Create_HDF5_matrix (File IException)" );
    //..// return(wrap(-1));
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception Create_HDF5_matrix (DataSet IException)" );
    //..// return(wrap(-1));
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception Create_HDF5_matrix (Group IException)" );
    //..// return(wrap(-1));
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception Create_HDF5_matrix (DataSpace IException)" );
    //..// return(wrap(-1));
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception Create_HDF5_matrix (Data TypeIException)" );
    //..// return(wrap(-1));
  }
  
  file->close();
  Rcpp::Rcout<<"\nDataset has been added to the file\n";
  //..// return(wrap(0));
  
}





  
//' Remove element group or dataset from  hdf5 file
//'
//' Remove group or dataset from  hdf5 file
//' 
//' @param filename, character array indicating the name of the file to create
//' @param element path to element, character array indicating the complete route to the element to be removed (folder or dataset). 
//' @return none
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdRemove_hdf5_element(std::string filename, std::string element)
{
  
  H5File* file;
  
  try
  {
    // int res;
    
    if(!ResFileExist(filename))
      throw std::range_error("File not exits, create file before add new dataset");
    
    file = new H5File( filename, H5F_ACC_RDWR );
    
    if(!exists_HDF5_element_ptr(file, element)) {
      file->close();
      throw std::range_error("Element not exits");
    } else{
      remove_HDF5_element_ptr(file, element);
    }
      
    file->close();
    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception Remove_hdf5_element (File IException)" );
    return(wrap(-1));
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception Remove_hdf5_element (DataSet IException)" );
    return(wrap(-1));
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception Remove_hdf5_element (Group IException)" );
    return(wrap(-1));
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception Remove_hdf5_element (DataSpace IException)" );
    return(wrap(-1));
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception Remove_hdf5_element (Data TypeIException)" );
    return(wrap(-1));
  }
  
  return(wrap(0));
  
}
  
  
  
  // Get dataset name inside a group
  StringVector get_dataset_names_from_group( H5File* file, std::string strgroup, std::string strprefix)
  {
    
    StringVector datasetnames;
    
    try{
      
      Exception::dontPrint();
      
      herr_t err;
      ssize_t len;
      hsize_t nobj;
      int otype;
      char memb_name[MAX_NAME];
      
      // get file id
      Group grp = file->openGroup(strgroup);
      hid_t gid = grp.getId();
      
      // get dataset names inside group
      err = H5Gget_num_objs(gid, &nobj);
      for (int i = 0; i < nobj; i++) 
      {
        len = H5Gget_objname_by_idx(gid, (hsize_t)i, memb_name, (size_t)MAX_NAME );
        
        otype =  H5Gget_objtype_by_idx(gid, (size_t)i );
        
        // 202109
        if( strprefix.compare("")!=0 ){
          if(otype == H5G_DATASET && (memb_name[0] == strprefix[0])) {
            datasetnames.push_back(memb_name);
          }
        } else {
          if(otype == H5G_DATASET ) {
            datasetnames.push_back(memb_name);
          }
        }
        
      }
    } catch(FileIException& error) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception get_dataset_names_from_group (File IException)" );
      return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception get_dataset_names_from_group (DataSet IException)" );
      return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
      ::Rf_error( "c++ exception get_dataset_names_from_group (Group IException)" );
      return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception get_dataset_names_from_group (DataSpace IException)" );
      return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      ::Rf_error( "c++ exception get_dataset_names_from_group (Data TypeIException)" );
      return -1;
    }
    
    return(datasetnames);
    
  }



// Get dataset name from pointer dataset
StringVector get_dataset_names_from_dataset_ptr( DataSet* dataset)
{
  
  StringVector datasetnames;
  
  try{
    
    Exception::dontPrint();
    
    std::string s = dataset->getObjName();
    std::size_t found = s.rfind("/");
    
    if (found!=std::string::npos)
      s.replace (0,(found+1),"");
    
    datasetnames.push_back(s);
    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception get_dataset_names_from_dataset_ptr (File IException)" );
    return -1;
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception get_dataset_names_from_dataset_ptr (DataSet IException)" );
    return -1;
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception get_dataset_names_from_dataset_ptr (Group IException)" );
    return -1;
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception get_dataset_names_from_dataset_ptr (DataSpace IException)" );
    return -1;
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception get_dataset_names_from_dataset_ptr (Data TypeIException)" );
    return -1;
  }
  
  return(datasetnames);
  
}




IntegerVector get_HDF5_dataset_size(DataSet dataset)
{
  // Get dataspace from dataset
  DataSpace dataspace = dataset.getSpace();
  IntegerVector dims;
  int ndims;
  
  // Get the number of dimensions in the dataspace.
  int rank = dataspace.getSimpleExtentNdims();
  
  // Get the dimension size of each dimension in the dataspace and
  // display them.
  hsize_t dims_out[2];
  ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
  
  if(rank==1)
    dims = IntegerVector::create( static_cast<int>(dims_out[0]), static_cast<int>(1));
  else
    dims = IntegerVector::create(static_cast<int>(dims_out[0]), static_cast<int>(dims_out[1]));
  
  return(dims);
}  
  
  
  
  
  
/*** TODO : 
 * Write slots
 *  Slot : path
 *  Objects ? Matrixs inside file??
 *  Size : Size of each object
 *  
 * Write all data inside a list for each object ??!!!
*/
  

/***R

*/