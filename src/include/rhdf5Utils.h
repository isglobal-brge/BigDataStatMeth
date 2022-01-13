#ifndef rhdf5Utils
#define rhdf5Utils

    #include <RcppEigen.h>
    #include <iostream>
    #include <string>
    #include <boost/algorithm/string.hpp>
    #include <fstream>
    #include <sys/stat.h>
    #include "H5Cpp.h"
    #include "ReadDelayedData.h"
    #include "hdf5_to_Eigen.h"
  

    // [[Rcpp::depends(RcppEigen)]]
    
    #define MAX_NAME 1024
    
    using namespace H5;
    using namespace Rcpp;
    
    const int	 RANK1 = 1;
    const int	 RANK2 = 2;
    const int	 RANK3 = 3;
    const int	DIM1 = 1;
    const int	DIM2 = 2;
    const int	DIM3 = 3;
    const int	MAXSTRING = 20;
    const hsize_t MAXSTRBLOCK = 100000;
    const hsize_t MAXELEMSINBLOCK = 250000;
    
    // a typedef for our managed H5File pointer
    typedef std::shared_ptr<H5::H5File> H5FilePtr;
    
    bool ResFileExist(const std::string& name);
    bool ResFileExist_filestream(std::string name);
    bool RemoveFile(std::string filename);
    
    //..// extern "C" StringVector get_dataset_names_from_group( H5File* file, std::string strgroup);
    StringVector get_dataset_names_from_group( H5File* file, std::string strgroup, std::string strprefix);
    StringVector get_dataset_names_from_dataset_ptr( DataSet* dataset);
    extern "C" int join_datasets(H5File* file, std::string strsubgroup, StringVector strinput, std::string strasout);
    
    extern "C" bool remove_HDF5_element_ptr(H5File* file, const H5std_string element);
    extern "C" bool remove_HDF5_multiple_elements_ptr(H5File* file, std::string strgroup, StringVector elements);
    extern "C" bool exists_HDF5_element_ptr(H5File* file, const H5std_string element);
    
    
    
    H5FilePtr Open_hdf5_file(const std::string& fname);
    extern "C" int create_HDF5_dataset(H5std_string filename, const std::string CDatasetName,
                                      const size_t rows, const size_t cols, std::string strdatatype);
    extern "C" int create_HDF5_dataset_ptr(H5File* file, const std::string CDatasetName, 
                                const size_t rows, const size_t cols, std::string strdatatype);
    extern "C" int create_HDF5_unlimited_matrix_dataset_ptr(H5File* file, const std::string CDatasetName, 
                                          const size_t rows, const size_t cols, std::string strdatatype);
    extern "C" int create_HDF5_unlimited_vector_dataset_ptr(H5File* file, const std::string CDatasetName, 
                                                           const size_t length, std::string strdatatype);
    
    extern "C" int extend_HDF5_matrix_subset_ptr(H5File* file, DataSet* dataset, const size_t rows, const size_t cols);
    extern "C" int extend_HDF5_vector_subset_ptr(H5File* file, DataSet* dataset, const size_t length);
      
    extern "C" int create_HDF5_group(H5std_string filename, const H5std_string hiCGroup);
    extern "C" int create_HDF5_group_ptr( H5File* file, const H5std_string mGroup);
    extern "C" int create_HDF5_groups_ptr( H5File* file, const H5std_string mGroup);
    
    extern "C" int get_HDF5_mean_sd_by_column_ptr(H5File* file, DataSet* dataset, Eigen::MatrixXd& normalize );
    
    extern "C" int Create_hdf5_file(std::string filename);
    extern "C" int create_HDF5_matrix(H5std_string filename, const std::string DatasetName, RObject DatasetValues);
    /*extern "C" int write_HDF5_matrix(H5std_string filename, const std::string CDatasetName, RObject DatasetValues);*/
    extern "C" int write_HDF5_matrix_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues);
    extern "C" int write_HDF5_matrix_transposed_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues);
    extern "C" int write_HDF5_matrix_from_R_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues, bool transposed);
    extern "C" int write_HDF5_matrix_by_blocks_from_R_ptr(H5File* file, const std::string CDatasetName, RObject DatasetValues, bool transposed);
    
    
    extern "C" int write_HDF5_matrix_subset(H5std_string filename, const std::string CDatasetName, 
                                           IntegerVector ivoffset, IntegerVector ivcount,
                                           IntegerVector ivstride, IntegerVector ivblock,
                                           RObject DatasetValues);
    
    extern "C" int write_HDF5_matrix_subset_v2( H5File* file, DataSet* dataset,
                                               IntegerVector ivoffset, IntegerVector ivcount,
                                               IntegerVector ivstride, IntegerVector ivblock,
                                               RObject DatasetValues);
    
    extern "C" int write_hdf5_string_vector(H5File* file, std::string datasetname, StringVector DatasetValues);
    StringVector get_hdf5_matrix_dimnames(H5File* file, std::string groupname, std::string datasetname, int idim );
    extern "C" int write_hdf5_matrix_dimnames(H5File* file, std::string groupname, std::string datasetname, 
                                             StringVector rownames, StringVector colnames );
    
    extern "C" int read_HDF5_matrix_subset(H5File* file, DataSet* dataset,
                                           IntegerVector ivoffset, IntegerVector ivcount,
                                           IntegerVector ivstride, IntegerVector ivblock,
                                           double* rdatablock);
    
    IntegerVector get_HDF5_dataset_size(DataSet dataset);
    
    void create_symLink( H5File* file, std::string original, std::string link);
    void create_hardLink( H5File* file, std::string original, std::string link);


#endif