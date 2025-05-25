/**
 * @file hdf5Datasets.hpp
 * @brief HDF5 dataset management class and utilities
 * @details This header file provides a comprehensive class for managing HDF5 datasets
 * and related operations. The implementation includes:
 * 
 * Key features:
 * - Dataset creation and management
 * - Data reading and writing operations
 * - Block-based data access
 * - Attribute handling
 * - Support for various data types
 * - Unlimited dataset support
 * 
 * The class supports:
 * - Integer, numeric (double), and string datasets
 * - Row and column-major data layouts
 * - Block-based read/write operations
 * - Dataset attributes
 * - Automatic resource management
 */

#ifndef BIGDATASTATMETH_HDF5_DATASETS_HPP
#define BIGDATASTATMETH_HDF5_DATASETS_HPP

#include "Utilities/Utilities.hpp"
#include <hdf5.h>


namespace BigDataStatMeth {

/**
 * @class hdf5Dataset
 * @brief Class for managing HDF5 datasets
 * @details Provides comprehensive functionality for creating, reading, writing,
 * and managing HDF5 datasets. Inherits from hdf5Group to handle group operations.
 * 
 * Key capabilities:
 * - Dataset creation with various data types
 * - Block-based data access
 * - Attribute management
 * - Resource cleanup
 * - Dimension handling
 */
class hdf5Dataset : public hdf5Group
{
    
public:
    /**
     * @brief Constructor with file, group, and dataset name
     * @param file HDF5 file pointer
     * @param group Group path
     * @param datasetname Dataset name
     * @param overwrite Whether to overwrite existing dataset
     */
    hdf5Dataset(H5::H5File* file, std::string group, std::string datasetname, bool overwrite) : 
    hdf5Group(file, group)
    {
        name = datasetname;
        boverwrite = overwrite;
    }
    
    hdf5Dataset(H5::H5File* file, std::string dataset, bool overwrite) : 
        hdf5Group(file, SplitElementName(dataset).path)
    {
        fullpath datasetroute = SplitElementName(dataset);
        name = datasetroute.filename;
        boverwrite = overwrite;
    }
    
    hdf5Dataset(BigDataStatMeth::hdf5File* oFile, std::string group, std::string datasetname, bool overwrite) : 
    hdf5Group(oFile, group)
    {
        groupname = group;
        name = datasetname;
        boverwrite = overwrite;
    }
    
    hdf5Dataset(std::string filename, std::string dataset, bool overwrite) : 
    hdf5Group(filename, SplitElementName(dataset).path)
    {
        fullpath datasetroute = SplitElementName(dataset);
        name = datasetroute.filename;
        boverwrite = overwrite;
        
    }
    
    
    hdf5Dataset(std::string filename, std::string group, std::string datasetname, bool overwrite) : 
    hdf5Group(filename, group)
    {
        name = datasetname;
        boverwrite = overwrite;
    }
    
    
    /**
     * @brief Remove the dataset
     * @details Deletes the dataset from the HDF5 file
     */
    virtual void remove() {
        remove_elements(pfile, getGroupName(), {name}); 
    }
    
    
    /**
     * @brief Create a new dataset
     * @details Creates a new HDF5 dataset with specified dimensions and data type.
     * Supported data types:
     * - "int": Integer dataset
     * - "numeric" or "real": Double dataset
     * - "string": String dataset
     * 
     * @param rows Number of rows
     * @param cols Number of columns
     * @param strdatatype Data type ("int", "numeric", "real", or "string")
     */
    virtual void createDataset(size_t rows, size_t cols, std::string strdatatype) 
    {
        
        try
        {

            H5::Exception::dontPrint();
            std::string fullDatasetPath = groupname + "/" + name;
            bool bRemoved = false;
            
            // dataset dimensions
            dimDataset[0] = cols;
            dimDataset[1] = rows;
            
            dimDatasetinFile[0] = rows;
            dimDatasetinFile[1] = cols;
            
            if( !exists_HDF5_element(pfile, groupname) ) {
                create_HDF5_groups(groupname);
            }
            
            H5::DataSpace dataspace( RANK2, dimDataset );
            
            bool bexists = exists_HDF5_element(pfile, fullDatasetPath);
            if( bexists == true && boverwrite == false) {
                Rcpp::Rcerr<<"\nDataset exits, please set overwrite = true to overwrite the existing dataset (DataSet IException)";
                return void();
            } else {
                
                if( boverwrite == true && bexists == true) {
                    remove_elements(pfile, getGroupName(), {name}); 
                    bRemoved = true;
                }
                
                type = strdatatype;
                
                if( type == "string") {
                    H5::CompType strtype(sizeof(names));
                    strtype.insertMember("chr", HOFFSET(names, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
                    pdataset = new H5::DataSet(pfile->createDataSet(fullDatasetPath, strtype, dataspace));
                } else if( type == "int" || type == "logic" || type == "factor") {
                    H5::IntType datatype( H5::PredType::NATIVE_INT );
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace ));
                    if(bRemoved == true) {
                        writeDataset(Rcpp::wrap(Eigen::MatrixXd::Zero(dimDatasetinFile[0], dimDatasetinFile[1]) ));
                    }
                } else if( type == "numeric" || type == "real") {
                    H5::IntType datatype( H5::PredType::NATIVE_DOUBLE ); 
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace ));
                    if(bRemoved == true) {
                        writeDataset(Rcpp::wrap(Eigen::MatrixXd::Zero(dimDatasetinFile[0], dimDatasetinFile[1]) ));
                    }
                } else {
                    close_file();
                    Rcpp::Rcerr<<"\nDataset data type not allowed or no matrix defined (createDataset)";
                }
            }
            
            dataspace.close();
            
            addAttribute( "internal", Rcpp::wrap("0") );
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createDataset (File IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createDataset (Group IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createDataset (DataSet IException)";
        } 
        return void();
    }
    
    
    /**
     * @brief Create a dataset based on another dataset's dimensions
     * @details Creates a new dataset with the same dimensions as the reference dataset
     * but with a specified data type.
     * 
     * @param dsLike Reference dataset to copy dimensions from
     * @param strdatatype Data type for the new dataset
     */
    virtual void createDataset(BigDataStatMeth::hdf5Dataset* dsLike, std::string strdatatype) 
    {
        try{
            
            createDataset( dsLike->ncols(), dsLike->nrows(), strdatatype);
        } catch(H5::FileIException& error) {
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createDataset (File IException)";
        } catch(H5::GroupIException& error) {
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createDataset (Group IException)";
        } catch(H5::DataSetIException& error) {
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createDataset (DataSet IException)";
        } 
        
        return void();
    }
    
    /**
     * @brief Create an unlimited dataset
     * @details Creates a new HDF5 dataset with unlimited dimensions, allowing it to
     * grow in size. The dataset is created with initial dimensions but can be extended.
     * 
     * Features:
     * - Unlimited dimensions in both directions
     * - Chunked storage for efficiency
     * - Configurable initial size
     * 
     * @param rows Initial number of rows
     * @param cols Initial number of columns
     * @param strdatatype Data type for the dataset
     */
    virtual void createUnlimitedDataset(size_t rows, size_t cols, std::string strdatatype) 
    {
        try
        {
            H5::Exception::dontPrint();
            
            herr_t status;
            hid_t cparms; 
            std::string fullDatasetPath = groupname + "/" + name;
            
            
            // dataset dimensions
            dimDataset[0] = cols;
            dimDataset[1] = rows;
            
            // file dimensions
            dimDatasetinFile[0] = rows;
            dimDatasetinFile[1] = cols;
            
            // set dataset as unlimited;
            unlimited = true;
            
            // Declare unlimited dimensions
            hsize_t  maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
            H5::DataSpace dataspace ( RANK2, dimDataset, maxdims );
            
            // Enabling chunking
            hsize_t chunk_dims[2];
            chunk_dims[0] = rows;
            chunk_dims[1] = cols;
            
            cparms = H5Pcreate(H5P_DATASET_CREATE);
            status = H5Pset_chunk( cparms, RANK2, chunk_dims);
            
            if(status<0) {
                Rcpp::Rcerr<<"\nc++ exception createUnlimitedDataset (setting chunk IException)";
                return void();
            } 
            
            if( !exists_HDF5_element(pfile, groupname) ) {
                create_HDF5_groups(groupname);
            }
            
            bool bexists = exists_HDF5_element(pfile, fullDatasetPath);
            if( bexists == true && boverwrite == false) {
                Rcpp::Rcout<<"\n Dataset exits, please set overwrite = true to overwrite data \n";
            } else {
                
                if(  bexists == true && boverwrite == true) {
                    remove_elements(pfile, fullDatasetPath); 
                }
                
                // Create dataset
                if( strdatatype == "int") {
                    H5::IntType datatype( H5::PredType::NATIVE_INT );
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace, cparms));
                } else {
                    H5::IntType datatype( H5::PredType::NATIVE_DOUBLE ); 
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace, cparms));
                }
            }
            
            dataspace.close();
            addAttribute( "internal", Rcpp::wrap("0") );
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createUnlimitedDataset (File IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createUnlimitedDataset (File IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception createUnlimitedDataset (File IException)";
        } 
        return void();
    }
    
    
    /**
     * @brief Extend an unlimited dataset
     * @details Increases the dimensions of an unlimited dataset to accommodate
     * more data. Only works with datasets created as unlimited.
     * 
     * @param rows New number of rows
     * @param cols New number of columns
     * @throws H5::DataSetIException if dataset is not unlimited
     */
    virtual void extendUnlimitedDataset(const size_t rows, const size_t cols)
    {
        try
        {
            if(unlimited == true) {
                H5::Exception::dontPrint();
                
                // Extend dataset size to:  oldDims + newDims
                hsize_t newdims[2];
                newdims[0] = cols;
                newdims[1] = rows;
                
                // hsize_t size[2];
                dimDataset[0] = dimDataset[0] + newdims[0];
                dimDataset[1] = dimDataset[1] + newdims[1];
                
                pdataset->extend( dimDataset );    
            } else {
                Rcpp::Rcerr<<"\n Dataset is not an unlimited dataset, fixed datasets can't be extended\n";
                return void();
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception extend_HDF5_matrix_subset_ptr (File IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception extend_HDF5_matrix_subset_ptr (DataSet IException)";
        }
        return void();
    }
    
    
    // Open an existing hdf5 dataSet
    H5::DataSet* openDataset()
    {
        try
        {
            H5::Exception::dontPrint();
            std::string fullPath = groupname + "/" + name;
            
            // Check if file pointer != nullptr
            if( !pfile)  {
                Rcpp::Rcerr<<"\nc++ exception Please create file before proceed";
            } else { 
                bool bexists = exists_HDF5_element(pfile, fullPath);
                if( bexists ) {
                    if ((pdataset == NULL) == TRUE) {
                        pdataset = new H5::DataSet(pfile->openDataSet(fullPath));    
                    }
                    getDimensExistingDataset();
                    getAttribute("internal");
                } else {
                    close_file();
                    // std::cerr<<"\nc++ exception, please create Dataset before proceed\n";
                    Rcpp::Rcerr<<"\nc++ exception, please create Dataset before proceed";
                    // return(pdataset);
                }
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception openDataset (File IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception openDataset (File GroupIException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the H5File operations
            close_file();
            Rcpp::Rcerr<<"\nc++ exception openDataset (File DataSetIException)";
        } 
        
        return(pdataset);
    }
    
    
    
    /**
     * @brief Write data to the dataset
     * @details Writes data to the dataset, handling different data types and
     * formats. Supports:
     * - Numeric matrices
     * - Integer matrices
     * - String data
     * - Vector data
     * 
     * @param DatasetValues R object containing the data to write
     * @throws H5::DataSetIException if write operation fails
     */
    virtual void writeDataset( Rcpp::RObject DatasetValues )
    {
        
        try
        {
            H5::Exception::dontPrint();
            std::vector<int> dims;
            
            if(Rcpp::is<Rcpp::NumericMatrix>(DatasetValues) || Rcpp::is<Rcpp::IntegerMatrix>(DatasetValues) ) {
                
                hsize_t ncolRObject = Rcpp::as<Rcpp::NumericMatrix>(DatasetValues).ncol();
                hsize_t nrowRObject = Rcpp::as<Rcpp::NumericMatrix>(DatasetValues).nrow();
                
                if( (nrowRObject == dimDataset[1] && ncolRObject == dimDataset[0]) || ncolRObject == 1 || nrowRObject == 1 ) {
                    if( nrowRObject * ncolRObject > dimDataset[0] * dimDataset[1]) {
                        Rcpp::Rcout<<"\n Data you are trying to write is bigger than existing hdf5 dataset size\n";
                        return void();
                    } else {
                        H5::DataSpace dataspace(RANK2, dimDataset);
                        
                        std::vector<double> matHiCValues = Rcpp::as<std::vector<double> >(Rcpp::as<Rcpp::NumericMatrix>(DatasetValues));
                        pdataset->write( &matHiCValues[0] , H5::PredType::NATIVE_DOUBLE);
                        dataspace.close();    
                    }
                } else {
                    Rcpp::Rcout<<"\n Data you are trying to write (a complete dataset) differs from existing hdf5 dataset size\n";
                    return void();
                }
                
            } else if(Rcpp::is<Rcpp::NumericVector>(DatasetValues) || Rcpp::is<Rcpp::IntegerVector>(DatasetValues)) {
                
                hsize_t dims[] = {dimDataset[1]};
                H5::DataSpace dataspace(RANK1, dims);
                
                if(Rcpp::is<Rcpp::IntegerVector>(DatasetValues) || Rcpp::is<Rcpp::LogicalVector>(DatasetValues) ) {
                    std::vector<int> vectHiCValues = Rcpp::as<std::vector<int> >(DatasetValues);
                    pdataset->write( vectHiCValues.data() , H5::PredType::NATIVE_INT);
                } else if(Rcpp::is<Rcpp::NumericVector>(DatasetValues) )  {
                    std::vector<double> vectHiCValues = Rcpp::as<std::vector<double> >(DatasetValues);
                    pdataset->write( vectHiCValues.data() , H5::PredType::NATIVE_DOUBLE);
                } 
                dataspace.close();
            } else if(Rcpp::is<Rcpp::StringVector >(DatasetValues) || Rcpp::is<Rcpp::StringMatrix >(DatasetValues)) {
                
                // Create the memory datatype.
                H5::CompType strtype(sizeof(names));
                strtype.insertMember("chr", HOFFSET(names, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
                
                pdataset->write(convert_DataFrame_to_RangeList(DatasetValues, true), strtype);
                    
            } else {
                Rcpp::Rcerr<<"\nMatrix data type not allowed (writeDataset)";
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (File IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (DataSet IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (Group IException)";
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (DataSpace IException)";
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (Data TypeIException)";
        }
        return void();
    }
    
    
    
    
    
    /**
     * @brief Write raw double data to dataset
     * @details Direct write of double array to dataset without type conversion.
     * 
     * @param mdata Pointer to double array
     * @throws H5::DataSetIException if write operation fails
     */
    void writeDataset ( double* mdata)
    {
        try
        {
            
            H5::Exception::dontPrint();

            std::vector<int> dims;
            
            if( sizeof(mdata) > dimDataset[0] * dimDataset[1]) {
                Rcpp::Rcout<<"\n Data you are trying to write is bigger than existing hdf5 dataset size\n";
                return void();
            } else {
                
                H5::DataSpace dataspace(RANK2, dimDataset);
                // std::vector<double> matHiCValues = Rcpp::as<std::vector<double> >(Rcpp::as<Rcpp::NumericMatrix>(DatasetValues));
                pdataset->write( &mdata[0] , H5::PredType::NATIVE_DOUBLE);
                dataspace.close();
                
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (File IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (DataSet IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (Group IException)";
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (DataSpace IException)";
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDataset (Data TypeIException)";
        }
        return void();
    }
    
    
    

    // SPECIFIC FUNCTIONS TO WORK WITH EIGEN OBJECTS -- RowMajor Matrix
    /**
     * @brief Write block of data in row-major order
     * @details Writes a block of data to the dataset using row-major memory layout.
     * Supports:
     * - Partial dataset updates
     * - Strided access
     * - Block-based operations
     * 
     * @param DatasetValues Matrix data in row-major format
     * @param vOffset Starting position for write
     * @param vCount Number of elements to write
     * @param vStride Stride between elements
     * @param vBlock Size of blocks
     */
    virtual void writeRowMajorDatasetBlock( Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> DatasetValues, 
                                    std::vector<hsize_t> vOffset,  std::vector<hsize_t> vCount, 
                                    std::vector<hsize_t> vStride, std::vector<hsize_t> vBlock)
    {
        try
        {
            // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
            H5::Exception::dontPrint();
            
            hsize_t hsOffset[2], hsCount[2], hsStride[2], hsBlock[2];
                
            if(vOffset[0] + vCount[0] <= dimDataset[0] || vOffset[1] + vCount[1] <= dimDataset[1]) {
                
                // Specify size and shape of subset to write
                hsOffset[0] = vOffset[0]; hsOffset[1] = vOffset[1];
                hsStride[0] = vStride[0]; hsStride[1] = vStride[1]; // default 1
                hsBlock[0] = vBlock[0]; hsBlock[1] = vBlock[1]; // default 1
                
                hsCount[0] = vCount[0]; hsCount[1] = vCount[1];
                
                H5::DataSpace dataspace(RANK2, hsCount);
                H5::DataSpace memspace(RANK2, hsCount, NULL);
                
                dataspace = pdataset->getSpace();
                dataspace.selectHyperslab( H5S_SELECT_SET, hsCount, hsOffset, hsStride, hsBlock);
                
                pdataset->write(&DatasetValues.data()[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
                memspace.close();
                dataspace.close();
                    
            } else {
                close_dataset_file();
                Rcpp::Rcerr<<"\nIt is not possible to write block in current position (writeRowMajorDatasetBlock)";
            }
                
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeRowMajorDatasetBlock (File IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeRowMajorDatasetBlock (DataSet IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeRowMajorDatasetBlock (Group IException)";
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeRowMajorDatasetBlock (DataSpace IException)";
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeRowMajorDatasetBlock (Data TypeIException)";
        }
        return void();
    }
    
    
    // SPECIFIC FUNCTIONS TO WORK WITH EIGEN OBJECTS -- RowMajor Matrix
    /**
     * @brief Write block of data in column-major order
     * @details Writes a block of data to the dataset using column-major memory layout.
     * Supports:
     * - Partial dataset updates
     * - Strided access
     * - Block-based operations
     * 
     * @param DatasetValues Matrix data in column-major format
     * @param vOffset Starting position for write
     * @param vStride Stride between elements
     * @param vBlock Size of blocks
     */
    virtual void writeColMajorDatasetBlock( Eigen::MatrixXd DatasetValues, 
                                            std::vector<hsize_t> vOffset, std::vector<hsize_t> vStride, std::vector<hsize_t> vBlock)
    {
        try
        {
            // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
            H5::Exception::dontPrint();
            
            hsize_t hsOffset[2], hsCount[2], hsStride[2], hsBlock[2];
            
            if(vOffset[0] + DatasetValues.rows() <= dimDataset[0] || vOffset[1] + DatasetValues.cols() <= dimDataset[1]) {
                
                // Specify size and shape of subset to write
                hsOffset[0] = vOffset[0]; hsOffset[1] = vOffset[1];
                hsStride[0] = vStride[0]; hsStride[1] = vStride[1]; // default 1
                hsBlock[0] = vBlock[0]; hsBlock[1] = vBlock[1]; // default 1
                
                hsCount[0] = DatasetValues.rows(); hsCount[1] = DatasetValues.cols();
                
                H5::DataSpace dataspace(RANK2, hsCount);
                H5::DataSpace memspace(RANK2, hsCount, NULL);
                
                dataspace = pdataset->getSpace();
                dataspace.selectHyperslab( H5S_SELECT_SET, hsCount, hsOffset, hsStride, hsBlock);
                
                pdataset->write(&DatasetValues.data()[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
                memspace.close();
                dataspace.close();
                
            } else {
                Rcpp::Rcerr<<"\nIt is not possible to write block in current position (writeColMajorDatasetBlock)";
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeColMajorDatasetBlock (File IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeColMajorDatasetBlock (DataSet IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeColMajorDatasetBlock (Group IException)";
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeColMajorDatasetBlock (DataSpace IException)";
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeColMajorDatasetBlock (Data TypeIException)";
        }
        return void();
    }
    
    /**
     * @brief Write block of data with optional transposition
     * @details Writes a block of data to the dataset with configurable layout and
     * optional transposition. Supports:
     * - Row/column-major formats
     * - Data transposition
     * - Block-based access
     * 
     * @param DatasetValues Data to write
     * @param vOffset Starting position
     * @param vCount Number of elements
     * @param vStride Stride between elements
     * @param vBlock Block size
     * @param bTranspose Whether to transpose data
     */
    virtual void writeDatasetBlock( Rcpp::RObject DatasetValues, std::vector<hsize_t> vOffset, 
                            std::vector<hsize_t> vCount, std::vector<hsize_t> vStride,
                            std::vector<hsize_t> vBlock, bool bTranspose)
    {
        try
        {
            // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
            H5::Exception::dontPrint();

            hsize_t hsOffset[2], hsCount[2], hsStride[2], hsBlock[2];

            // Specify size and shape of subset to write
            hsOffset[0] = vOffset[0]; hsOffset[1] = vOffset[1];
            hsStride[0] = vStride[0]; hsStride[1] = vStride[1]; // default 1
            hsBlock[0] = vBlock[0]; hsBlock[1] = vBlock[1]; // default 1
            
            if(Rcpp::is<Rcpp::NumericMatrix>(DatasetValues) || Rcpp::is<Rcpp::IntegerMatrix>(DatasetValues) ||
               Rcpp::is<Rcpp::NumericVector>(DatasetValues) || Rcpp::is<Rcpp::IntegerVector>(DatasetValues) ) {
                
                if(Rcpp::is<Rcpp::NumericMatrix>(DatasetValues) || Rcpp::is<Rcpp::IntegerMatrix>(DatasetValues)) {
                    if(bTranspose == false) {
                        hsCount[0] = vCount[0];
                        hsCount[1] = vCount[1];
                    } else {
                        hsCount[0] = vCount[1];
                        hsCount[1] = vCount[0];
                    }
                    
                } else if(Rcpp::is<Rcpp::NumericVector>(DatasetValues) || Rcpp::is<Rcpp::IntegerVector>(DatasetValues)) {
                        hsCount[0] = 1;
                        hsCount[1] = vCount[0];
                }

                if(vOffset[0] + hsCount[0] <= dimDataset[0] || vOffset[1] + hsCount[1] <= dimDataset[1]) {
                    H5::DataSpace dataspace(RANK2, hsCount);
                    H5::DataSpace memspace(RANK2, hsCount, NULL);

                    dataspace = pdataset->getSpace();
                    dataspace.selectHyperslab( H5S_SELECT_SET, hsCount, hsOffset, hsStride, hsBlock);
                    
                    if(Rcpp::is<Rcpp::NumericMatrix>(DatasetValues) || Rcpp::is<Rcpp::IntegerMatrix>(DatasetValues)) {
                        std::vector<double> matdata(hsCount[0]*hsCount[1]);

                        if (bTranspose == true) {
                            matdata = Rcpp::as<std::vector<double> >(transpose(Rcpp::as<Rcpp::NumericMatrix>(DatasetValues)));
                        } else {
                            matdata = Rcpp::as<std::vector<double> >(Rcpp::as<Rcpp::NumericMatrix>(DatasetValues));
                        }
                        
                        pdataset->write(&matdata[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
                        memspace.close();
                        dataspace.close();
                    } else if(Rcpp::is<Rcpp::NumericVector>(DatasetValues) || Rcpp::is<Rcpp::IntegerVector>(DatasetValues)) {
                        std::vector<double> matdata = Rcpp::as<std::vector<double> >(Rcpp::as<Rcpp::NumericVector>(DatasetValues));
                        pdataset->write(&matdata[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
                        memspace.close();
                        dataspace.close();
                    }
                } else {
                    Rcpp::Rcerr<<"\nIt is not possible to write block in current position (writeDatasetBlock)";
                }

            } else if(Rcpp::is<Rcpp::StringMatrix>(DatasetValues) || Rcpp::is<Rcpp::StringVector>(DatasetValues) ) {
                if(Rcpp::is<Rcpp::StringMatrix>(DatasetValues)) {
                    hsCount[0] = Rcpp::as<Rcpp::StringMatrix>(DatasetValues).rows();
                    hsCount[1] = Rcpp::as<Rcpp::StringMatrix>(DatasetValues).cols();
                } else if(Rcpp::is<Rcpp::StringVector>(DatasetValues)) {
                    hsCount[0] = 1;
                    hsCount[1] = Rcpp::as<Rcpp::StringVector>(DatasetValues).length();
                }

                if(vOffset[0] + hsCount[0] <= dimDataset[0] || vOffset[1] + hsCount[1] <= dimDataset[1]) {

                    // Create the memory datatype.
                    H5::CompType strtype(sizeof(names));
                    strtype.insertMember("chr", HOFFSET(names, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));

                    H5::DataSpace dataspace(RANK2, hsCount);
                    H5::DataSpace memspace(RANK2, hsCount, NULL);

                    dataspace = pdataset->getSpace();
                    dataspace.selectHyperslab( H5S_SELECT_SET, hsCount, hsOffset, hsStride, hsBlock);

                    pdataset->write(convert_DataFrame_to_RangeList(DatasetValues, false), strtype, memspace, dataspace);
                    memspace.close();
                    dataspace.close();

                } else {
                    Rcpp::Rcerr<<"\nIt is not possible to write block in current position (writeDatasetBlock)";
                }
            } else {
                Rcpp::Rcerr<<"\nMatrix data type not allowed (writeDatasetBlock)";
            }

        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock (File IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock (DataSet IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock (Group IException)";
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock (DataSpace IException)";
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock (Data TypeIException)";
        }
        return void();
    }

    /**
     * @brief Write block of vector data
     * @details Writes a block of vector data to the dataset. Useful for
     * column/row-wise operations.
     * 
     * @param DatasetValues Vector data to write
     * @param vOffset Starting position
     * @param vCount Number of elements
     * @param vStride Stride between elements
     * @param vBlock Block size
     */
    virtual void writeDatasetBlock( std::vector<double> DatasetValues, std::vector<hsize_t> vOffset, 
                                    std::vector<hsize_t> vCount, std::vector<hsize_t> vStride,
                                    std::vector<hsize_t> vBlock)
    {
        try
        {
            // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
            H5::Exception::dontPrint();
            
            hsize_t hsOffset[2], hsCount[2], hsStride[2], hsBlock[2];
            
            // Specify size and shape of subset to write
            hsStride[0] = vStride[0]; hsStride[1] = vStride[1]; // default 1
            hsBlock[0] = vBlock[0]; hsBlock[1] = vBlock[1]; // default 1
                
            hsOffset[0] = vOffset[0]; hsOffset[1] = vOffset[1];    
            hsCount[0] = vCount[0]; hsCount[1] = vCount[1];
            
            
            if(vOffset[0] + hsCount[0] <= dimDataset[0] && vOffset[1] + hsCount[1] <= dimDataset[1] && 
               DatasetValues.size()<= dimDatasetinFile[0] * dimDatasetinFile[1]) {
                
                H5::DataSpace dataspace(RANK2, hsCount);
                H5::DataSpace memspace(RANK2, hsCount, NULL);
                
                dataspace = pdataset->getSpace();
                dataspace.selectHyperslab( H5S_SELECT_SET, hsCount, hsOffset, hsStride, hsBlock);
                
                pdataset->write(&DatasetValues[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
                memspace.close();
                dataspace.close();
                
            } else {
                
                Rcpp::Rcerr<<"\nIt is not possible to write block in current position (writeDatasetBlock)";
            }
                
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock std::vector (File IException)";
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock std::vector (DataSet IException)";
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock std::vector (Group IException)";
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock std::vector (DataSpace IException)";
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            close_dataset();
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception writeDatasetBlock std::vector (Data TypeIException)";
        }
        return void();
    }
    

    
    // Read rhdf5 data matrix subset, 
    // input : 
    //      ivoffset : start position
    //      ivcount : block size
    //      ivstride :(1,1) by default.
    //      ivblock : (1,1) by default.
    // output : 
    //    rdatablock : matrix block
    void readDatasetBlock(std::vector<hsize_t> ivoffset, std::vector<hsize_t> ivcount,
                             std::vector<hsize_t> ivstride, std::vector<hsize_t> ivblock,
                             double* rdatablock)
    {
        
        try
        {
            H5::Exception::dontPrint();
            
            hsize_t offset[2], count[2], stride[2], block[2];
            
            offset[0] = ivoffset[0]; offset[1] = ivoffset[1];
            count[0] = ivcount[0]; count[1] = ivcount[1];
            stride[0] = ivstride[0]; stride[1] = ivstride[1];
            block[0] = ivblock[0]; block[1] = ivblock[1];
            
            // Define Memory Dataspace. Get file dataspace and select a subset from the file dataspace.
            hsize_t dimsm[2];
            dimsm[0] = count[0]; 
            dimsm[1] = count[1];
            
            H5::DataSpace memspace(RANK2, dimsm, NULL);
            
            //  Get dataset dataspace.
            H5::DataSpace dataspace = pdataset->getSpace();
            dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
            
            H5T_class_t type_class = pdataset->getTypeClass();
            
            // Get class of datatype and print message if it's an integer.
            if( type_class == H5T_INTEGER || type_class == H5T_FLOAT ) {
                pdataset->read( rdatablock, H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
            } else {
                Rcpp::Rcerr<<"\nc++ exception readDatasetBlock (Data type not allowed, maybe are trying to read string matrix?)";
                return void();
                // return(nullptr);
            }// else if (type_class == H5T_FLOAT) {
            //     pdataset->read( rdatablock, H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
            // } 

            memspace.close();
            dataspace.close();
            
        } catch( H5::FileIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception readDatasetBlock (File IException)";
            // return(nullptr);
            return void();
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception readDatasetBlock (DataSet IException)";
            // return(nullptr);
            return void();
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception readDatasetBlock (Group IException)";
            // return(nullptr);
            return void();
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception readDatasetBlock (DataSpace IException)";
            // return(nullptr);
            return void();
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception readDatasetBlock (Data TypeIException)";
            // return(nullptr);
            return void();
        } catch(std::exception &ex) {
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception readDatasetBlock \n"<< ex.what();
            return void();
        } catch (...) {
            close_dataset_file();
            Rcpp::Rcerr<<"\nC++ exception readDatasetBlock (unknown reason)";
            return void();
        }
        // return(rdatablock);
        return void();
    }

    
    
    
    /* Create single attribute to hdf5 group */
    /**
     * @brief Add an attribute to the dataset
     * @details Adds or updates an attribute associated with the dataset.
     * Supports various R data types.
     * 
     * @param attrName Name of the attribute
     * @param attr_data Attribute data
     * @return EXEC_OK on success, EXEC_ERROR on failure
     */
    int addAttribute( std::string attrName, Rcpp::RObject attr_data)
    {
        
        hsize_t dims[1] = { DIM1 };
        
        try
        {
            H5::Exception::dontPrint();
            
            // // Open an existing file and dataset.
            // H5File file( HiCfilename, H5F_ACC_RDWR );
            
            // Create the data space for the attribute.
            H5::DataSpace attr_dataspace = H5::DataSpace (1, dims );
            
            // Group group = file.openGroup( HiCObject);
            
            if( Rf_isString(attr_data) )
            {
                // Prepare string data
                int strlength = Rcpp::as<std::string>(attr_data).size();
                //..// char stringdata[strlength+1];
                //..// 
                //..// std::string word = Rcpp::as<std::string>(attr_data).c_str();
                //..// 
                //..// int j = 0;
                //..// for( j=0; (unsigned)j < word.size() && j < (strlength); j++ )
                //..//     stringdata[j] = word[j];
                //..// 
                //..// stringdata[j] = '\0'; // insert hdf5 end of string
                
                // Create group attribute and write
                H5::Attribute attribute =  pdataset->createAttribute(attrName, 
                                                            H5::StrType(H5::PredType::C_S1, strlength ), 
                                                            attr_dataspace);
                
                //..// attribute.write(H5::StrType(H5::PredType::C_S1, strlength), stringdata);
                attribute.write(H5::StrType(H5::PredType::C_S1, strlength), Rcpp::as<std::string>(attr_data).data());
            } 
            else if( Rf_isInteger(attr_data) ) 
            {
                int data[] =  {Rcpp::as<int>(attr_data)};
                // Create group attribute and write
                H5::Attribute attribute = pdataset->createAttribute(attrName, H5::PredType::NATIVE_INT,attr_dataspace);
                attribute.write(H5::PredType::NATIVE_INT, data);
            }
            else if( Rf_isNumeric(attr_data) ) 
            {
                double data[] = {Rcpp::as<double>(attr_data)};
                // Create group attribute and write
                H5::Attribute attribute = pdataset->createAttribute(attrName, H5::PredType::NATIVE_DOUBLE,attr_dataspace);
                attribute.write(H5::PredType::NATIVE_DOUBLE, data);
            }
            
            attr_dataspace.close();
        }  // end of try block
        
        catch( H5::DataSpaceIException& error ) {
            error.printErrorStack();
            return -1;
        } catch( H5::AttributeIException& error ) { // catch failure caused by the H5File operations
            error.printErrorStack();
            return -1;
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            error.printErrorStack();
            return -1;
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            error.printErrorStack();
            return -1;
        }
        
        return 0;  // successfully terminated
        
    }
    
    
    /**
     * @brief Get an attribute from the dataset
     * @details Retrieves the value of a named attribute from the dataset.
     * 
     * @param strAtribute Name of the attribute to retrieve
     */
    void getAttribute(std::string strAtribute)
    {
        try
        {
            // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
            H5::Exception::dontPrint();
            
            if(  strAtribute == "internal" && !exists_HDF5_element(pfile, strAtribute)) {
                internalDataset = false;
            } else {
                
                H5::Attribute *attr = new H5::Attribute(pdataset->openAttribute(strAtribute));
                H5::DataType  *type = new H5::DataType(attr->getDataType());
                H5::StrType stype = attr->getStrType();
                
                if( type->getClass() == H5T_INTEGER) {
                    Rcpp::Rcout<<"\nInteger\n";
                } else if (type->getClass() == H5T_FLOAT) {
                    Rcpp::Rcout<<"\nDouble\n";
                } else if (type->getClass() == H5T_STRING) {
                    
                    std::string strAttrVal;
                    attr->read(stype, strAttrVal);
                    
                    if( strAtribute == "internal" ) {
                        if( strAttrVal == "0" ) {
                            internalDataset = false;
                        } else {
                            internalDataset = true; }
                    }
                } else {
                    Rcpp::Rcout<<"\nPos No SEP"<<type->getClass()<<"\n";
                }
                
                attr->close();
                
            }
            
        } catch( H5::FileIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception getAttribute (File IException)";
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception getAttribute (DataSet IException)";
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception getAttribute (Group IException)";
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception getAttribute (DataSpace IException)";
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            Rcpp::Rcerr<<"\nc++ exception getAttribute (Data TypeIException)";
        }
        return void();
    }
    
    /**
     * @brief Get dataset pointer
     * @return Pointer of the dataset
     */
    H5::DataSet* getDatasetptr() { return(pdataset); }  // Return dataset pointer
    
    /**
     * @brief Get dataset name
     * @return Name of the dataset
     */
    std::string getDatasetName() { return(name); }  // Return dataset name
    
    /**
     * @brief Get group name
     * @return Name of the group containing the dataset
     */
    std::string getGroup() { return(getGroupName()); }  // Return group name
    
    /**
     * @brief Get file name
     * @return Name of the file containing the dataset
     */
    std::string getFileName() { return(getFilename()); }  // Return file name
    
    /**
     * @brief Get number of rows in internal format
     * @return Number of rows in the dataset's internal storage
     */
    hsize_t nrows() { return(dimDataset[0]); }  // Return number of rows
    
    /**
     * @brief Get number of columns in internal format
     * @return Number of columns in the dataset's internal storage
     */
    hsize_t ncols() { return(dimDataset[1]); }  // Return number of columns
    
    /**
     * @brief Get number of rows in R format
     * @return Number of rows in R's representation (transposed)
     */
    hsize_t nrows_r() { return(dimDatasetinFile[1]); }  // Return number of rows in R (transposed values)
    
    /**
     * @brief Get number of columns in R format
     * @return Number of columns in R's representation (transposed)
     */
    hsize_t ncols_r() { return(dimDatasetinFile[0]); }  // Return number of columns in R (transposed values)
    
    /**
     * @brief Get dataset dimension
     * @return Pointer to dataset dimension (rows x columns)
     */
    hsize_t* dim() { return(dimDataset); }  // Return dataset dimension (rows x columns)
    
    /**
     * @brief Get dataset file dimension
     * @return Pointer to dataset file dimension (rows x columns)
     */
    hsize_t* dimFile() { return(dimDatasetinFile); } // Return dataset file dimensions (rows x columns)
    
    /**
     * @brief Get number of rows in file
     * @return Number of rows in file storage
     */
    hsize_t nrows_file() { return(dimDatasetinFile[0]); }  // Return number of rows
    
    /**
     * @brief Get number of columns in file
     * @return Number of columns in file storage
     */
    hsize_t ncols_file() { return(dimDatasetinFile[1]); }  // Return number of columns
    
    
    /**
     * @brief Check if dataset is unlimited
     * @return True if dataset has unlimited dimensions
     */
    bool isUnlimited() { return(unlimited); }  // Return if dataset is an unlimited dataset
    
    /**
     * @brief Check if dataset is internal
     * @return True if dataset is used for internal storage
     */
    bool isInternal() { return(internalDataset); }  // Return if dataset is an unlimited dataset
    
    /**
     * @brief Check if dataset is open
     * @return True if dataset is currently open
     */
    bool isOpen() { return( isDatasetOpen()); } // Return if dataset is open and exists or not
    
    
    // Destructor
    
    /**
     * @brief Destructor
     * @details Closes the dataset and releases resources
     */
    virtual ~hdf5Dataset(){
        if(pdataset){
            pdataset->close();
        }
    }
    
    
protected:
    
    // ------------------------
    //   Struct declaration
    // -----------------------
    
    /**
     * @brief Structure for string data
     * @details Used for storing fixed-length string data in HDF5
     */
    typedef struct names {
        char chr[MAXSTRING];
    } names;
    
    
    // ------------------------
    //   Variables declaration
    // ------------------------
    
     /**
     * @brief Dataset pointer
     * @details Pointer to the HDF5 dataset object
     */
    H5::DataSet* pdataset = nullptr;

    /**
     * @brief Dataset dimensions
     * @details Array storing the dimensions of the dataset [rows, columns]
     */
    hsize_t dimDataset[2];

    /**
     * @brief File storage dimensions
     * @details Array storing the dimensions as stored in file [rows, columns]
     */
    hsize_t dimDatasetinFile[2];

    /**
     * @brief Dataset name
     */
    std::string name;

    /**
     * @brief Dataset type
     * @details String indicating the data type ("int", "numeric", "string", etc.)
     */
    std::string type;

    /**
     * @brief Column names dataset
     * @details Name of associated dataset storing column names
     */
    std::string colnamesDataset;

    /**
     * @brief Row names dataset
     * @details Name of associated dataset storing row names
     */
    std::string rownamesDataset;

    /**
     * @brief Overwrite flag
     * @details Whether to overwrite existing datasets
     */
    bool boverwrite;

    /**
     * @brief Unlimited dimensions flag
     * @details Whether the dataset has unlimited dimensions
     */
    bool unlimited = false;

    /**
     * @brief Internal dataset flag
     * @details Whether the dataset is used for internal storage
     */
    bool internalDataset = false;
    
    
    // ------------------------
    //   Function declarations
    // ------------------------
    
    /**
     * @brief Close the dataset
     * @details Closes the dataset and releases associated resources
     */
    void close_dataset() 
    {
        pdataset->close();
    }
    
    /**
     * @brief Close dataset and file
     * @details Closes both the dataset and its associated file
     */
    void close_dataset_file()
    {
        pdataset->close();
        pfile->close();
    }
    
    /**
     * @brief Convert DataFrame to range list
     * @details Converts R DataFrame to HDF5-compatible range list format
     * 
     * @param DatasetValues R DataFrame to convert
     * @param bFullDataset Whether to convert entire dataset
     * @return Array of name structures
     */
    names* convert_DataFrame_to_RangeList(Rcpp::RObject DatasetValues, bool bFullDataset) 
    {
        
        int datarows = Rcpp::as<Rcpp::StringVector>(DatasetValues).size();
        int isizetoWrite = datarows;
        
        if(bFullDataset) {
            isizetoWrite = dimDataset[0] * dimDataset[1]; } 
        
        names *names_list = new names[isizetoWrite];  // Convert to range list
        
        // Write data to dataset, if data to write is smaller than dataset then empty positions='\0' 
        for( int i = 0; i < isizetoWrite; i++ ) {
            int j = 0;
            if(i< datarows) {
                Rcpp::String wchrom = Rcpp::as<Rcpp::StringVector>(DatasetValues)(i);
                std::string word = wchrom.get_cstring();
                
                for( j = 0; (unsigned)j < word.size() && j < (MAXSTRING-1); j++ ) {
                    names_list[i].chr[j] = word[j]; 
                }
            }
            names_list[i].chr[j] = '\0'; // insert hdf5 end of string
        }
        return(names_list);
    }
    
    
    /**
     * @brief Get dimensions of existing dataset
     * @details Retrieves and stores the dimensions of an existing dataset
     */
    void getDimensExistingDataset()
    {
        try
        {
            H5::Exception::dontPrint();
            // Get dataspace from dataset
            H5::DataSpace dataspace = pdataset->getSpace();
            
            // Get the number of dimensions in the dataspace.
            int rank = dataspace.getSimpleExtentNdims();
            
            // Get the dimension size of each dimension in the dataspace
            hsize_t dims_out[2];
            hsize_t dims_max[2];
            int ndims = dataspace.getSimpleExtentDims( dims_out, dims_max);
            
            if( ndims>0 ) {
                
                if(dims_max[0] == H5S_UNLIMITED) {
                    unlimited = true;
                }
                
                if( rank == 1) {
                    // dims = IntegerVector::create( static_cast<int>(dims_out[0]), static_cast<int>(1));
                    dimDataset[0] = dims_out[0];
                    dimDataset[1] = 1;
                } else if( rank == 2 ){
                    // dims = IntegerVector::create(static_cast<int>(dims_out[0]), static_cast<int>(dims_out[1]));
                    dimDataset[0] = dims_out[0];
                    dimDataset[1] = dims_out[1];
                }
                
                dimDatasetinFile[0] = dims_out[0];
                dimDatasetinFile[1] = dims_out[1];    
            }
            
        } catch( H5::FileIException& error) { 
            Rcpp::Rcerr<<"\nc++ exception getDimensExistingDataset (File IException)";
        } catch(H5::DataSetIException& error) { 
            Rcpp::Rcerr<<"\nc++ exception getDimensExistingDataset (DataSet IException)";
        } catch(H5::GroupIException& error) { 
            Rcpp::Rcerr<<"\nc++ exception getDimensExistingDataset (Group IException)";
        } catch(H5::DataSpaceIException& error) { 
            Rcpp::Rcerr<<"\nc++ exception getDimensExistingDataset (DataSpace IException)";
        } 
        
        return void();
    }
    
    
    /**
     * @brief Check if dataset is open
     * @details Verifies if the dataset is currently open and accessible
     * 
     * @return True if dataset is open, false otherwise
     */
    bool isDatasetOpen() {
        hid_t datasetId = pdataset->getId();
        // Check ID
        if (H5Iis_valid(datasetId) > 0) {
            return true; 
        } else {
            return false; 
        }
    }
    
};
}

#endif // BIGDATASTATMETH_HDF5_DATASETS_HPP