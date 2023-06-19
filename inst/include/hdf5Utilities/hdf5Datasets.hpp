#ifndef BIGDATASTATMETH_HDF5_DATASETS_HPP
#define BIGDATASTATMETH_HDF5_DATASETS_HPP

#include "BigDataStatMeth.hpp"


namespace BigDataStatMeth {

class hdf5Dataset : public hdf5Group 
{
    
public:

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
        // pdataset = openDataset();
        // getDimensExistingDataset();
        
    }
    
    
    hdf5Dataset(std::string filename, std::string group, std::string datasetname, bool overwrite) : 
    hdf5Group(filename, group)
    {
        name = datasetname;
        boverwrite = overwrite;
        // pdataset = openDataset();
        // getDimensExistingDataset();
    }
    
    
    
    
    // Create empty hdf5 DataSet, strdatatype can be: 
    //  . "int": integer dataset
    //  . "numeric" or "real": double dataset
    //  . "string": string dataset
    void createDataset(size_t rows, size_t cols, std::string strdatatype) 
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
                ::Rf_error("Dataset exits, please set overwrite = true to overwrite the existing dataset (DataSet IException)");
            } else {
                
                if( boverwrite == true && bexists == true) {
                    remove_elements(pfile, getGroupName(), {name}); 
                    bRemoved = true;
                }
                
                type = strdatatype;
                
                if( type == "string") {
                    // Create the memory datatype.
                    H5::CompType strtype(sizeof(names));
                    strtype.insertMember("chr", HOFFSET(names, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
                    pdataset = new H5::DataSet(pfile->createDataSet(fullDatasetPath, strtype, dataspace));
                } else if( type == "int" || type == "logic" || type == "factor") {
                    H5::IntType datatype( H5::PredType::NATIVE_INT );
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace ));
                    if(bRemoved == true) {
                        writeDataset(Rcpp::wrap(Eigen::MatrixXd::Zero(dimDataset[0], dimDataset[1]) ));    
                    }
                } else if( type == "numeric" || type == "real") {
                    H5::IntType datatype( H5::PredType::NATIVE_DOUBLE ); 
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace ));
                    if(bRemoved == true) {
                        writeDataset(Rcpp::wrap(Eigen::MatrixXd::Zero(dimDataset[0], dimDataset[1]) ));
                    }
                } else {
                    ::Rf_error( "Dataset data type not allowed or no matrix defined (createDataset)" );
                }
            }
            
            dataspace.close();
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createDataset (File IException) " );
        } catch(H5::GroupIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createDataset (Group IException)");
        } catch(H5::DataSetIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createDataset (DataSet IException)");
        } 
        return void();
    }
    
    
    void createDataset(BigDataStatMeth::hdf5Dataset* dsLike, std::string strdatatype) 
    {
        try{
            
            createDataset( dsLike->ncols(), dsLike->nrows(), strdatatype);
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createDataset (File IException) " );
        } catch(H5::GroupIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createDataset (Group IException)");
        } catch(H5::DataSetIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createDataset (DataSet IException)");
        } 
        
        return void();
    }
    
    // Create empty hdf5 data file
    void createUnlimitedDataset(size_t rows, size_t cols, std::string strdatatype) 
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
            
            if( !exists_HDF5_element(pfile, groupname) ) {
                create_HDF5_groups(groupname);
            }
            
            if( exists_HDF5_element(pfile, fullDatasetPath) && boverwrite == false) {
                Rcpp::Rcout<<"\n Dataset exits, please set overwrite = true to overwrite data \n";
            } else {
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
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createUnlimitedDataset (File IException) " );
        } catch(H5::GroupIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createUnlimitedDataset (File IException) " );
        } catch(H5::DataSetIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception createUnlimitedDataset (File IException) " );
        } 
        return void();
    }
    
    
    void extendUnlimitedDataset(const size_t rows, const size_t cols)
    {
        try
        {
            if(unlimited == true) {
                H5::Exception::dontPrint();
                
                // Extend dataset size to:  oldDims + newDims
                hsize_t newdims[2];
                newdims[0] = cols;
                newdims[1] = rows;
                
                hsize_t size[2];
                dimDataset[0] = dimDataset[0] + newdims[0];
                dimDataset[1] = dimDataset[1] + newdims[1];
                
                pdataset->extend( dimDataset );    
            } else {
                Rcpp::Rcout<<"\n Dataset is not an unlimited dataset, fixed datasets can't be extended\n";
                return void();
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (File IException)" );
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            ::Rf_error( "c++ exception extend_HDF5_matrix_subset_ptr (DataSet IException)" );
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
                ::Rf_error( "c++ exception Please create file before proceed" );
            } else { 
                bool bexists = exists_HDF5_element(pfile, fullPath);
                if( bexists ) {
                    if ((pdataset == NULL) == TRUE) {
                        pdataset = new H5::DataSet(pfile->openDataSet(fullPath));    
                    }
                    getDimensExistingDataset();
                } else {
                    ::Rf_error( "c++ exception Please create Dataset before proceed" );
                }
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception hdf5File (File IException) " );
        } catch(H5::GroupIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception hdf5Dataset (File IException) " );
        } catch(H5::DataSetIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception hdf5Dataset (File IException) " );
        } 
        return(pdataset);
    }
    
    
    
    // Write data in existing hdf5 dataSet - Only writes data if the size of existing hdf5 dataset
    // is equal to input data (size)
    void writeDataset( Rcpp::RObject DatasetValues )
    {
        
        try
        {
            H5::Exception::dontPrint();
            std::vector<int> dims;
            
            if(Rcpp::is<Rcpp::NumericMatrix>(DatasetValues) || Rcpp::is<Rcpp::IntegerMatrix>(DatasetValues) ) {
                hsize_t ncolRObject = Rcpp::as<Rcpp::NumericMatrix>(DatasetValues).ncol();
                hsize_t nrowRObject = Rcpp::as<Rcpp::NumericMatrix>(DatasetValues).nrow();
                
                if( (nrowRObject == dimDataset[1] && ncolRObject == dimDataset[0]) || ncolRObject == 1 || nrowRObject == 1 ) {
                    H5::DataSpace dataspace(RANK2, dimDataset);
                    
                    std::vector<double> matHiCValues = Rcpp::as<std::vector<double> >(Rcpp::as<Rcpp::NumericMatrix>(DatasetValues));
                    pdataset->write( &matHiCValues[0] , H5::PredType::NATIVE_DOUBLE);
                    dataspace.close();
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
                ::Rf_error( "Matrix data type not allowed (writeDataset)" );
            }
            
        } catch(H5::FileIException error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception writeDataset (File IException)" );
        } catch(H5::DataSetIException error) { // catch failure caused by the DataSet operations
            ::Rf_error( "c++ exception writeDataset (DataSet IException)" );
        } catch(H5::GroupIException error) { // catch failure caused by the Group operations
            ::Rf_error( "c++ exception writeDataset (Group IException)" );
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception writeDataset (DataSpace IException)" );
        } catch(H5::DataTypeIException error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception writeDataset (Data TypeIException)" );
        }
        return void();
    }
    
    
    
    
    // SPECIFIC FUNCTIONS TO WORK WITH EIGEN OBJECTS -- RowMajor Matrix
    // Write data in existing hdf5 dataSet - Writes a block inside a hdf5
    // dataset starting at offset position x-y
    void writeRowMajorDatasetBlock( Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> DatasetValues, 
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
                ::Rf_error( "It is not possible to write block in current position (writeRowMajorDatasetBlock)" );
            }
                
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception writeRowMajorDatasetBlock (File IException)" );
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            ::Rf_error( "c++ exception writeRowMajorDatasetBlock (DataSet IException)" );
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            ::Rf_error( "c++ exception writeRowMajorDatasetBlock (Group IException)" );
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception writeRowMajorDatasetBlock (DataSpace IException)" );
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception writeRowMajorDatasetBlock (Data TypeIException)" );
        }
        return void();
    }
    
    
    // GENERIC FUNCTION USE WITH R-OBJECTS 
    // Write data in existing hdf5 dataSet - Writes a block inside a hdf5
    // dataset starting at offset position x-y
    void writeDatasetBlock( Rcpp::RObject DatasetValues, std::vector<hsize_t> vOffset, 
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
                    ::Rf_error( "It is not possible to write block in current position (writeDatasetBlock)" );
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
                    ::Rf_error( "It is not possible to write block in current position (writeDatasetBlock)" );
                }
            } else {
                ::Rf_error( "Matrix data type not allowed (writeDatasetBlock)" );
            }

        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception writeDatasetBlock (File IException)" );
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            ::Rf_error( "c++ exception writeDatasetBlock (DataSet IException)" );
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            ::Rf_error( "c++ exception writeDatasetBlock (Group IException)" );
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception writeDatasetBlock (DataSpace IException)" );
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception writeDatasetBlock (Data TypeIException)" );
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
    double* readDatasetBlock(std::vector<hsize_t> ivoffset, std::vector<hsize_t> ivcount,
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
            
            //  Get dataspace of the dataset.
            H5::DataSpace dataspace = pdataset->getSpace();
            dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
            
            H5T_class_t type_class = pdataset->getTypeClass();
            
            // Get class of datatype and print message if it's an integer.
            if( type_class == H5T_INTEGER || type_class == H5T_FLOAT ) {
                pdataset->read( rdatablock, H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
            } else {
                ::Rf_error( "c++ exception readDatasetBlock (Data type not allowed, maybe are trying to read string matrix?)" );
                return(nullptr);
            }// else if (type_class == H5T_FLOAT) {
            //     pdataset->read( rdatablock, H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
            // } 

            memspace.close();
            dataspace.close();
            
        } catch( H5::FileIException& error) { 
            close_dataset_file();
            ::Rf_error( "c++ exception readDatasetBlock (File IException)" );
            return(nullptr);
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            ::Rf_error( "c++ exception readDatasetBlock (DataSet IException)" );
            return(nullptr);
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            ::Rf_error( "c++ exception readDatasetBlock (Group IException)" );
            return(nullptr);
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            ::Rf_error( "c++ exception readDatasetBlock (DataSpace IException)" );
            return(nullptr);
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            ::Rf_error( "c++ exception readDatasetBlock (Data TypeIException)" );
            return(nullptr);
        }
        return(rdatablock);
    }
    
    
    
    H5::DataSet* getDatasetptr() { return(pdataset); }  // Return dataset pointer
    std::string getDatasetName() { return(name); }  // Return dataset name
    std::string getGroup() { return(getGroupName()); }  // Return group name
    std::string getFile() { return(getFilename()); }  // Return group name
    hsize_t nrows() { return(dimDataset[0]); }  // Return number of rows
    hsize_t ncols() { return(dimDataset[1]); }  // Return number of columns
    hsize_t* dim() { return(dimDataset); }  // Return dataset dimension (rows x columns)
    
    // Destructor
    ~hdf5Dataset(){
        if(pdataset){
            pdataset->close();
        }
    }
    
    
private:
    
    // ------------------------
    //   Struct declaration
    // ------------------------
    
    typedef struct names {
        char chr[MAXSTRING];
    } names;
    
    
    // ------------------------
    //   Variables declaration
    // ------------------------
    
    std::string name;
    std::string type;
    bool boverwrite;
    H5::DataSet* pdataset = nullptr;
    hsize_t dimDataset[2];
    hsize_t dimDatasetinFile[2];
    bool unlimited = false;
    
    
    // ------------------------
    //   Function declarations
    // ------------------------
    
    void close_dataset() 
    {
        pdataset->close();
    }
    
    void close_dataset_file()
    {
        pdataset->close();
        pfile->close();
    }

    
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
                
                for( j = 0; j < word.size() && j < (MAXSTRING-1); j++ ) {
                    names_list[i].chr[j] = word[j]; 
                }
            }
            names_list[i].chr[j] = '\0'; // insert hdf5 end of string
        }
        return(names_list);
    }
    
    
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
            int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
            
            if( rank == 1) {
                // dims = IntegerVector::create( static_cast<int>(dims_out[0]), static_cast<int>(1));
                dimDataset[0] = dims_out[0];
                dimDataset[1] = 1;
            } else if( rank == 2 ){
                // dims = IntegerVector::create(static_cast<int>(dims_out[0]), static_cast<int>(dims_out[1]));
                dimDataset[0] = dims_out[0];
                dimDataset[1] = dims_out[1];
            }
        } catch( H5::FileIException& error) { 
            ::Rf_error( "c++ exception getDimensExistingDataset (File IException)" );
        } catch(H5::DataSetIException& error) { 
            ::Rf_error( "c++ exception getDimensExistingDataset (DataSet IException)" );
        } catch(H5::GroupIException& error) { 
            ::Rf_error( "c++ exception getDimensExistingDataset (Group IException)" );
        } catch(H5::DataSpaceIException& error) { 
            ::Rf_error( "c++ exception getDimensExistingDataset (DataSpace IException)" );
        } 
        
        return void();
    }  
    
    
    
};
}

#endif // BIGDATASTATMETH_HDF5_DATASETS_HPP