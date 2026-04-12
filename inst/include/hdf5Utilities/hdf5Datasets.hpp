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
 * 
 * Compression level management:
 * - Each dataset instance tracks its own compression level (default = 6, gzip balanced)
 * - When creating a new dataset, compression_level controls gzip compression (0 = none, 1-9)
 * - When opening an existing dataset via openDataset(), the actual compression level is
 *   read from the HDF5 dataset creation property list (dcpl) and stored in compression_level.
 *   If no deflate filter is present, compression_level is set to 0.
 * - Use setCompressionLevel() to override before createDataset().
 * - Use getCompressionLevel() to query the current level (useful when propagating
 *   compression from an input dataset to output datasets).
 */

#ifndef BIGDATASTATMETH_HDF5_DATASETS_HPP
#define BIGDATASTATMETH_HDF5_DATASETS_HPP

// #include "Utilities/Utilities.hpp"
// #include <hdf5.h>

#include <RcppEigen.h>
#include "H5Cpp.h"


namespace BigDataStatMeth {

#ifdef _OPENMP
    class HDF5ThreadSafety {
    private:
        static omp_lock_t& getHDF5Lock() {
            static omp_lock_t hdf5_lock;
            return hdf5_lock;
        }
        
        static bool& getLockInitialized() {
            static bool lock_initialized = false;
            return lock_initialized;
        }
        
    public:
        static void initLock() {
            if (!getLockInitialized()) {
                omp_init_lock(&getHDF5Lock());
                getLockInitialized() = true;
            }
        }
        
        // RAII Lock Guard - automatically unlocks on scope exit or exception
        class LockGuard {
        public:
            LockGuard() {
                initLock();
                omp_set_lock(&getHDF5Lock());
            }
            
            ~LockGuard() {
                omp_unset_lock(&getHDF5Lock());
            }
            
            // Non-copyable
            LockGuard(const LockGuard&) = delete;
            LockGuard& operator=(const LockGuard&) = delete;
        };
    };
#endif


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
     * @brief Set the compression level for this dataset instance
     * @details Sets the compression level used by createDataset() when no explicit
     * compression_level argument is passed. Call this before createDataset() to
     * propagate compression settings from an input dataset to an output dataset.
     * @param level Compression level 0-9 (0 = no compression, 6 = balanced default)
     */
    void setCompressionLevel(int level) {
        compression_level = level;
    }
    
    
    /**
     * @brief Inherit the compression level from an input dataset if not already set
     * @details Sets the compression level only when no explicit value has been assigned
     * (i.e., compression_level < 0, the initial sentinel value). This allows header-level
     * algorithms to propagate compression from input to output datasets while respecting
     * any prior call to setCompressionLevel() made by a higher-level wrapper (R6 method,
     * C++ API caller, or global option).
     *
     * Priority order respected by this design:
     *   1. Explicit user value via setCompressionLevel()
     *   2. Global option passed by the R6 wrapper via setCompressionLevel()
     *   3. Inherited from input dataset via inheritCompressionLevel()  ← this function
     *   4. Built-in default (6) applied in createDataset() when level is still -1
     *
     * Use this instead of setCompressionLevel() inside algorithm headers:
     * @code
     *   dsOut->inheritCompressionLevel(dsIn->getCompressionLevel());
     *   dsOut->createDataset(rows, cols, "numeric");
     * @endcode
     * @param level Compression level 0-9 read from an input dataset via getCompressionLevel()
     */
    void inheritCompressionLevel(int level) {
        if (compression_level < 0)
            compression_level = level;
    }
    
    /**
     * @brief Get the current compression level for this dataset instance
     * @details Returns the compression level that will be used by createDataset()
     * when no explicit argument is passed. After openDataset() is called, this
     * reflects the actual compression level stored in the HDF5 file (0 if the
     * dataset was created without compression).
     * Use this to propagate compression from an input dataset to output datasets:
     * @code
     *   dsOut->setCompressionLevel(dsIn->getCompressionLevel());
     *   dsOut->createDataset(rows, cols, "numeric");
     * @endcode
     * @return Compression level 0-9
     */
    int getCompressionLevel() const {
        return compression_level;
    }
    
    
    /**
    * @brief Create a new dataset with optional compression
    * @details Creates a new HDF5 dataset with specified dimensions, data type, and compression.
    * Supported data types:
    * - "int": Integer dataset
    * - "numeric" or "real": Double dataset
    * - "string": String dataset
    * 
    * Compression features:
    * - Automatic chunking when compression is enabled
    * - gzip compression with configurable level
    * - Shuffle filter for improved compression ratio
    * - Completely transparent - no changes needed for read/write operations
    * 
    * @param rows Number of rows
    * @param cols Number of columns  
    * @param strdatatype Data type ("int", "numeric", "real", or "string")
    * @param compression_level Compression level (0=none, 1-9=gzip level, 6=balanced default)
    *                         Higher values = better compression, more CPU usage
    * 
    * @throws H5::FileIException if file operations fail
    * @throws H5::GroupIException if group operations fail
    * @throws H5::DataSetIException if dataset operations fail
    * 
    * @note Compression is applied only at dataset creation time
    * @note Chunking is automatically configured when compression > 0
    * @note All subsequent read/write operations are transparent regardless of compression
    * @note Typical space savings: 60-80% with compression_level=6
    * 
    * @since Added compression support in version X.X.X
    */
    virtual void createDataset(size_t rows, size_t cols, std::string strdatatype, int compression_level = -1) 
    {
        try
        {
            
            H5::Exception::dontPrint();
            std::string fullDatasetPath = groupname + "/" + name;
            bool bRemoved = false;
            
            int eff_compression = (compression_level < 0) ? this->compression_level : compression_level;
            if (eff_compression < 0) eff_compression = 6;  // default if eff_compression not setted
            
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
                Rf_error("Dataset already exits, please set overwrite = true to overwrite the existing dataset (DataSet IException)");
                return void();
            } else {
                
                if( boverwrite == true && bexists == true) {
                    remove_elements(pfile, getGroupName(), {name}); 
                    bRemoved = true;
                }
                
                type = strdatatype;
                
                // Configure compression and chunking
                hid_t cparms = H5P_DEFAULT;
                // if (compression_level > 0) {
                if (eff_compression > 0) {                    
                    cparms = H5Pcreate(H5P_DATASET_CREATE);
                    
                    // Configure chunking (required for compression)
                    hsize_t chunk_dims[2] = {
                        std::min((hsize_t)1024, dimDataset[0]),
                        std::min((hsize_t)1024, dimDataset[1])
                    };
                    H5Pset_chunk(cparms, RANK2, chunk_dims);
                    
                    // Configure compression
                    // H5Pset_deflate(cparms, compression_level);  // gzip compression
                    H5Pset_deflate(cparms, eff_compression);      // gzip compression
                    H5Pset_shuffle(cparms);                     // shuffle filter for better compression
                }
                
                if( type == "string") {
                    H5::CompType strtype(sizeof(names));
                    strtype.insertMember("chr", HOFFSET(names, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
                    pdataset = new H5::DataSet(pfile->createDataSet(fullDatasetPath, strtype, dataspace, cparms));
                } else if( type == "int" || type == "logic" || type == "factor") {
                    H5::IntType datatype( H5::PredType::NATIVE_INT );
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace, cparms ));
                    if(bRemoved == true) {
                        writeDataset(Rcpp::wrap(Eigen::MatrixXd::Zero(dimDatasetinFile[0], dimDatasetinFile[1]) ));
                    }
                } else if( type == "numeric" || type == "real" || type == "dataframe") {
                    H5::FloatType datatype( H5::PredType::NATIVE_DOUBLE ); 
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace, cparms ));
                    if(bRemoved == true) {
                        writeDataset(Rcpp::wrap(Eigen::MatrixXd::Zero(dimDatasetinFile[0], dimDatasetinFile[1]) ));
                    }
                } else {
                    if (cparms != H5P_DEFAULT) H5Pclose(cparms);
                    close_file();
                    Rf_error("Dataset data type not allowed or no matrix defined (createDataset)");
                }
                
                // Clean up compression properties
                if (cparms != H5P_DEFAULT) H5Pclose(cparms);
            }
            
            dataspace.close();
            
            addAttribute( "internal", Rcpp::wrap("0") );
            
        } catch(H5::FileIException& error) { 
            close_file();
            Rf_error("c++ exception createDataset (File IException)");
        } catch(H5::GroupIException& error) { 
            close_file();
            Rf_error("c++ exception createDataset (Group IException)");
        } catch(H5::DataSetIException& error) { 
            close_file();
            Rf_error("c++ exception createDataset (DataSet IException)");
        } 
        return void();
    }
    
 
    /**
     * @brief Create a dataset based on another dataset's dimensions with optional compression
     * @details Creates a new dataset with the same dimensions as the reference dataset
     * but with a specified data type and compression settings. This is useful for creating
     * derived datasets or transformations while maintaining dimensional consistency.
     * 
     * Features:
     * - Inherits dimensions from reference dataset
     * - Independent data type specification
     * - Same compression capabilities as primary createDataset()
     * - Maintains dimensional relationships for related datasets
     * 
     * @param dsLike Reference dataset to copy dimensions from
     * @param strdatatype Data type for the new dataset ("int", "numeric", "real", "string")
     * @param compression_level Compression level (0=none, 1-9=gzip level, 6=balanced default)
     *                         0 = No compression
     *                         1-3 = Light compression (fast)
     *                         4-6 = Balanced compression (recommended)
     *                         7-9 = Maximum compression (slower)
     * 
     * @throws H5::FileIException if file operations fail
     * @throws H5::GroupIException if group operations fail  
     * @throws H5::DataSetIException if dataset operations fail
     * 
     * @note Reference dataset must be properly initialized with valid dimensions
     * @note New dataset is independent - changes don't affect reference dataset
     * @note Compression settings are applied to new dataset only
     * 
     * @see createDataset(size_t, size_t, std::string, int) for detailed compression documentation
     * @since Added compression support in version X.X.X
     */
    virtual void createDataset(BigDataStatMeth::hdf5Dataset* dsLike, std::string strdatatype, int compression_level = -1) 
    {
        try{
            createDataset( dsLike->ncols(), dsLike->nrows(), strdatatype, compression_level);
        } catch(H5::FileIException& error) {
            close_file();
            Rf_error("c++ exception createDataset (File IException)");
        } catch(H5::GroupIException& error) {
            close_file();
            Rf_error("c++ exception createDataset (Group IException)");
        } catch(H5::DataSetIException& error) {
            close_file();
            Rf_error("c++ exception createDataset (DataSet IException)");
        } 
        
        return void();
    }
    

    /**
     * @brief Create an unlimited dataset with optional compression
     * @details Creates a new HDF5 dataset with unlimited dimensions, allowing it to
     * grow in size dynamically. The dataset is created with initial dimensions but can be 
     * extended using extendUnlimitedDataset(). Compression is fully compatible with
     * unlimited datasets.
     * 
     * Features:
     * - Unlimited dimensions in both directions (rows and columns)
     * - Chunked storage for efficiency (required for unlimited datasets)
     * - Configurable initial size with dynamic growth capability
     * - Full compression support with automatic chunk-based compression
     * - Optimal for datasets with unknown final size
     * 
     * Compression considerations:
     * - Chunking is always enabled (required for unlimited datasets)
     * - Chunk size equals initial dimensions for optimal performance
     * - Compression is applied per chunk, maintaining efficiency during growth
     * - Best compression achieved when chunks are reasonably sized (>1KB)
     * 
     * @param rows Initial number of rows (used as chunk size)
     * @param cols Initial number of columns (used as chunk size)  
     * @param strdatatype Data type for the dataset ("int", "numeric")
     *                    Note: String type not supported for unlimited datasets
     * @param compression_level Compression level (0=none, 1-9=gzip level, 6=balanced default)
     *                         Recommended: 4-6 for unlimited datasets to balance 
     *                         compression ratio with extension performance
     * 
     * @throws H5::FileIException if file operations fail
     * @throws H5::GroupIException if group operations fail
     * @throws H5::DataSetIException if dataset operations fail
     * 
     * @note Only numeric data types supported ("int", "numeric")
     * @note Use extendUnlimitedDataset() to grow dataset size
     * @note Chunking configuration affects both performance and compression efficiency
     * @note Initial dimensions should represent typical data access patterns
     * 
     * @see extendUnlimitedDataset() for growing dataset dimensions
     * @see createDataset() for fixed-size datasets with compression
     * @since Added compression support in version X.X.X
     */
    virtual void createUnlimitedDataset(size_t rows, size_t cols, std::string strdatatype, int compression_level = -1) 
    {
        try
        {
            
            H5::Exception::dontPrint();
            
            herr_t status;
            hid_t cparms; 
            std::string fullDatasetPath = groupname + "/" + name;
            
            int eff_compression = (compression_level < 0) ? this->compression_level : compression_level;
            if (eff_compression < 0) eff_compression = 6;  // default si nadie lo fijó
            
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
                Rf_error("c++ exception createUnlimitedDataset (setting chunk IException)");
                return void();
            }
            
            // Configure compression
            // if (compression_level > 0) {
            if (eff_compression > 0) {
                // H5Pset_deflate(cparms, compression_level);  // gzip compression
                H5Pset_deflate(cparms, eff_compression);  // gzip compression
                H5Pset_shuffle(cparms);                     // shuffle filter for better compression
            }
            
            if( !exists_HDF5_element(pfile, groupname) ) {
                create_HDF5_groups(groupname);
            }
            
            bool bexists = exists_HDF5_element(pfile, fullDatasetPath);
            if( bexists == true && boverwrite == false) {
                close_file();
                Rf_error("Dataset exits, please set overwrite = true to overwrite the existing dataset (createUnlimitedDataset)");
            } else {
                
                if(  bexists == true && boverwrite == true) {
                    remove_elements(pfile, fullDatasetPath); 
                }
                
                // Create dataset
                if( strdatatype == "int") {
                    H5::IntType datatype( H5::PredType::NATIVE_INT );
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace, cparms));
                } else {
                    H5::FloatType datatype( H5::PredType::NATIVE_DOUBLE ); 
                    pdataset = new H5::DataSet(pfile->createDataSet( fullDatasetPath, datatype, dataspace, cparms));
                }
            }
            
            dataspace.close();
            addAttribute( "internal", Rcpp::wrap("0") );
            
        } catch(H5::FileIException& error) {
            close_file();
            Rf_error("c++ exception createUnlimitedDataset (File IException)");
        } catch(H5::GroupIException& error) {
            close_file();
            Rf_error("c++ exception createUnlimitedDataset (Group IException)");
        } catch(H5::DataSetIException& error) {
            close_file();
            Rf_error("c++ exception createUnlimitedDataset (DataSet IException)");
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
                
#ifdef _OPENMP
                HDF5ThreadSafety::LockGuard lock_guard;
#endif
                
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
                Rf_error(" Dataset is not an unlimited dataset, fixed datasets can't be extended");
                return void();
            }
            
        } catch(H5::FileIException& error) {
            close_dataset_file();
            Rf_error("c++ exception extendUnlimitedDataset (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception extendUnlimitedDataset (DataSet IException)");
        }
        return void();
    }
    
    
    // Open an existing hdf5 dataSet
    H5::DataSet* openDataset()
    {
        try
        {
            
            std::string fullPath = groupname + "/" + name;
            
            // Check if file pointer != nullptr
            if( !pfile)  {
                Rf_error("c++ exception Please create file before proceed");
            } else { 
                
                H5::Exception::dontPrint();
                
                bool bexists = exists_HDF5_element(pfile, fullPath);
                if( bexists ) {
                    if ((pdataset == NULL) == TRUE) {
                        pdataset = new H5::DataSet(pfile->openDataSet(fullPath));    
                    }
                    
                    try {
                        H5T_class_t type_class = pdataset->getTypeClass();
                        if (type_class == H5T_INTEGER || type_class == H5T_FLOAT) {
                            type = "real";
                        } else if (type_class == H5T_STRING) {
                            type = "string";
                        } else {
                            type = "unknown";
                        }
                    } catch (const H5::Exception& e) {
                        type = "unknown";
                    }
                    
                    getDimensExistingDataset();
                    getAttribute("internal");
                    
                    // Read actual compression level from the dataset creation property list.
                    // This keeps compression_level consistent with what is on disk:
                    // after openDataset() getCompressionLevel() always reflects reality.
                    {
                        hid_t dcpl_id = H5Dget_create_plist(pdataset->getId());
                        if (dcpl_id >= 0) {
                            compression_level = 0;  // default: no compression found
                            int nfilters = H5Pget_nfilters(dcpl_id);
                            for (int fi = 0; fi < nfilters; fi++) {
                                unsigned int flags = 0, filter_config = 0;
                                size_t cd_nelmts = 1;
                                unsigned int cd_values[1] = {0};
                                char filter_name[64] = {'\0'};
                                H5Z_filter_t ftype = H5Pget_filter2(
                                    dcpl_id, fi, &flags, &cd_nelmts,
                                    cd_values, sizeof(filter_name),
                                    filter_name, &filter_config);
                                if (ftype == H5Z_FILTER_DEFLATE && cd_nelmts > 0) {
                                    compression_level = static_cast<int>(cd_values[0]);
                                    break;
                                }
                            }
                            H5Pclose(dcpl_id);
                        }
                    }
                    
                } else {
                    close_file();
                    // std::cerr<<"\nc++ exception, please create Dataset before proceed\n";
                    Rf_error("c++ exception, please create Dataset before proceed");
                    // return(pdataset);
                }
            }
            
        } catch(H5::FileIException& error) {
            close_file();
            Rf_error("c++ exception openDataset (File IException)");
        } catch(H5::GroupIException& error) {
            close_file();
            Rf_error("c++ exception openDataset (File GroupIException)");
        } catch(H5::DataSetIException& error) {
            close_file();
            Rf_error("c++ exception openDataset (File DataSetIException)");
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
            // H5::Exception::dontPrint();
            std::vector<int> dims;
            
            if(Rcpp::is<Rcpp::NumericMatrix>(DatasetValues) || Rcpp::is<Rcpp::IntegerMatrix>(DatasetValues) ) {
                
                hsize_t ncolRObject = Rcpp::as<Rcpp::NumericMatrix>(DatasetValues).ncol();
                hsize_t nrowRObject = Rcpp::as<Rcpp::NumericMatrix>(DatasetValues).nrow();
                
                if( (nrowRObject == dimDataset[1] && ncolRObject == dimDataset[0]) || ncolRObject == 1 || nrowRObject == 1 ) {
                    
                    if( nrowRObject * ncolRObject > dimDataset[0] * dimDataset[1]) {
                        Rcpp::Rcout<<"\n Data you are trying to write is bigger than existing hdf5 dataset size\n";
                        return void();
                    } else {
                        
                        // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                        HDF5ThreadSafety::LockGuard lock_guard;
#endif
                        
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
                
                // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                HDF5ThreadSafety::LockGuard lock_guard;
#endif
                
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
                
                // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                HDF5ThreadSafety::LockGuard lock_guard;
#endif
                
                // Create the memory datatype.
                H5::CompType strtype(sizeof(names));
                strtype.insertMember("chr", HOFFSET(names, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
                
                pdataset->write(convert_DataFrame_to_RangeList(DatasetValues, true), strtype);
                    
            } else {
                Rf_error("Matrix data type not allowed (writeDataset)");
            }
            
        } catch(H5::FileIException& error) {
            close_dataset_file();
            Rf_error("c++ exception writeDataset (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeDataset (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeDataset (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeDataset (DataSpace IException)");
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeDataset (Data TypeIException)");
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
            
            // H5::Exception::dontPrint();

            std::vector<int> dims;
            
            if( sizeof(mdata) > dimDataset[0] * dimDataset[1]) {
                Rcpp::Rcout<<"\n Data you are trying to write is bigger than existing hdf5 dataset size\n";
                return void();
            } else {
                
                // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                HDF5ThreadSafety::LockGuard lock_guard;
#endif
                
                H5::DataSpace dataspace(RANK2, dimDataset);
                // std::vector<double> matHiCValues = Rcpp::as<std::vector<double> >(Rcpp::as<Rcpp::NumericMatrix>(DatasetValues));
                pdataset->write( &mdata[0] , H5::PredType::NATIVE_DOUBLE);
                dataspace.close();
                
            }
            
        } catch(H5::FileIException& error) {
            close_dataset_file();
            Rf_error("c++ exception writeDataset (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeDataset (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeDataset (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeDataset (DataSpace IException)");
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeDataset (Data TypeIException)");
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
            // H5::Exception::dontPrint();
            
            hsize_t hsOffset[2], hsCount[2], hsStride[2], hsBlock[2];
                
            if(vOffset[0] + vCount[0] <= dimDataset[0] || vOffset[1] + vCount[1] <= dimDataset[1]) {
                
                // Specify size and shape of subset to write
                hsOffset[0] = vOffset[0]; hsOffset[1] = vOffset[1];
                hsStride[0] = vStride[0]; hsStride[1] = vStride[1]; // default 1
                hsBlock[0] = vBlock[0]; hsBlock[1] = vBlock[1]; // default 1
                
                hsCount[0] = vCount[0]; hsCount[1] = vCount[1];
                
                // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                HDF5ThreadSafety::LockGuard lock_guard;
#endif
                
                H5::DataSpace dataspace(RANK2, hsCount);
                H5::DataSpace memspace(RANK2, hsCount, NULL);
                
                dataspace = pdataset->getSpace();
                dataspace.selectHyperslab( H5S_SELECT_SET, hsCount, hsOffset, hsStride, hsBlock);
                
                pdataset->write(&DatasetValues.data()[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
                memspace.close();
                dataspace.close();
                    
            } else {
#ifdef _OPENMP
                HDF5ThreadSafety::LockGuard lock_guard;
#endif
                close_dataset_file();
                Rf_error("It is not possible to write block in current position (writeRowMajorDatasetBlock)");
            }
                
        } catch(H5::FileIException& error) {
            close_dataset_file();
            Rf_error("c++ exception writeRowMajorDatasetBlock (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeRowMajorDatasetBlock (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeRowMajorDatasetBlock (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeRowMajorDatasetBlock (DataSpace IException)");
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeRowMajorDatasetBlock (Data TypeIException)");
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
            // H5::Exception::dontPrint();
            
            hsize_t hsOffset[2], hsCount[2], hsStride[2], hsBlock[2];
            
            if(vOffset[0] + DatasetValues.rows() <= dimDataset[0] || vOffset[1] + DatasetValues.cols() <= dimDataset[1]) {
                
                // Specify size and shape of subset to write
                hsOffset[0] = vOffset[0]; hsOffset[1] = vOffset[1];
                hsStride[0] = vStride[0]; hsStride[1] = vStride[1]; // default 1
                hsBlock[0] = vBlock[0]; hsBlock[1] = vBlock[1]; // default 1
                
                hsCount[0] = DatasetValues.rows(); hsCount[1] = DatasetValues.cols();
                
                // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                HDF5ThreadSafety::LockGuard lock_guard;
#endif
                
                H5::DataSpace dataspace(RANK2, hsCount);
                H5::DataSpace memspace(RANK2, hsCount, NULL);
                
                dataspace = pdataset->getSpace();
                dataspace.selectHyperslab( H5S_SELECT_SET, hsCount, hsOffset, hsStride, hsBlock);
                
                pdataset->write(&DatasetValues.data()[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
                memspace.close();
                dataspace.close();
                
            } else {
                Rf_error("It is not possible to write block in current position (writeColMajorDatasetBlock)");
            }
            
        } catch(H5::FileIException& error) {
            close_dataset_file();
            Rf_error("c++ exception writeColMajorDatasetBlock (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeColMajorDatasetBlock (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeColMajorDatasetBlock (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeColMajorDatasetBlock (DataSpace IException)");
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception writeColMajorDatasetBlock (Data TypeIException)");
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
            // H5::Exception::dontPrint();

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
                    
                    // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                    HDF5ThreadSafety::LockGuard lock_guard;
#endif
                    
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
                    Rf_error("It is not possible to write block in current position (writeDatasetBlock)");
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

                    // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                    HDF5ThreadSafety::LockGuard lock_guard;
#endif
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
                    Rf_error("It is not possible to write block in current position (writeDatasetBlock)");
                }
            } else {
                Rf_error("Matrix data type not allowed (writeDatasetBlock)");
            }

        } catch(H5::FileIException& error) {
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock (DataSpace IException)");
        } catch(H5::DataTypeIException& error) { 
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock (Data TypeIException)");
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
            // H5::Exception::dontPrint();
            
            hsize_t hsOffset[2], hsCount[2], hsStride[2], hsBlock[2];
            
            // Specify size and shape of subset to write
            hsStride[0] = vStride[0]; hsStride[1] = vStride[1]; // default 1
            hsBlock[0] = vBlock[0]; hsBlock[1] = vBlock[1]; // default 1
                
            hsOffset[0] = vOffset[0]; hsOffset[1] = vOffset[1];    
            hsCount[0] = vCount[0]; hsCount[1] = vCount[1];
            
            
            if(vOffset[0] + hsCount[0] <= dimDataset[0] && vOffset[1] + hsCount[1] <= dimDataset[1] && 
               DatasetValues.size()<= dimDatasetinFile[0] * dimDatasetinFile[1]) {
                
                // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
                HDF5ThreadSafety::LockGuard lock_guard;
#endif
                
                H5::DataSpace dataspace(RANK2, hsCount);
                H5::DataSpace memspace(RANK2, hsCount, NULL);
                
                dataspace = pdataset->getSpace();
                dataspace.selectHyperslab( H5S_SELECT_SET, hsCount, hsOffset, hsStride, hsBlock);
                
                pdataset->write(&DatasetValues[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
                memspace.close();
                dataspace.close();
                
            } else {
                
                Rf_error("It is not possible to write block in current position (writeDatasetBlock)");
            }
                
        } catch(H5::FileIException& error) {
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock std::vector (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock std::vector (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock std::vector (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock std::vector (DataSpace IException)");
        } catch(H5::DataTypeIException& error) { 
            close_dataset();
            close_dataset_file();
            Rf_error("c++ exception writeDatasetBlock std::vector (Data TypeIException)");
        }
        return void();
    }
    
    
    /**
     * @brief Writes dimension names alongside this dataset
     *
     * @details Stores row and column names in the hidden group
     * @c group/.<dataset>_dimnames following the BigDataStatMeth convention:
     * - @c /1 stores row names
     * - @c /2 stores column names
     *
     * Unlike @c hdf5Dims::writeDimnames(), this method uses the existing open
     * file handle (@c this->pfile) instead of opening a second independent
     * handle.  This prevents HDF5 metadata cache staleness for long-lived
     * objects such as those managed by the R6 class.
     *
     * If the @c _dimnames group already exists it is removed and recreated so
     * that names are always fully replaced.
     *
     * @param rownames Character vector of row names. Pass @c character(0) to
     *                 skip writing row names.
     * @param colnames Character vector of column names. Pass @c character(0)
     *                 to skip writing column names.
     *
     * @note Length of @p rownames is validated against @c ncols() (= nrow_R).
     *       Length of @p colnames is validated against @c nrows() (= ncol_R).
     *       Mismatched lengths are silently ignored.
     *
     * @throws Rf_error on HDF5 I/O failures
     */
    void writeDimnames(Rcpp::StringVector rownames, Rcpp::StringVector colnames)
    {
        const std::string orig_groupname = groupname;       // ← AÑADIR: guarda el original
        const std::string dimgrp   = groupname + "/." + name + "_dimnames";
        const std::string path_row = dimgrp + "/1";
        const std::string path_col = dimgrp + "/2";
        
        try {
            H5::Exception::dontPrint();
            
            if (exists_HDF5_element(pfile, dimgrp))
                remove_elements(pfile, groupname, {"." + name + "_dimnames"});
            
            create_HDF5_groups(dimgrp);
            groupname = orig_groupname;                     // ← AÑADIR: restaura tras la corrupción
            
            if (rownames.length() > 0 &&
                static_cast<hsize_t>(rownames.length()) == ncols())
                writeStringDataset(path_row, rownames);
            
            if (colnames.length() > 0 &&
                static_cast<hsize_t>(colnames.length()) == nrows())
                writeStringDataset(path_col, colnames);
            
        } catch (H5::FileIException& e) {
            groupname = orig_groupname;                     // ← AÑADIR: restaura en error también
            Rf_error("c++ exception hdf5Dataset::writeDimnames (File IException)");
        } catch (H5::DataSetIException& e) {
            groupname = orig_groupname;
            Rf_error("c++ exception hdf5Dataset::writeDimnames (DataSet IException)");
        } catch (H5::GroupIException& e) {
            groupname = orig_groupname;
            Rf_error("c++ exception hdf5Dataset::writeDimnames (Group IException)");
        } catch (std::exception& e) {
            groupname = orig_groupname;
            Rf_error("c++ exception hdf5Dataset::writeDimnames: %s", e.what());
        }
    }
    
    
    /**
     * @brief Reads dimension names stored alongside this dataset
     *
     * @details Reads row and column names from the hidden group
     * @c group/.<dataset>_dimnames created by @c writeDimnames().  Uses the
     * existing open file handle (@c this->pfile) — no new handle is opened —
     * so the result is always consistent with the current file state regardless
     * of when the dataset was opened.
     *
     * Returns @c character(0) for any component that has not been written.
     *
     * @return Named @c Rcpp::List with two elements:
     * \describe{
     *   \item{rownames}{CharacterVector of row names, or @c character(0) if absent}
     *   \item{colnames}{CharacterVector of col names, or @c character(0) if absent}
     * }
     *
     * @throws Rf_error on HDF5 I/O failures
     */
    Rcpp::List readDimnames()
    {
        const std::string path_row =
            groupname + "/." + name + "_dimnames/1";
        const std::string path_col =
            groupname + "/." + name + "_dimnames/2";
        
        try {
            H5::Exception::dontPrint();
            
            return Rcpp::List::create(
                Rcpp::Named("rownames") = readStringDataset(path_row),
                Rcpp::Named("colnames") = readStringDataset(path_col));
            
        } catch (H5::FileIException& e) {
            Rf_error("c++ exception hdf5Dataset::readDimnames (File IException)");
        } catch (H5::DataSetIException& e) {
            Rf_error("c++ exception hdf5Dataset::readDimnames (DataSet IException)");
        } catch (std::exception& e) {
            Rf_error("c++ exception hdf5Dataset::readDimnames: %s", e.what());
        }
        
        // Fallback — only reached if Rf_error somehow doesn't throw
        return Rcpp::List::create(
            Rcpp::Named("rownames") = Rcpp::CharacterVector(0),
            Rcpp::Named("colnames") = Rcpp::CharacterVector(0));
    }
    

    /**
     * @brief Create subset dataset with selected rows/columns (memory efficient)
     * @details Creates a new dataset containing only the specified rows or columns
     * using HDF5's hyperslab selection for direct disk-to-disk copy without loading
     * data into memory. Ideal for big datasets.
     * 
     * @param indices Vector of row/column indices to include (0-based)
     * @param select_rows If true, selects rows; if false, selects columns
     * @param new_group Target group for the new dataset (default: same group)
     * @param new_name Name for the new dataset (default: original_name + "_subset")
     * 
     * @since 0.99.0
     */
    virtual void createSubsetDataset(const std::vector<int>& indices, 
                                     bool select_rows = true,
                                     const std::string& new_group = "",
                                     const std::string& new_name = "")
    {
        try {
            // H5::Exception::dontPrint();
            
            if (!pdataset || indices.empty()) {
                Rf_error("Dataset not open or indices empty");
                return void();
            }
            
            // Get dimensions from R perspective (siguiendo lÃ³gica de la clase)
            hsize_t rows_from_r = nrows_r();  // dimDatasetinFile[1] - filas como las ve R
            hsize_t cols_from_r = ncols_r();  // dimDatasetinFile[0] - columnas como las ve R
            
            // Get actual file dimensions for hyperslab operations
            hsize_t file_rows = nrows_file();  // dimDatasetinFile[0]
            hsize_t file_cols = ncols_file();  // dimDatasetinFile[1]
            
            // Validate indices against R perspective
            hsize_t max_index = select_rows ? rows_from_r : cols_from_r;
            for (int idx : indices) {
                if (idx < 0 || idx >= max_index) {
                    Rf_error("Index %d out of bounds", idx);
                    return void();
                }
            }
            
            // Set target group and name
            std::string target_group = new_group.empty() ? groupname : new_group;
            std::string target_name = new_name.empty() ? name + "_subset" : new_name;
            
            // Calculate new dimensions from R perspective
            hsize_t new_rows_r = select_rows ? indices.size() : rows_from_r;
            hsize_t new_cols_r = select_rows ? cols_from_r : indices.size();
            
#ifdef _OPENMP
            HDF5ThreadSafety::LockGuard lock_guard;
#endif
            
            // Create destination dataset (usando dimensiones como las espera R)
            BigDataStatMeth::hdf5Dataset dest_dataset(pfile, target_group, target_name, true);
            dest_dataset.createDataset(new_rows_r, new_cols_r, type);
            
            // Get source and destination dataspaces
            H5::DataSpace src_space = pdataset->getSpace();
            H5::DataSpace dest_space = dest_dataset.getDatasetptr()->getSpace();
            
            if (select_rows) {
                // Usuario quiere filas de R = columnas en archivo HDF5
                for (size_t i = 0; i < indices.size(); ++i) {
                    // Select source column in file (filas de R son columnas en archivo)
                    hsize_t src_offset[2] = {0, static_cast<hsize_t>(indices[i])};
                    hsize_t count[2] = {file_rows, 1};
                    src_space.selectHyperslab(H5S_SELECT_SET, count, src_offset);
                    // Select destination column
                    hsize_t dest_offset[2] = {0, i};
                    dest_space.selectHyperslab(H5S_SELECT_SET, count, dest_offset);
                    // Create memory space and copy
                    H5::DataSpace mem_space(2, count);
                    std::vector<double> buffer(file_rows);
                    pdataset->read(buffer.data(), H5::PredType::NATIVE_DOUBLE, mem_space, src_space);
                    dest_dataset.getDatasetptr()->write(buffer.data(), H5::PredType::NATIVE_DOUBLE, mem_space, dest_space);
                }
            } else {
                // Usuario quiere columnas de R = filas en archivo HDF5
                for (size_t i = 0; i < indices.size(); ++i) {
                    // Select source row in file (columnas de R son filas en archivo)
                    hsize_t src_offset[2] = {static_cast<hsize_t>(indices[i]), 0};
                    hsize_t count[2] = {1, file_cols};
                    src_space.selectHyperslab(H5S_SELECT_SET, count, src_offset);
                    
                    // Select destination row
                    hsize_t dest_offset[2] = {i, 0};
                    dest_space.selectHyperslab(H5S_SELECT_SET, count, dest_offset);
                    
                    // Create memory space and copy
                    H5::DataSpace mem_space(2, count);
                    std::vector<double> buffer(file_cols);
                    pdataset->read(buffer.data(), H5::PredType::NATIVE_DOUBLE, mem_space, src_space);
                    dest_dataset.getDatasetptr()->write(buffer.data(), H5::PredType::NATIVE_DOUBLE, mem_space, dest_space);
                }
            }
            
            // Rcpp::Rcout << "Created subset dataset: " << target_group << "/" << target_name 
            //             << " (" << new_rows_r << "x" << new_cols_r << " from R perspective)" << std::endl;
            
        } catch (H5::DataSetIException& error) {
            close_dataset_file();
            Rf_error("c++ exception createSubsetDataset (DataSet IException)");
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
            // H5::Exception::dontPrint();
            
            hsize_t offset[2], count[2], stride[2], block[2];
            
            offset[0] = ivoffset[0]; offset[1] = ivoffset[1];
            count[0] = ivcount[0]; count[1] = ivcount[1];
            stride[0] = ivstride[0]; stride[1] = ivstride[1];
            block[0] = ivblock[0]; block[1] = ivblock[1];
            
            // Define Memory Dataspace. Get file dataspace and select a subset from the file dataspace.
            hsize_t dimsm[2];
            dimsm[0] = count[0]; 
            dimsm[1] = count[1];
            
            // Thread-safe HDF5 I/O operations
#ifdef _OPENMP
            HDF5ThreadSafety::LockGuard lock_guard;
#endif
            
            H5::DataSpace memspace(RANK2, dimsm, NULL);
            
            //  Get dataset dataspace.
            H5::DataSpace dataspace = pdataset->getSpace();
            dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block); 
            
            H5T_class_t type_class = pdataset->getTypeClass();
            
            // Get class of datatype and print message if it's an integer.
            if( type_class == H5T_INTEGER || type_class == H5T_FLOAT ) {
                pdataset->read( rdatablock, H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
            } else {
                Rf_error("c++ exception readDatasetBlock (Data type not allowed, maybe are trying to read string matrix?)");
                return void();
            }
            
            memspace.close();
            dataspace.close();
            
        } catch( H5::FileIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception readDatasetBlock (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception readDatasetBlock (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception readDatasetBlock (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception readDatasetBlock (DataSpace IException)");
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception readDatasetBlock (Data TypeIException)");
        } catch(std::exception &ex) {
            close_dataset_file();
            Rf_error("c++ exception readDatasetBlock: %s", ex.what());
        } catch (...) {
            close_dataset_file();
            Rf_error("C++ exception readDatasetBlock (unknown reason)");
        }
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
            // error.printErrorStack();
             Rf_error("c++ exception addAttribute (DataSpace IException)");
            return -1;
        } catch( H5::AttributeIException& error ) {
            // error.printErrorStack();
            Rf_error("c++ exception addAttribute (Attribute IException)");
            return -1;
        } catch( H5::DataSetIException& error ) { 
            // error.printErrorStack();
            Rf_error("c++ exception addAttribute (DataSet IException)");
            return -1;
        } catch( H5::FileIException& error ) {
            // error.printErrorStack();
            Rf_error("c++ exception addAttribute (File IException)");
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
                    Rf_error("Error: Unknown data type");
                    // Rcpp::Rcout<<"\nPos No SEP"<<type->getClass()<<"\n";
                }
                
                attr->close();
                
            }
            
        } catch( H5::FileIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception getAttribute (File IException)");
        } catch(H5::DataSetIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception getAttribute (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception getAttribute (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception getAttribute (DataSpace IException)");
        } catch(H5::DataTypeIException& error) { 
            close_dataset_file();
            Rf_error("c++ exception getAttribute (Data TypeIException)");
        }
        return void();
    }
    
    
    /**
     * @brief Move dataset to a new location within the HDF5 file
     * @details Moves the dataset and its associated rownames/colnames datasets to a new
     * location within the same HDF5 file. The operation also updates the internal class
     * variables to reflect the new location and reopens the dataset at the new path.
     * 
     * Features:
     * - Moves main dataset using H5Lmove for efficiency
     * - Automatically moves associated rownames and colnames datasets
     * - Creates parent groups automatically if they don't exist
     * - Updates internal class variables (name, groupname)
     * - Reopens dataset at new location
     * - Optional overwrite of destination
     * - Preserves all dataset attributes and properties
     * 
     * @param new_path New complete path for the dataset (e.g., "/new_group/new_name")
     * @param overwrite Whether to overwrite destination if it exists (default: false)
     * 
     * @throws H5::Exception if HDF5 operations fail
     * @throws H5::DataSetIException if dataset operations fail
     * @throws H5::GroupIException if group operations fail
     * 
     * @note The dataset must be currently open
     * @note Parent groups will be created automatically if they don't exist
     * @note Associated rownames/colnames datasets are moved to the same new group
     * @note All internal class variables are updated to reflect the new location
     * 
     * @example
     * ```cpp
     * dataset.moveDataset("/new_group/renamed_dataset", true);
     * ```
     * 
     * @since 0.99.0
     */
    virtual void moveDataset(const std::string& new_path, bool overwrite = false)
    {
        try {
            H5::Exception::dontPrint();
            
            // Validate input
            if (new_path.empty()) {
                Rf_error("New path cannot be empty");
                return void();
            }
            
            if (!pdataset) {
                Rf_error("Dataset is not open");
                return void();
            }
            
            // Parse new path to extract group and dataset name
            fullpath newLocation = SplitElementName(new_path);
            std::string new_group = newLocation.path;
            std::string new_name = newLocation.filename;
            
            // Build current and destination paths
            std::string current_path = groupname + "/" + name;
            std::string dest_path = new_group + "/" + new_name;
            
            // Check if moving to same location
            if (current_path == dest_path) {
                Rcpp::Rcout << "Dataset is already at the specified location" << std::endl;
                return void();
            }
            
            // Create parent groups if they don't exist
            if (new_group != "/" && !exists_HDF5_element(pfile, new_group)) {
                
                create_HDF5_groups(new_group);
            }
            
            // Check if destination exists
            bool dest_exists = exists_HDF5_element(pfile, dest_path);
            
            if (dest_exists && !overwrite) {
                Rcpp::Rcout << "Destination dataset already exists: " << dest_path 
                            << ". Set overwrite = TRUE to replace it." << std::endl;
                return void();
            }
            
            // Remove destination if it exists and overwrite is true
            if (dest_exists && overwrite) {
                // Parse destination path to extract group and element name for remove_elements
                std::string dest_group = "/";
                std::string dest_element = dest_path;
                
                size_t last_slash = dest_path.find_last_of('/');
                if (last_slash != std::string::npos && last_slash > 0) {
                    dest_group = dest_path.substr(0, last_slash);
                    dest_element = dest_path.substr(last_slash + 1);
                } else if (last_slash == 0) {
                    dest_group = "/";
                    dest_element = dest_path.substr(1);
                }
                
                remove_elements(pfile, dest_group, {dest_element});
            }
            
            // Close current dataset before moving
            if (pdataset) {
                pdataset->close();
                delete pdataset;
                pdataset = nullptr;
            }
            
            // Get file ID for H5Lmove operations
            hid_t file_id = pfile->getId();
            
            // Move main dataset
            herr_t move_status = H5Lmove(file_id, current_path.c_str(), 
                                         file_id, dest_path.c_str(), 
                                         H5P_DEFAULT, H5P_DEFAULT);
            
            if (move_status < 0) {
                Rf_error("Failed to move dataset from %s to %s", current_path.c_str(), dest_path.c_str());
                return void();
            }
            
            // Move associated rownames dataset if it exists
            if (!rownamesDataset.empty()) {
                std::string current_rownames = groupname + "/" + rownamesDataset;
                std::string new_rownames = new_group + "/" + rownamesDataset;
                
                if (exists_HDF5_element(pfile, current_rownames)) {
                    herr_t rownames_status = H5Lmove(file_id, current_rownames.c_str(), 
                                                     file_id, new_rownames.c_str(), 
                                                     H5P_DEFAULT, H5P_DEFAULT);
                    if (rownames_status < 0) {
                        Rcpp::Rcerr << "Warning: Failed to move rownames dataset" << std::endl;
                    }
                } else {
                    Rcpp::Rcout<<"\nNo rownames to move";
                }
            } else {
                Rcpp::Rcout<<"\nNo rownames to move";
            }
            
            // Move associated colnames dataset if it exists
            if (!colnamesDataset.empty()) {
                std::string current_colnames = groupname + "/" + colnamesDataset;
                std::string new_colnames = new_group + "/" + colnamesDataset;
                
                if (exists_HDF5_element(pfile, current_colnames)) {
                    herr_t colnames_status = H5Lmove(file_id, current_colnames.c_str(), 
                                                     file_id, new_colnames.c_str(), 
                                                     H5P_DEFAULT, H5P_DEFAULT);
                    if (colnames_status < 0) {
                        Rcpp::Rcerr << "Warning: Failed to move colnames dataset" << std::endl;
                    }
                } else {
                    Rcpp::Rcout<<"\nNo colnames to move";
                }
            } else {
                Rcpp::Rcout<<"\nNo colnames to move";
            }
            
            // Update internal class variables
            groupname = new_group;
            name = new_name;
            
            // Reopen dataset at new location
            pdataset = new H5::DataSet(pfile->openDataSet(dest_path));
            
            // Update dimensions in case they need refresh
            getDimensExistingDataset();
            
        } catch (H5::FileIException& error) {
            close_dataset_file();
            Rf_error("c++ exception moveDataset (File IException)");
        } catch (H5::DataSetIException& error) {
            close_dataset_file();
            Rf_error("c++ exception moveDataset (DataSet IException)");
        } catch (H5::GroupIException& error) {
            close_dataset_file();
            Rf_error("c++ exception moveDataset (Group IException)");
        } catch (H5::DataSpaceIException& error) {
            close_dataset_file();
            Rf_error("c++ exception moveDataset (DataSpace IException)");
        } catch (H5::DataTypeIException& error) {
            close_dataset_file();
            Rf_error("c++ exception moveDataset (Data TypeIException)");
        }
        
        return void();
    }
    
    
    /**
     * @brief Copy this dataset to a dataset in a different (or same) HDF5 file.
     * @details Reads the source block-by-block and writes to the destination,
     * so memory usage is O(blockrows * ncols) regardless of matrix size.
     * The destination file is created if it does not exist yet.
     * If source and destination file are the same, prefer moveDataset() or
     * a simple HDF5 group copy instead.
     *
     * @param dst_file   Path to the destination HDF5 file.
     * @param dst_group  Group path in the destination file.
     * @param dst_name   Dataset name in the destination file.
     * @param blockrows  Rows per copy block (default 500).
     * @param overwrite  If true, overwrite existing destination dataset.
     */
    virtual void copyToFile(const std::string& dst_file,
                            const std::string& dst_group,
                            const std::string& dst_name,
                            hsize_t blockrows = 500,
                            bool overwrite = true)
    {
        
        try {
            
            H5::Exception::dontPrint();
            
            if (!pdataset)
                Rf_error("copyToFile: source dataset is not open");
            
            const hsize_t nrows = this->nrows();
            const hsize_t ncols = this->ncols();
            
            // Pre-create destination file if it does not exist yet
            {
                std::ifstream probe(dst_file);
                if (!probe.good()) {
                    H5::H5File newf(dst_file, H5F_ACC_TRUNC);
                    newf.close();
                }
            }
            
            // Create destination dataset with same dimensions
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsDst(
                    new BigDataStatMeth::hdf5Dataset(dst_file, dst_group, dst_name, overwrite)
            );
            dsDst->createDataset(nrows, ncols, "real");
            if (dsDst->getDatasetptr() == nullptr)
                Rf_error("copyToFile: cannot create destination '%s/%s' in '%s'",
                         dst_group.c_str(), dst_name.c_str(), dst_file.c_str());
            
            // Block-by-block copy
            std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
            for (hsize_t row = 0; row < nrows; row += blockrows) {
                const hsize_t toread = (row + blockrows <= nrows)
                ? blockrows : (nrows - row);
                std::vector<double> buf(toread * ncols);
                this->readDatasetBlock({row, 0}, {toread, ncols},
                                       stride, block, buf.data());
                Eigen::Map<Eigen::Matrix<double,
                                         Eigen::Dynamic, Eigen::Dynamic,
                                         Eigen::RowMajor>> X(buf.data(), toread, ncols);
                dsDst->writeRowMajorDatasetBlock(X,
                                                 {row, 0}, {toread, ncols},
                                                 stride, block);
            }
            
        } catch (H5::FileIException& error) {
            Rf_error("c++ exception copyToFile (File IException): %s",
                     error.getCDetailMsg());
        } catch (H5::DataSetIException& error) {
            Rf_error("c++ exception copyToFile (DataSet IException): %s",
                     error.getCDetailMsg());
        } catch (H5::DataSpaceIException& error) {
            Rf_error("c++ exception copyToFile (DataSpace IException): %s",
                     error.getCDetailMsg());
        } catch (std::exception& ex) {
            Rf_error("c++ exception copyToFile: %s", ex.what());
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
    
    /**
     * @brief Set rownames path inside hdf5 data file
     * @return void
     */
    void setRownamesDatasetPath(std::string fullpath) { rownamesDataset = fullpath; } // set rownames route 
    
    /**
     * @brief Set colnames path inside hdf5 data file
     * @return void
     */
    void setColnamesDatasetPath(std::string fullpath) { colnamesDataset = fullpath; } // set colnames route
    
    /**
     * @brief Get dataset type
     * @return String representing the dataset type (e.g., "int", "numeric", "string")
     */
    std::string getDatatype() { return(type); }
    

    // /**
    //  * @brief Reads dimension names stored alongside an HDF5 dataset
    //  *
    //  * Opens the hidden @c _dimnames group written by @c writeDimnames() and
    //  * reads back the row-names vector (@c /1) and the column-names vector
    //  * (@c /2).  If a vector does not exist an empty @c Rcpp::CharacterVector
    //  * is returned for that component.
    //  *
    //  * @return Named @c Rcpp::List with two elements:
    //  *   - @c "rownames" : @c CharacterVector of row names (length 0 if absent)
    //  *   - @c "colnames" : @c CharacterVector of col names (length 0 if absent)
    //  *
    //  * @pre The @c hdf5Dims object must have been constructed with
    //  *      @c bWrite = false (read mode).
    //  * @pre The underlying hdf5Dataset must be open.
    //  *
    //  * @throws std::runtime_error on HDF5 I/O failures
    //  */
    // Rcpp::List readDimnames()
    // {
    //     // Path convention: group/.<dataset>_dimnames/1 (rows) and /2 (cols)
    //     const std::string dimgrp   = groupname + "/." + name + "_dimnames";
    //     const std::string path_row = dimgrp + "/1";
    //     const std::string path_col = dimgrp + "/2";
    //     
    //     return Rcpp::List::create(
    //         Rcpp::Named("rownames") = readStringDataset(path_row),
    //         Rcpp::Named("colnames") = readStringDataset(path_col));
    // }
    
    // Destructor
    
    /**
     * @brief Destructor
     * @details Closes the dataset and releases resources
     */
    virtual ~hdf5Dataset() noexcept {
        try {
            if (pdataset) {
                if (H5Iis_valid(pdataset->getId()) > 0)
                    pdataset->close();
                #ifdef OWNS_PDATASET
                    delete pdataset;
                #endif
                pdataset = nullptr;
            }
        } catch (...) {}
    }
    
    // virtual ~hdf5Dataset() noexcept {
    //     try {
    //         if (pdataset) {
    //             pdataset->close();
    //             #ifdef OWNS_PDATASET
    //                 delete pdataset;
    //             #endif
    //             pdataset = nullptr;
    //         }
    //     } catch (...) {}
    // }
    
    
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
    
    /**
     * @brief Compression level for this dataset instance
     * @details Controls gzip compression for createDataset() and createUnlimitedDataset().
     * - Default value: 6 (gzip balanced, ~60-80% space savings)
     * - Range: 0 (no compression) to 9 (maximum compression, most CPU)
     * - After openDataset() is called, updated to reflect the actual compression level
     *   stored on disk (0 if the dataset was written without compression).
     * - Set via setCompressionLevel() before createDataset() to override.
     * - Read via getCompressionLevel() to propagate to output datasets.
     */
    int compression_level = -1;
    
    
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
            Rf_error("c++ exception getDimensExistingDataset (File IException)");
        } catch(H5::DataSetIException& error) { 
            Rf_error("c++ exception getDimensExistingDataset (DataSet IException)");
        } catch(H5::GroupIException& error) { 
            Rf_error("c++ exception getDimensExistingDataset (Group IException)");
        } catch(H5::DataSpaceIException& error) { 
            Rf_error("c++ exception getDimensExistingDataset (DataSpace IException)");
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
    
private:
    
    /**
     * @brief Reads a fixed-length string CompType dataset from the HDF5 file
     *
     * @details Reads a dataset stored with the BigDataStatMeth string convention
     * (CompType with a single @c chr member of length @c MAXSTRING) using the
     * existing open file handle.  Returns an empty vector when the dataset does
     * not exist, so callers do not need to check existence beforehand.
     *
     * @param path  Full HDF5 path of the string dataset to read
     *              (e.g. @c group/.<name>_dimnames/1)
     * @return      @c Rcpp::CharacterVector with the stored strings, or
     *              @c CharacterVector(0) when the path does not exist
     *
     * @throws std::runtime_error on HDF5 I/O failures
     */
    Rcpp::CharacterVector readStringDataset(const std::string& path)
    {
        try {
            H5::Exception::dontPrint();
            
            if (!exists_HDF5_element(pfile, path))
                return Rcpp::CharacterVector(0);
            
            typedef struct name_t { char chr[MAXSTRING]; } name_t;
            
            H5::DataSet   hds  = pfile->openDataSet(path);
            H5::DataSpace spc  = hds.getSpace();
            hsize_t       dims[1] = {0};
            spc.getSimpleExtentDims(dims);
            
            if (dims[0] == 0) {
                hds.close();
                return Rcpp::CharacterVector(0);
            }
            
            H5::CompType mtype(sizeof(name_t));
            mtype.insertMember("chr", HOFFSET(name_t, chr),
                               H5::StrType(H5::PredType::C_S1, MAXSTRING));
            
            std::vector<name_t> buf(dims[0]);
            hds.read(buf.data(), mtype);
            hds.close();
            
            Rcpp::CharacterVector cv(dims[0]);
            for (hsize_t i = 0; i < dims[0]; ++i)
                cv[static_cast<int>(i)] = std::string(buf[i].chr);
            
            return cv;
            
        } catch (H5::FileIException& e) {
            throw std::runtime_error( "c++ exception readStringDataset (File IException): " + e.getDetailMsg());
        } catch (H5::DataSetIException& e) {
            throw std::runtime_error( "c++ exception readStringDataset (DataSet IException): " + e.getDetailMsg());
        } catch (H5::DataSpaceIException& e) {
            throw std::runtime_error( "c++ exception readStringDataset (DataSpace IException): " + e.getDetailMsg());
        } catch (std::exception& e) {
            throw std::runtime_error( std::string("c++ exception readStringDataset: ") + e.what());
        }
    }
    
    
    /**
     * @brief Writes a fixed-length string CompType dataset to the HDF5 file
     *
     * @details Creates a new dataset at @p path using the BigDataStatMeth string
     * convention (CompType with a single @c chr member of length @c MAXSTRING).
     * Strings longer than @c MAXSTRING-1 characters are silently truncated.
     * The dataset is written in a single block; for very large string vectors
     * block-wise writing could be added following the pattern in
     * @c hdf5Dims::writeStringVector().
     *
     * @param path    Full HDF5 path where the dataset will be created.
     *                The parent group must already exist.
     * @param values  String vector to write
     *
     * @pre The parent group of @p path must exist in the file.
     * @pre @p values must be non-empty.
     *
     * @throws std::runtime_error on HDF5 I/O failures
     */
    void writeStringDataset(const std::string& path, Rcpp::StringVector values)
    {
        try {
            H5::Exception::dontPrint();
            
            typedef struct name_t { char chr[MAXSTRING]; } name_t;
            
            hsize_t       n = static_cast<hsize_t>(values.length());
            hsize_t       dims[1] = {n};
            H5::DataSpace spc(1, dims);
            
            H5::CompType mtype(sizeof(name_t));
            mtype.insertMember("chr", HOFFSET(name_t, chr),
                               H5::StrType(H5::PredType::C_S1, MAXSTRING));
            
            H5::DataSet hds = pfile->createDataSet(path, mtype, spc);
            
            std::vector<name_t> buf(n);
            for (hsize_t i = 0; i < n; ++i) {
                std::string word = Rcpp::as<std::string>(
                    values[static_cast<int>(i)]);
                hsize_t j = 0;
                for (; j < word.size() && j < MAXSTRING - 1; ++j)
                    buf[i].chr[j] = word[j];
                buf[i].chr[j] = '\0';
            }
            
            hds.write(buf.data(), mtype);
            hds.close();
            spc.close();
            
        } catch (H5::FileIException& e) {
            throw std::runtime_error(  "c++ exception writeStringDataset (File IException): " + e.getDetailMsg());
        } catch (H5::DataSetIException& e) {
            throw std::runtime_error( "c++ exception writeStringDataset (DataSet IException): " + e.getDetailMsg());
        } catch (H5::DataSpaceIException& e) {
            throw std::runtime_error( "c++ exception writeStringDataset (DataSpace IException): " + e.getDetailMsg());
        } catch (H5::DataTypeIException& e) {
            throw std::runtime_error( "c++ exception writeStringDataset (DataType IException): " + e.getDetailMsg());
        } catch (std::exception& e) {
            throw std::runtime_error( std::string("c++ exception writeStringDataset: ") + e.what());
        }
    }
    
};
}

#endif // BIGDATASTATMETH_HDF5_DATASETS_HPP
