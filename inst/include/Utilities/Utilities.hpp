/**
 * @file Utilities.hpp
 * @brief Core utility functions and data structures for BigDataStatMeth
 *
 * This file provides essential utility functions and data structures used
 * throughout the BigDataStatMeth package. It includes functionality for
 * file path handling, R object type detection, matrix operations,
 * and block size calculations for efficient memory management.
 *
 * Key features:
 * - Path and filename manipulation
 * - R object type detection and dimension extraction
 * - Matrix block size optimization
 * - SVD and QR decomposition structures
 * - Memory-efficient block operations
 * - System-aware thread optimization
 * - HPC/Linux server specific configurations
 * - Performance monitoring capabilities
 * - Intelligent block size optimization
 *
 * @note This is a core utility file that many other components depend on.
 */

#ifndef BIGDATASTATMETH_UTILITIES_HPP
#define BIGDATASTATMETH_UTILITIES_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"
#include <chrono> // For performance monitoring

namespace BigDataStatMeth {

    // Structs
    
    /**
     * @brief Structure for storing file path components
     */
    struct fullpath {
        std::string path ;
        std::string filename ;
    };

    /**
     * @brief Structure for storing SVD decomposition results
     */
    struct svdeig {
        Eigen::VectorXd d;
        Eigen::MatrixXd u;
        Eigen::MatrixXd v;
        bool bokuv = false;
        bool bokd = false;
        std::string hdf5file = "";
    };

    /**
     * @brief Structure for storing QR decomposition results
     */
    struct strQR {
        Eigen::MatrixXd Q;
        Eigen::MatrixXd R;
    };

    // Performance monitoring struct
    struct PerformanceInfo {
        double execution_time;
        int threads_used;
        std::string system_type;
        size_t memory_used;
        bool was_optimized;
    };

    
    // Functions
    
    /**
     * @brief Splits a full file path into directory path and filename
     *
     * @param str Full file path
     * @return fullpath Structure containing separated path and filename
     *
     * @details Handles both forward and backward slashes for cross-platform compatibility
     */
    inline fullpath SplitElementName (std::string str)
    {
        fullpath currentpath;
        std::size_t found = str.find_last_of("/\\");
        
        if( found< str.length() ) {
            currentpath.filename =  str.substr(found+1);
            currentpath.path = str.substr(0,found);
        }else {
            currentpath.filename = str;
            currentpath.path = "";
        }
        return(currentpath);
    }
    
    
    /**
     * @brief Determines the data type of an R object
     *
     * @param obj R object to analyze
     * @return std::string Type identifier string
     *
     * @details Supported types:
     * - "numeric": Numeric vectors/matrices
     * - "int": Integer vectors/matrices
     * - "char": Character vectors/matrices
     * - "logic": Logical vectors/matrices
     * - "dataframe": Data frames
     * - "list": Lists
     * - "S4": S4 objects
     * - "NULL": NULL objects
     * - "unknown": Unrecognized types
     *
     * @throws std::exception on type detection errors
     */
    inline std::string getObjecDataType(Rcpp::RObject obj) 
    {
        std::string strtype = "";
        
        try 
        {
            if( Rcpp::is<Rcpp::NumericVector>(obj) ) {
                strtype = "numeric";
            } else if( Rcpp::is<Rcpp::IntegerVector>(obj) ) {
                strtype = "int";
            } else if( Rcpp::is<Rcpp::CharacterVector>(obj) ) {
                strtype = "char";
            } else if( Rcpp::is<Rcpp::LogicalVector>(obj) ) {
                strtype = "logic";
            } else if( Rcpp::is<Rcpp::DataFrame>(obj) ) {
                strtype = "dataframe";
            } else if( Rcpp::is<Rcpp::List>(obj) ) {
                strtype = "list";
            } else if( obj.isS4() ) {
                strtype = "S4";
            } else if( obj.isNULL() ) {
                strtype = "NULL";
            } else {
                strtype = "unknown";
            }
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getObjecDataType: "<< ex.what() << "\n";
        }
        
        return(strtype);
    }
    
    
    /**
     * @brief Gets dimensions of an R object
     *
     * @param obj R object to measure
     * @param strtype Optional pre-determined object type
     * @return Rcpp::IntegerVector Two-element vector with dimensions
     *
     * @details
     * - Returns [1, length] for vectors
     * - Returns [rows, cols] for matrices
     * - Returns [nrow, ncol] for data frames
     * - Returns [0, 0] for unsupported types
     *
     * @throws std::exception on dimension extraction errors
     */
    inline Rcpp::IntegerVector getObjectDims(Rcpp::RObject obj, std::string strtype = "") 
    {
        Rcpp::IntegerVector dims(2);
        
        try 
        {
            if(strtype =="") {
                strtype = getObjecDataType(obj);    
            }
            
            if( strtype == "numeric"  || strtype == "int" || strtype == "factor" ){
                if( Rf_isMatrix(obj)) {
                    dims[0] = Rcpp::as<Rcpp::NumericMatrix>(obj).rows();
                    dims[1] = Rcpp::as<Rcpp::NumericMatrix>(obj).cols();
                } else {
                    dims[0] = 1;
                    dims[1] = Rcpp::as<Rcpp::NumericVector>(obj).length();
                }
            }else if( strtype == "logic" ) {
                if( Rf_isMatrix(obj)) {
                    dims[0] = Rcpp::as<Rcpp::LogicalMatrix>(obj).rows();
                    dims[1] = Rcpp::as<Rcpp::LogicalMatrix>(obj).cols();
                } else {
                    dims[0] = 1;
                    dims[1] = Rcpp::as<Rcpp::LogicalVector>(obj).length();
                }
            } else if( strtype == "char" ){
                if( Rf_isMatrix(obj)) {
                    dims[0] = Rcpp::as<Rcpp::CharacterMatrix>(obj).rows();
                    dims[1] = Rcpp::as<Rcpp::CharacterMatrix>(obj).cols();
                } else {
                    dims[0] = 1;
                    dims[1] = Rcpp::as<Rcpp::CharacterVector>(obj).length();
                }
            } else if(strtype == "dataframe"){
                dims[0] = Rcpp::as<Rcpp::DataFrame>(obj).nrows();
                dims[1] = Rcpp::as<Rcpp::DataFrame>(obj).length();
                
            } else {
                dims[0] = 0;
                dims[1] = 0;
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getObjectDims: "<<ex.what()<< " \n";
        }
        
        return(dims);
    }
    
    
    
    /**
     * @brief Calculates optimal block size for matrix operations
     *
     * @param nRowsA Number of rows in matrix A
     * @param nColsA Number of columns in matrix A
     * @param nRowsB Number of rows in matrix B
     * @param nColsB Number of columns in matrix B
     * @param ifactor Block size scaling factor
     * @param block_size Optional user-specified block size
     * @return int Optimal block size
     *
     * @details Block size determination:
     * - Considers matrix dimensions
     * - Respects maximum block size limits
     * - Allows user override with warnings
     * - Optimizes for memory usage
     *
     * @throws std::exception on calculation errors
     */
    inline int getMaxBlockSize ( int nRowsA, int nColsA, int nRowsB, int nColsB, int ifactor, Rcpp::Nullable<int> block_size = R_NilValue) 
    {
        int iblock_size;
        
        try
        {
            // Get system information for intelligent optimization
            const SystemInfo& sys_info = get_system_info();
            
            iblock_size = std::min( std::min( nRowsA, nColsA), std::min( nRowsB, nColsB) );

            if (block_size.isNotNull()) {
                iblock_size = Rcpp::as<int> (block_size);
                if( (unsigned)iblock_size > (MAXBLOCKSIZE / ifactor) ) {
                    Rcpp::warning("Warning: block size %i is bigger than the maximum recommended %i.", iblock_size, (MAXBLOCKSIZE / ifactor));
                }
            } else {
                // System-aware block size optimization
                if ((unsigned)iblock_size > (MAXBLOCKSIZE / ifactor)) {
                    iblock_size = MAXBLOCKSIZE / ifactor;
                }
                
                // Apply system-specific optimizations
                if (sys_info.is_hpc) {
                    // For HPC, use slightly smaller blocks for better cache efficiency
                    iblock_size = std::min(iblock_size, (int)(MAXBLOCKSIZE / (ifactor * 1.2)));
                } else if (sys_info.is_linux_server && sys_info.num_cores > 64) {
                    // For large servers, optimize for NUMA
                    iblock_size = std::min(iblock_size, (int)(MAXBLOCKSIZE / (ifactor * 1.1)));
                }
            }
                
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getMaxBlockSize: "<<ex.what()<< " \n";
        }
        
        return(iblock_size);
    }
    
    
    /**
     * @brief Optimizes block size for edge cases
     *
     * @param fullSize Total size of dimension
     * @param blockSize Current block size
     * @param iDesp Current displacement
     * @param currentSize Current size being processed
     * @return size_t Optimized block size
     *
     * @details Optimization strategy:
     * - Handles last block to avoid single element operations
     * - Adjusts for matrix boundaries
     * - Prevents overflow
     * - Optimizes for performance
     *
     * @throws std::exception on optimization errors
     */
    inline size_t getOptimBlockSize( size_t fullSize, size_t blockSize, size_t iDesp, size_t currentSize ) 
    {
        try
        {
            if( iDesp + blockSize == fullSize - 1) {
                currentSize = blockSize + 1;
            } else if( iDesp + blockSize > fullSize ) { 
                currentSize = fullSize - iDesp; 
            } else {
                currentSize = blockSize;
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getOptimBlockSize: "<<ex.what()<< " \n";
        }
        
        return(currentSize);
    }

    
    
    // Get the number of rows to read taking in to account the maximum elements per block
    // util when we have rectangular matrices, especially in omics data where we have
    // few samples and thousands of variables
    
    /**
     * @brief Calculates optimal block size for matrix operations
     *
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @return std::vector<hsize_t> Vector containing optimized dimensions
     *
     * @details Particularly useful for:
     * - Omics data with few samples and many variables
     * - Memory-efficient processing of rectangular matrices
     * - Optimizing I/O operations
     *
     * @note Considers both memory constraints and processing efficiency
     */
    inline std::vector<hsize_t> getMatrixBlockSize( int nrows, int ncols ) 
    {
        size_t  maxRows = nrows,
                maxCols = ncols;
        
        std::vector<hsize_t> blockSize = {0, 0};
        
        try
        {
            // Get system information for optimization
            const SystemInfo& sys_info = get_system_info();
            
            // Base calculation (existing logic)
            if( nrows < ncols ) {
                if( maxRows < MAXBLOCKSIZE ){
                    maxRows = nrows;
                } else{
                    maxRows = MAXBLOCKSIZE;
                }
                
                maxCols = std::floor( MAXELEMSINBLOCK / maxRows );
                if( maxCols> (unsigned)ncols || maxCols + 1 == (unsigned)ncols) {
                    maxCols = ncols;
                }
            } else {
                if( maxCols < MAXBLOCKSIZE ){
                    maxCols = ncols;
                } else{
                    maxCols = MAXBLOCKSIZE;
                }
                maxRows = std::floor( MAXELEMSINBLOCK / maxCols );
                if( maxRows> (unsigned)nrows || maxRows + 1 == (unsigned)nrows) {
                    maxRows = nrows;
                }
            }    

            // Apply system-specific optimizations
            if (sys_info.is_hpc) {
                // HPC systems: optimize for cache efficiency
                maxRows = std::min(maxRows, (size_t)(MAXBLOCKSIZE * 0.8));
                maxCols = std::min(maxCols, (size_t)(MAXBLOCKSIZE * 0.8));
            } else if (sys_info.is_linux_server && sys_info.num_cores > 64) {
                // Large servers: NUMA-aware optimization
                maxRows = std::min(maxRows, (size_t)(MAXBLOCKSIZE * 0.9));
                maxCols = std::min(maxCols, (size_t)(MAXBLOCKSIZE * 0.9));
            }
            
            blockSize[0] = maxRows;
            blockSize[1] = maxCols;
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getMatrixBlockSize: "<<ex.what()<< " \n";
        }
        
        return(blockSize);
    }
    
    
    
    // Get the number of rows to read taking in to account the maximum elements per block
    // util when we have rectangular matrices, especially in omics data where we have
    // few samples and thousands of variables
    /**
     * @brief Calculates optimal block size for vector operations
     *
     * @param maxSize Maximum size of the vector
     * @return hsize_t Optimal block size for vector processing
     *
     * @details Block size optimization:
     * - Considers memory constraints
     * - Optimizes for vector operations
     * - Ensures efficient I/O
     *
     * @note Particularly useful for long vectors that don't fit in memory
     */
    inline hsize_t getVectorBlockSize(int maxSize) 
    {
        hsize_t blockSize = 0;
        
        try
        {
            const SystemInfo& sys_info = get_system_info();
            
            if( (unsigned)maxSize > MAXELEMSINBLOCK) {
                blockSize = MAXELEMSINBLOCK;
            } else {
                blockSize = maxSize;
            }
            
            // Apply system-specific limits
            if (sys_info.is_hpc) {
                blockSize = std::min(blockSize, (hsize_t)(MAXELEMSINBLOCK * 0.8));
            } else if (sys_info.is_linux_server && sys_info.num_cores > 64) {
                blockSize = std::min(blockSize, (hsize_t)(MAXELEMSINBLOCK * 0.9));
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getVectorBlockSize: "<<ex.what()<< " \n";
        }
        
        return(blockSize);
    }
    
    
    
    /**
     * @brief Determines initial position for matrix operations
     *
     * @param transp Whether to use transposed coordinates
     * @param desp Displacement from start
     * @return std::vector<hsize_t> Initial position coordinates
     *
     * @details Position calculation:
     * - Handles both normal and transposed matrices
     * - Accounts for displacement
     * - Returns [row, col] coordinates
     */
    inline std::vector<hsize_t> getInitialPosition(bool transp, int desp)
    {
        std::vector<hsize_t> voffset = {0, 0};
        
        if(transp == true) {
            voffset[1] = desp;
        } else {
            voffset[0] = desp;
        }
        
        return(voffset);
    }
    
    
    
    /**
     * @brief Checks if a file exists at the specified path
     *
     * @param fullPath Complete file path to check
     * @return bool True if file exists, false otherwise
     *
     * @details File check:
     * - Uses system calls to verify existence
     * - Handles both relative and absolute paths
     * - Thread-safe implementation
     */
    inline bool Rcpp_FileExist(std::string fullPath) 
    {
        bool exists = false;
        
        std::fstream fileStream;
        fileStream.open(fullPath);
        
        if (fileStream.good()) {
            exists = true;
        } else {
            exists = false;
        }
        return(exists);
    }
    
    
    
    /**
     * @brief Calculates dimensions for reading matrix blocks
     *
     * @param transp Whether to transpose dimensions
     * @param count Number of elements to read
     * @param rows Total rows in matrix
     * @param cols Total columns in matrix
     * @return std::vector<hsize_t> Dimensions for reading
     *
     * @details Dimension calculation:
     * - Handles transposed matrices
     * - Respects matrix boundaries
     * - Optimizes for memory usage
     * - Returns [rows_to_read, cols_to_read]
     */
    inline std::vector<hsize_t> getSizetoRead(bool transp, int count, int rows, int cols)
    {
        std::vector<hsize_t> vcount = {0, 0};

        if(transp == true)
        {
            vcount[0] = cols;
            vcount[1] = count;
            
        } else {
            vcount[0] = count;
            vcount[1] = cols;
        }
        
        return(vcount);
    }
    
    
    
    /**
     * @brief Intelligent thread determination with system awareness
     */
    inline unsigned int get_threads(bool bparal, Rcpp::Nullable<int> threads = R_NilValue) 
    {
        if(bparal == false) {
            return 1;
        }
        
        // Use the smart thread detection
        return get_optimal_threads(threads);
    }

    
    
    /**
     * @brief Start performance monitoring for an operation
     * @return PerformanceInfo Structure to track performance
     */
    inline PerformanceInfo start_performance_monitoring() {
        PerformanceInfo perf;
        const SystemInfo& sys_info = get_system_info();
        
        perf.execution_time = 0.0;
        perf.threads_used = getDTthreads(INT_MAX, false);
        perf.system_type = sys_info.system_type;
        perf.memory_used = 0;
        perf.was_optimized = sys_info.is_hpc || sys_info.is_linux_server;
        
        return perf;
    }
    
    /**
     * @brief Check if current system should use HPC optimizations
     * @return bool True if HPC optimizations should be used
     */
    inline bool should_use_hpc_optimizations() {
        const SystemInfo& sys_info = get_system_info();
        return sys_info.is_hpc || sys_info.is_linux_server;
    }
    
    /**
     * @brief Get recommended OpenMP schedule for current system
     * @return std::string OpenMP schedule recommendation
     */
    inline std::string get_recommended_omp_schedule() {
        const SystemInfo& sys_info = get_system_info();
        
        if (sys_info.is_hpc || sys_info.is_linux_server) {
            return "static";  // Better for HPC/servers
        } else {
            return "dynamic"; // Better for desktops
        }
    }
    
    /**
     * @brief Apply system-optimized OpenMP pragma settings
     * @param num_threads Number of threads to use
     */
    inline void apply_omp_settings(int num_threads) {
        #ifdef _OPENMP
                apply_hpc_optimized_settings(num_threads);
        #endif
    }
    
    /**
     * @brief Test OpenMP performance with current settings
     * @return double Time taken for a simple parallel operation
     */
    inline double test_openmp_performance() {
        const int test_size = 10000;
        std::vector<double> test_data(test_size, 1.0);
        double sum = 0.0;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        #ifdef _OPENMP
                #pragma omp parallel for reduction(+:sum)
                for (int i = 0; i < test_size; ++i) {
                    sum += test_data[i] * test_data[i];
                }
        #else
                for (int i = 0; i < test_size; ++i) {
                    sum += test_data[i] * test_data[i];
                }
        #endif
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        return duration.count() / 1000.0; // Return milliseconds
    }

    /**
     * @brief Print system diagnostic information
     */
    inline void print_system_diagnostics() {
        const SystemInfo& sys_info = get_system_info();
        
        Rcpp::Rcout << "\n=== BigDataStatMeth System Diagnostics ===\n";
        Rcpp::Rcout << "System Type: " << sys_info.system_type << "\n";
        Rcpp::Rcout << "OS: " << sys_info.os_name << "\n";
        Rcpp::Rcout << "Cores: " << sys_info.num_cores << "\n";
        Rcpp::Rcout << "Detection: " << sys_info.detection_method << "\n";
        Rcpp::Rcout << "Conservative Mode: " << (sys_info.is_conservative_mode ? "YES" : "NO") << "\n";
        Rcpp::Rcout << "HPC Environment: " << (sys_info.is_hpc ? "YES" : "NO") << "\n";
        Rcpp::Rcout << "Linux Server: " << (sys_info.is_linux_server ? "YES" : "NO") << "\n";
        Rcpp::Rcout << "Recommended Threads: " << sys_info.recommended_threads << "\n";
        Rcpp::Rcout << "Current Threads: " << getDTthreads(INT_MAX, false) << "\n";
        Rcpp::Rcout << "Recommended Schedule: " << get_recommended_omp_schedule() << "\n";
        
        if (sys_info.is_hpc) {
            Rcpp::Rcout << "\nHPC OPTIMIZATION: Smart settings applied\n";
            Rcpp::Rcout << "- Uses all allocated threads\n";
            Rcpp::Rcout << "- Static scheduling\n";
            Rcpp::Rcout << "- Optimized block sizes\n";
            
            // Show job info if available
            const char* job_id = std::getenv("SLURM_JOB_ID");
            if (!job_id) job_id = std::getenv("PBS_JOBID");
            if (job_id) {
                Rcpp::Rcout << "- Job ID: " << job_id << "\n";
            }
            
            const char* slurm_cpus = std::getenv("SLURM_CPUS_PER_TASK");
            const char* pbs_cpus = std::getenv("PBS_NCPUS");
            if (slurm_cpus) {
                Rcpp::Rcout << "- SLURM allocation: " << slurm_cpus << " CPUs\n";
            }
            if (pbs_cpus) {
                Rcpp::Rcout << "- PBS allocation: " << pbs_cpus << " CPUs\n";
            }
            
        } else if (sys_info.is_linux_server) {
            Rcpp::Rcout << "\nLINUX SERVER OPTIMIZATION: Smart settings applied\n";
            Rcpp::Rcout << "- Uses all available threads\n";
            Rcpp::Rcout << "- Static scheduling\n";
            Rcpp::Rcout << "- NUMA-aware optimizations\n";
        }
        
        Rcpp::Rcout << "==========================================\n\n";
    }

}

#endif // BIGDATASTATMETH_UTILITIES_HPP