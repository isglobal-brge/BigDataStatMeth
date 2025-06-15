/**
 * @file BigDataStatMeth.hpp
 * @brief Main header file for the BigDataStatMeth package
 * @details This is the main header file for the BigDataStatMeth package, which provides
 * statistical methods and utilities for handling large-scale data using HDF5 and in-memory
 * operations. The package includes:
 * 
 * Key components:
 * - HDF5-based matrix operations
 * - Memory-efficient algorithms
 * - Parallel processing capabilities
 * - Statistical computations
 * - Omics data utilities
 * 
 * The package is organized into several modules:
 * - hdf5Algebra: Matrix operations on HDF5 files
 * - hdf5Utilities: HDF5 file handling utilities
 * - hdf5Omics: Specialized functions for omics data
 * - memAlgebra: In-memory matrix operations
 * - Utilities: General utility functions
 * 
 * @author Dolors Pelegrí
 * @version 0.63.0
 * @date 2025
 */



#ifndef BIGDATASTATMETH_HPP
#define BIGDATASTATMETH_HPP

// Version and feature macros
#define BIGDATASTATMETH_VERSION_MAJOR 1
#define BIGDATASTATMETH_VERSION_MINOR 1
#define BIGDATASTATMETH_VERSION_PATCH 0
#define BIGDATASTATMETH_VERSION_STRING "1.1.0"

#define BIGDATASTATMETH_HAS_IO_OPTIMIZATION 1
#define BIGDATASTATMETH_HAS_STORAGE_DETECTION 1
#define BIGDATASTATMETH_HAS_SYSTEM_ANALYSIS 1
#define BIGDATASTATMETH_HAS_PERFORMANCE_DIAGNOSTICS 1
#define BIGDATASTATMETH_HAS_INTEGRATION_HELPERS 1

// =============================================================================
// SYSTEM INCLUDES 
// =============================================================================

#include <RcppEigen.h>

// Common headers
#include "H5Cpp.h" //1.8.18
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <regex>
#include <chrono>
#include <map>
#include <utility>
#include <cmath>
#include <iomanip>

// Common headers - multi-threading
#include <thread>
#include "Utilities/pkg_omp.h" // first for clang-13-omp, #5122

#include <algorithm>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#define BIGDATASTATMETH_HAS_OPENMP 1
#else
#define BIGDATASTATMETH_HAS_OPENMP 0
#endif


#include <R.h>
#include <Rversion.h>
#if !defined(R_VERSION) || R_VERSION < R_Version(3, 5, 0)  // R-exts$6.14
#  define ALTREP(x) 0     // #2866
#  define USE_RINTERNALS  // #3301
#  define DATAPTR_RO(x) ((const void *)DATAPTR(x))
#endif
#include <Rinternals.h>
#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped STRING_ELT and VECTOR_ELT
#include <stdint.h>    // for uint64_t rather than unsigned long long
#include <stdbool.h>
#ifdef WIN32  // positional specifiers (%n$) used in translations; #4402
#  define snprintf dt_win_snprintf  // see our snprintf.c; tried and failed to link to _sprintf_p on Windows
#endif
#ifdef sprintf
#undef sprintf
#endif

#define IS_TRUE_OR_FALSE(x) (TYPEOF(x)==LGLSXP && LENGTH(x)==1 && LOGICAL(x)[0]!=NA_LOGICAL)

// for use with bit64::integer64
#define NA_INTEGER64  INT64_MIN
#define MAX_INTEGER64 INT64_MAX

// openme-utils.cpp functions
void initDTthreads();
int getDTthreads(const int64_t n, const bool throttle);
static const char *mygetenv(const char *name, const char *unset);
SEXP getDTthreads_R(SEXP verbose);


// =============================================================================
// CONSTANTS DEFINITION
// =============================================================================

/**
 * @brief Maximum length for names in HDF5 objects
 */
#define MAX_NAME 1024

const int RANK1 = 1;
const int RANK2 = 2;
const int RANK3 = 3;
const int DIM1 = 1;
const int DIM2 = 2;
const int DIM3 = 3;

const int MAXSTRING = 20;
const hsize_t MAXSTRBLOCK = 1 << 5;

/**
 * @brief Maximum number of elements in a block for HDF5 operations
 * @details This constant defines the maximum number of elements that can be processed
 * in a single block during HDF5 operations. The value is set to (2^29 - 1) to ensure
 * efficient memory usage while maintaining performance.
 */
const hsize_t MAXELEMSINBLOCK = ((2 << 29) - 1);

/**
 * @brief Maximum block size for matrix operations
 * @details Calculated as the square root of MAXELEMSINBLOCK to optimize matrix
 * block operations.
 */
const hsize_t MAXBLOCKSIZE = std::floor(std::sqrt(MAXELEMSINBLOCK));

/**
 * @brief Maximum block size for multiplication operations
 */
const hsize_t MAXMULTBLOCKSIZE = 1 << 13;

/**
 * @brief Execution status codes
 * @{
 */
const int EXEC_OK = 0;      /**< Successful execution */
const int EXEC_ERROR = 1;   /**< Execution error */
const int EXEC_WARNING = 2; /**< Execution warning */
/** @} */

namespace BigDataStatMeth {

    /**
     * @brief Storage type enumeration for precise hardware detection
     */
    enum class StorageType {
            HDD = 1,        ///< Traditional hard disk drive
            SATA_SSD = 4,   ///< SATA solid state drive  
            NVME_SSD = 8,   ///< NVMe solid state drive
            NETWORK = 6,    ///< Network/distributed storage
            UNKNOWN = 2     ///< Unknown storage type
    };
    
    /**
     * @brief System type enumeration for environment detection
     */
    enum class SystemType {
        DESKTOP,        ///< Desktop/laptop system
        SERVER,         ///< Server system (16+ cores)
        HPC_CLUSTER,    ///< HPC cluster with job scheduler
        CRAN_CHECK,     ///< CRAN checking environment
        CONTAINER       ///< Containerized environment
    };
    
    /**
     * @brief Performance test results structure
     */
    struct PerformanceResult {
        std::string operation_type;
        int test_size;
        int iterations;
        int threads_used;
        std::vector<double> times_ms;
        double mean_time_ms;
        double min_time_ms;
        double max_time_ms;
        double std_dev_ms;
        double ops_per_second;
        double relative_std_dev;
        
        PerformanceResult() : test_size(0), iterations(0), threads_used(0), 
        mean_time_ms(0), min_time_ms(0), max_time_ms(0), 
        std_dev_ms(0), ops_per_second(0), relative_std_dev(0) {}
    };
    
    /**
     * @brief System information structure
     */
    struct SystemInfo {
        std::string version;
        std::string system_type;
        std::string storage_type;
        int cpu_cores;
        int current_openmp_threads;
        int storage_io_capacity;
        int memory_safe_threads_5kx5k;
        double available_memory_gb;
        std::vector<std::string> features;
        bool openmp_enabled;
        bool is_hpc;
        bool is_cran;
        bool is_container;
    };
    
    /**
     * @brief Configuration information structure
     */
    struct ConfigInfo {
        std::string version;
        int current_threads;
        int available_cores;
        std::string system_type;
        std::string storage_type;
        bool is_hpc;
        bool is_cran;
        bool is_container;
        bool optimization_active;
        double thread_efficiency;
        double available_memory_kb;
        std::vector<std::string> environment_vars;
    };
    
    // =============================================================================
    // FORWARD DECLARATIONS AND EARLY UTILITY FUNCTIONS
    // =============================================================================
    
    /**
     * @brief Get BigDataStatMeth version information
     * 
     * @return String with version information
     */
    inline std::string get_version() {
        return BIGDATASTATMETH_VERSION_STRING;
    }
    
    /**
     * @brief Check if specific features are available
     * 
     * @param feature Feature name to check
     * @return true if feature is available, false otherwise
     */
    inline bool has_feature(const std::string& feature) {
        if (feature == "openmp") {
            return BIGDATASTATMETH_HAS_OPENMP;
        } else if (feature == "io_optimization") {
            return BIGDATASTATMETH_HAS_IO_OPTIMIZATION;
        } else if (feature == "storage_detection") {
            return BIGDATASTATMETH_HAS_STORAGE_DETECTION;
        } else if (feature == "system_analysis") {
            return BIGDATASTATMETH_HAS_SYSTEM_ANALYSIS;
        } else if (feature == "performance_diagnostics") {
            return BIGDATASTATMETH_HAS_PERFORMANCE_DIAGNOSTICS;
        } else if (feature == "integration_helpers") {
            return BIGDATASTATMETH_HAS_INTEGRATION_HELPERS;
        }
        return false;
    }
    
    // =============================================================================
    // FORWARD DECLARATIONS FOR FUNCTIONS USED IN SYSTEM ANALYSIS
    // =============================================================================
    
    // Core system detection functions
    StorageType detect_storage_type_advanced(bool enable_benchmark = false, size_t benchmark_size = 1024);
    SystemType detect_system_type();
    int calculate_memory_safe_threads(size_t matrix_elements, size_t element_size = 8, double safety_factor = 0.8);
    SystemInfo get_system_info_cpp(bool include_benchmark = false);
    ConfigInfo get_config_info_cpp();
    
    // I/O optimization core functions
    // int get_io_aware_threads(SEXP user_threads, int matrix_rows, int matrix_cols, const std::string& operation_type);
    int get_io_aware_threads(Rcpp::Nullable<int> user_threads, int matrix_rows, int matrix_cols, 
                             const std::string& operation_type = "multiplication");
    // int get_io_aware_threads_enhanced(SEXP user_threads, int matrix_rows, int matrix_cols, 
    //                                   const std::string& operation_type, bool enable_storage_benchmark);
    int get_io_aware_threads_enhanced(Rcpp::Nullable<int> user_threads, int matrix_rows, int matrix_cols,
                                      const std::string& operation_type, bool enable_storage_benchmark);
    void configure_bigdata_io_threads(int matrix_rows, int matrix_cols, const std::string& operation_type, 
                                      bool enable_storage_benchmark = false);
    
    // I/O intensity calculation functions
    double calculate_io_intensity_ratio(size_t matrix_elements, const std::string& operation_type);
    double calculate_io_intensity_ratio_enhanced(size_t matrix_elements, const std::string& operation_type, 
                                                 StorageType storage_type);
    
    // Diagnostic functions
    void print_io_thread_diagnosis(int matrix_rows, int matrix_cols, const std::string& operation_type);
    void diagnose_io_cpp(int matrix_rows, int matrix_cols, const std::string& operation_type = "multiplication",
                         bool enable_storage_benchmark = false, bool print_output = true);
    
    // Performance testing functions
    PerformanceResult test_performance_cpp(int test_size = 2000, const std::string& operation_type = "multiplication",
                                           int iterations = 3);
    std::pair<int, std::vector<PerformanceResult>> auto_tune_cpp(int sample_rows = 1000, int sample_cols = 1000,
                                                                 const std::string& operation_type = "multiplication",
                                                                 int max_threads = 0, int test_iterations = 2);
    std::map<std::string, bool> validate_optimization_cpp(bool quick_test = false, bool include_benchmarks = false);
    ConfigInfo check_config_cpp();
    std::pair<PerformanceResult, PerformanceResult> compare_performance_cpp(int matrix_rows = 2000, int matrix_cols = 2000,
                                                                            const std::string& operation_type = "multiplication",
                                                                            int iterations = 3);
    
    // Integration helper functions
    int auto_configure_threads(int matrix_rows, int matrix_cols, const std::string& operation_type, int user_threads = 0);
    // int integrate_io_optimization(int matrix_rows, int matrix_cols, const std::string& operation_type,
    //                               int user_threads = 0, bool verbose = false);
    int integrate_io_optimization(int matrix_rows, int matrix_cols, const std::string& operation_type,
                                  Rcpp::Nullable<int> user_threads = R_NilValue, bool verbose = false);
    int optimize_hdf5_operation(const std::string& file_path, int matrix_rows, int matrix_cols,
                                const std::string& operation_type, int user_threads = 0);
    int monitor_memory_pressure(size_t matrix_elements, int current_threads);
    int adaptive_thread_adjustment(std::chrono::high_resolution_clock::time_point operation_start_time,
                                   int initial_threads, double target_efficiency = 0.8);
    
    // Performance logging functions
    void enable_performance_logging(const std::string& log_file = "bigdatastatmeth_performance.log");
    void disable_performance_logging();

} // namespace BigDataStatMeth (temporary close for includes)

    // =============================================================================
    // MAIN SYS UTILITIES (now all forward declarations are available)
    // =============================================================================
    #include "Utilities/performance/sysAnalysis.hpp"
    #include "Utilities/performance/sysIODiagnostic.hpp"
    #include "Utilities/performance/sysPerformance.hpp"
    #include "Utilities/performance/sysOmpIO.hpp"
    // #include "Utilities/performance/io_aware.hpp"
    // #include "Utilities/performance/sysOmp.hpp"
    
    
    // =============================================================================
    // MAIN UTILITIES (HDF5 AND ON-MEMORY)
    // =============================================================================
    #include "Utilities/openme-utils.hpp"
    #include "Utilities/Utilities.hpp"
    #include "hdf5Utilities/hdf5Utilities.hpp"
    #include "hdf5Omics/hdf5OmicsUtils.hpp"
    
    
    // =============================================================================
    // CORE UTILITIES - Classes and related to classes
    // =============================================================================
    #include "hdf5Utilities/hdf5Files.hpp"
    #include "hdf5Utilities/hdf5Groups.hpp"
    #include "hdf5Utilities/hdf5Datasets.hpp"
    #include "hdf5Utilities/hdf5DatasetsInternal.hpp"
    #include "hdf5Utilities/hdf5Dims.hpp"
    #include "hdf5Utilities/hdf5CheckClose.hpp"
    
    
    // =============================================================================
    // ON-MEMORY ALBEGRA AND OTHER FUNCTIONS
    // =============================================================================
    #include "memAlgebra/memOtherFunctions.hpp"
    #include "memAlgebra/memOptimizedProducts.hpp"
    #include "memAlgebra/memMultiplication.hpp"
    #include "memAlgebra/memSubstract.hpp"
    #include "memAlgebra/memSum.hpp"
    
    
    // =============================================================================
    // HDF5  UTILITIES
    // =============================================================================
    #include "hdf5Utilities/hdf5Methods.hpp"
    #include "hdf5Utilities/hdf5ImportFiles.hpp"
    #include "hdf5Utilities/hdf5RemoveElements.hpp"
    
    
    
    // =============================================================================
    // HDF5 ALBEGRA
    // =============================================================================
    #include "hdf5Algebra/matrixTriangular.hpp"
    #include "hdf5Algebra/matrixDiagonal.hpp"
    #include "hdf5Algebra/vectormatrix.hpp"
    #include "hdf5Algebra/matrixSdMean.hpp"
    #include "hdf5Algebra/matrixNormalization.hpp"
    #include "hdf5Algebra/multiplication.hpp"
    #include "hdf5Algebra/matrixEquationSolver.hpp"
    #include "hdf5Algebra/matrixSvdBlock.hpp"
    #include "hdf5Algebra/matrixSvd.hpp"
    #include "hdf5Algebra/matrixPCA.hpp"
    #include "hdf5Algebra/multiplicationSparse.hpp"
    #include "hdf5Algebra/matrixQR.hpp"
    #include "hdf5Algebra/tcrossprod.hpp"
    #include "hdf5Algebra/matrixInvCholesky.hpp"
    #include "hdf5Algebra/matrixSum.hpp"
    #include "hdf5Algebra/matrixSubstract.hpp"
    #include "hdf5Algebra/matrixPseudoinverse.hpp"
    #include "hdf5Algebra/crossprod.hpp"
    
    
    // =============================================================================
    // HDF5 DATASET UTILITIES
    // =============================================================================
    #include "hdf5Utilities/hdf5SplitDataset.hpp"
    #include "hdf5Utilities/hdf5BindDatasets.hpp"
    #include "hdf5Utilities/hdf5ReduceDataset.hpp"
    #include "hdf5Utilities/hdf5SortDataset.hpp"
    #include "hdf5Utilities/hdf5ApplytoDatasets.hpp"
    
    
    // =============================================================================
    // HDF5 OMIC UTILITIES
    // =============================================================================
    #include "hdf5Omics/hdf5RemoveMAF.hpp"
    #include "hdf5Utilities/hdf5ImputeData.hpp"
    #include "hdf5Utilities/hdf5RemoveLowData.hpp"


    namespace BigDataStatMeth {
    
        /**
         * @brief Quick setup function for common use cases
         * 
         * One-line setup that applies optimal I/O-aware configuration
         * for most common usage scenarios.
         * 
         * @param matrix_rows Typical matrix rows in your workload
         * @param matrix_cols Typical matrix columns in your workload
         * @param operation_type Primary operation type
         * @param enable_diagnostics Whether to print diagnostic information
         * @param enable_storage_benchmark Whether to run storage benchmark
         */
        inline void quick_setup(int matrix_rows, int matrix_cols,
                                const std::string& operation_type = "multiplication",
                                bool enable_diagnostics = false,
                                bool enable_storage_benchmark = false) {
            
            if (enable_diagnostics) {
                print_io_thread_diagnosis(matrix_rows, matrix_cols, operation_type);
            }
            
            configure_bigdata_io_threads(matrix_rows, matrix_cols, operation_type, enable_storage_benchmark);
        }
    
    } // namespace BigDataStatMeth
 
#endif // BIGDATASTATMETH_HPP