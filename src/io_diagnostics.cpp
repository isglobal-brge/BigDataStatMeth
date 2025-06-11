/**
 * @file io_diagnostics.cpp
 * @brief R-accessible diagnostic functions for I/O-aware optimization
 * 
 * Simple functions that expose the I/O-aware optimization capabilities
 * to R users without complex dependencies.
 */

#include <BigDataStatMeth.hpp>

//' Diagnose I/O performance for BigDataStatMeth operations
//'
//' Analyzes your matrix size and operation to show the I/O-aware 
//' thread optimization strategy and expected performance impact.
//' This function provides detailed system analysis including storage
//' detection, I/O intensity calculation, and optimization reasoning.
//'
//' @param matrix_rows Number of rows in your matrix
//' @param matrix_cols Number of columns in your matrix
//' @param operation_type Type of BigDataStatMeth operation
//' @details
//' Supported operation types:
//' \itemize{
//'   \item "multiplication" - Matrix multiplication (A %*% B)
//'   \item "svd" - Singular Value Decomposition
//'   \item "crossprod" - Cross-product (t(A) %*% A)
//'   \item "tcrossprod" - Transpose cross-product (A %*% t(A))
//' }
//' 
//' The analysis includes:
//' \itemize{
//'   \item Storage type detection (HDD/SSD/NVMe/Network)
//'   \item I/O intensity ratio calculation
//'   \item System type identification (HPC/Server/Desktop)
//'   \item Thread recommendation with reasoning
//'   \item Expected performance improvement estimates
//' }
//' @examples
//' # Diagnose your specific workload
//' bdDiagnoseIO(5000, 5000, "multiplication")
//' 
//' # Check different operations
//' bdDiagnoseIO(10000, 8000, "svd")
//' bdDiagnoseIO(5000, 5000, "crossprod")
//' @export
// [[Rcpp::export]]
void bdDiagnoseIO(int matrix_rows, int matrix_cols, 
                  std::string operation_type = "multiplication") {
    try {
        BigDataStatMeth::print_io_thread_diagnosis(matrix_rows, matrix_cols, operation_type);
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdDiagnoseIO: " << ex.what() << std::endl;
    }
}

//' Get I/O-aware thread recommendation
//'
//' Returns the optimal number of threads for your specific matrix size
//' and operation, considering I/O bottlenecks, storage characteristics,
//' and system constraints. This function performs the same analysis as
//' the internal optimization but returns the result without applying it.
//'
//' @param matrix_rows Number of rows in your matrix
//' @param matrix_cols Number of columns in your matrix
//' @param operation_type Type of operation (default: "multiplication")
//' @param max_threads Maximum threads to consider (optional)
//' @return Integer with recommended thread count
//' @details
//' The recommendation considers:
//' \itemize{
//'   \item Storage I/O capacity (1 for HDD, 4 for SSD, 8 for NVMe)
//'   \item I/O to computation ratio for the specific operation
//'   \item System type (HPC cluster, server, desktop)
//'   \item Memory constraints on large systems (>64 cores)
//'   \item Job scheduler allocations (SLURM_CPUS_PER_TASK, PBS_NCPUS)
//'   \item CRAN compliance requirements
//' }
//' @note This function does not modify any OpenMP settings,
//'       use bdConfigureIO() to apply the recommendations
//' @examples
//' # Get recommendation for your workload
//' threads <- bdGetOptimalThreads(5000, 5000, "multiplication")
//' cat("Recommended threads:", threads, "\n")
//' 
//' # Limit maximum threads
//' threads <- bdGetOptimalThreads(10000, 10000, "svd", max_threads = 16)
//' 
//' # Use recommendation in manual configuration
//' # omp_set_num_threads(threads)
//' @export
// [[Rcpp::export]]
int bdGetOptimalThreads(int matrix_rows, int matrix_cols,
                       std::string operation_type = "multiplication",
                       Rcpp::Nullable<int> max_threads = R_NilValue) {
    try {
        return BigDataStatMeth::get_io_aware_threads(max_threads, matrix_rows, matrix_cols, operation_type);
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdGetOptimalThreads: " << ex.what() << std::endl;
        return 1;
    }
}

//' Configure BigDataStatMeth with I/O-aware optimization
//'
//' Applies optimal OpenMP configuration for your typical workload.
//' This sets the thread count and OpenMP parameters for the current R session
//' based on intelligent analysis of your matrix size, operation type,
//' storage characteristics, and system capabilities.
//'
//' @param matrix_rows Typical number of rows in your matrices
//' @param matrix_cols Typical number of columns in your matrices  
//' @param operation_type Primary operation type you'll be using (default: "multiplication")
//' @details
//' This function:
//' \itemize{
//'   \item Calculates optimal thread count using I/O-aware analysis
//'   \item Applies OpenMP thread configuration (omp_set_num_threads)
//'   \item Sets appropriate OpenMP scheduling policy
//'   \item Disables dynamic thread adjustment for consistent performance
//'   \item Configures thread affinity settings when beneficial
//' }
//' 
//' The configuration persists for the current R session and affects
//' all subsequent BigDataStatMeth operations and other OpenMP code.
//' @note Call this function once at the beginning of your analysis session
//' @warning This modifies global OpenMP settings which may affect
//'          other packages using OpenMP in the same R session
//' @examples
//' # Configure for matrix multiplication workload
//' bdConfigureIO(5000, 5000, "multiplication")
//' 
//' # Configure for SVD operations
//' bdConfigureIO(10000, 8000, "svd")
//' 
//' # Now your BigDataStatMeth operations will use optimal settings
//' result <- bdblockmult_hdf5("data.h5", "matrices", "A", "B")
//' @export
// [[Rcpp::export]]
void bdConfigureIO(int matrix_rows, int matrix_cols,
                   std::string operation_type = "multiplication") {
    try {
        BigDataStatMeth::configure_bigdata_io_threads(matrix_rows, matrix_cols, operation_type);
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdConfigureIO: " << ex.what() << std::endl;
    }
}

//' Check current BigDataStatMeth thread configuration
//'
//' Shows the current OpenMP thread configuration and system information.
//' Useful for verifying that optimization has been applied correctly,
//' understanding system characteristics, and debugging performance issues.
//'
//' @return List with current configuration details
//' @details
//' Returns a list containing:
//' \itemize{
//'   \item current_threads - Currently configured OpenMP threads
//'   \item available_cores - Total CPU cores available
//'   \item system_type - Detected system type (HPC/Server/Desktop/CRAN)
//'   \item is_hpc - Whether HPC environment was detected
//'   \item is_cran - Whether CRAN compliance mode is active
//'   \item optimization_active - Whether I/O-aware optimization appears active
//'   \item thread_efficiency - Ratio of current threads to available cores
//' }
//' @note This function only reads current settings and does not
//'       modify any configuration
//' @examples
//' # Check current configuration
//' config <- bdCheckConfig()
//' print(config)
//' 
//' # See if optimization is active
//' if (config$optimization_active) {
//'   message("I/O-aware optimization is active")
//' } else {
//'   message("Consider running bdConfigureIO() for optimization")
//' }
//' 
//' # Check thread efficiency
//' if (config$thread_efficiency < 0.5) {
//'   message("Conservative threading detected - likely I/O optimized")
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List bdCheckConfig() {
    try {
        int current_threads = 1;
        int available_cores = 1;
        
#ifdef _OPENMP
        current_threads = omp_get_max_threads();
        available_cores = omp_get_num_procs();
#endif
        
        // Check environment indicators
        bool is_hpc = (std::getenv("SLURM_JOB_ID") != nullptr) ||
                      (std::getenv("PBS_JOBID") != nullptr);
        bool is_cran = (std::getenv("OMP_THREAD_LIMIT") != nullptr) &&
                       (std::atoi(std::getenv("OMP_THREAD_LIMIT")) <= 2);
        
        // Determine if optimization appears active
        bool optimization_active = (current_threads < available_cores) || 
                                  (current_threads <= 32 && available_cores > 64);
        
        std::string system_type = "Desktop";
        if (is_hpc) {
            system_type = "HPC Cluster";
        } else if (available_cores > 16) {
            system_type = "Server";
        }
        
        if (is_cran) {
            system_type = "CRAN Environment";
        }
        
        return Rcpp::List::create(
            Rcpp::Named("current_threads") = current_threads,
            Rcpp::Named("available_cores") = available_cores,
            Rcpp::Named("system_type") = system_type,
            Rcpp::Named("is_hpc") = is_hpc,
            Rcpp::Named("is_cran") = is_cran,
            Rcpp::Named("optimization_active") = optimization_active,
            Rcpp::Named("thread_efficiency") = (double)current_threads / available_cores
        );
        
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdCheckConfig: " << ex.what() << std::endl;
        return R_NilValue;
    }
}

//' Quick performance test for current thread configuration
//'
//' Performs a simple matrix operation to test the performance
//' of the current thread configuration. Useful for measuring
//' the impact of I/O-aware optimization and comparing different
//' configurations.
//'
//' @param test_size Size of test matrix (default: 2000)
//' @return Numeric vector with timing in milliseconds
//' @details
//' The test performs a matrix multiplication of two random matrices
//' of size test_size × test_size. This operation uses the current
//' OpenMP configuration and provides a benchmark for computational
//' performance under the current thread settings.
//' 
//' The test is CPU-intensive and may not reflect I/O-bound BigDataStatMeth
//' operations, but provides a baseline for thread configuration effectiveness.
//' @note 
//' - Larger test sizes provide more stable timing but take longer
//' - Results may vary between runs due to system load
//' - This test uses in-memory matrices, not HDF5 I/O
//' @examples
//' # Test current performance
//' time_ms <- bdTestPerformance()
//' cat("Test completed in", time_ms, "milliseconds\n")
//' 
//' # Compare different configurations
//' time_before <- bdTestPerformance()
//' bdConfigureIO(5000, 5000, "multiplication")
//' time_after <- bdTestPerformance()
//' 
//' # Compare results
//' improvement <- time_before / time_after
//' cat("Performance improvement:", improvement, "x\n")
//' 
//' # Test with larger matrices for more stable results
//' time_stable <- bdTestPerformance(test_size = 3000)
//' @export
// [[Rcpp::export]]
double bdTestPerformance(int test_size = 2000) {
    try {
        auto start = std::chrono::high_resolution_clock::now();
        
        // Simple matrix operation that uses OpenMP
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(test_size, test_size);
        Eigen::MatrixXd B = Eigen::MatrixXd::Random(test_size, test_size);
        Eigen::MatrixXd C = A * B;
        
        // Force computation to complete
        volatile double sum = C.sum();
        (void)sum; // Suppress unused variable warning
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        return (double)duration.count();
        
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdTestPerformance: " << ex.what() << std::endl;
        return -1.0;
    }
}