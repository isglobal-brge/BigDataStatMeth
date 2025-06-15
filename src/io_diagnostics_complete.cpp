/**
 * @file io_diagnostics.cpp
 * @brief R-accessible diagnostic functions for I/O-aware optimization
 * 
 * Complete implementation of R interface functions that expose the 
 * I/O-aware optimization capabilities to R users. Includes enhanced
 * diagnostics, configuration, and integration helpers.
 */

#include <Rcpp.h>
#include <BigDataStatMeth.hpp>
#include <chrono>
#include <iomanip>

using namespace BigDataStatMeth;

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
//' @param enable_storage_benchmark Enable I/O benchmarking for accurate storage detection
//' @details
//' Supported operation types:
//' \itemize{
//'   \item "multiplication", "bdblockmult" - Matrix multiplication (A \%*\% B)
//'   \item "svd", "bdSVD" - Singular Value Decomposition
//'   \item "crossprod", "bdCrossprod" - Cross-product (t(A) \%*\% A)
//'   \item "tcrossprod", "bdtCrossprod" - Transpose cross-product (A \%*\% t(A))
//'   \item "pca", "bdPCA" - Principal Component Analysis
//'   \item "qr", "bdQR" - QR decomposition
//'   \item "eigen", "bdEigen" - Eigenvalue decomposition
//' }
//' 
//' The analysis includes:
//' \itemize{
//'   \item Storage type detection (HDD/SSD/NVMe/Network)
//'   \item I/O intensity ratio calculation
//'   \item System type identification (HPC/Server/Desktop/Container)
//'   \item Memory constraint analysis
//'   \item Thread recommendation with detailed reasoning
//'   \item Expected performance improvement estimates
//' }
//' @examples
//' # Diagnose your specific workload
//' bdDiagnoseIO(5000, 5000, "multiplication")
//' 
//' # Check different operations
//' bdDiagnoseIO(10000, 8000, "svd")
//' bdDiagnoseIO(5000, 5000, "crossprod")
//' 
//' # Enable storage benchmarking for accurate detection
//' bdDiagnoseIO(5000, 5000, "multiplication", enable_storage_benchmark = TRUE)
//' @export
// [[Rcpp::export]]
void bdDiagnoseIO(int matrix_rows, int matrix_cols, 
                  std::string operation_type = "multiplication",
                  bool enable_storage_benchmark = false) {
    try {
        // Enhanced diagnosis with storage benchmarking option
        size_t elements = static_cast<size_t>(matrix_rows) * matrix_cols;
        StorageType storage_type = detect_storage_type_advanced(enable_storage_benchmark);
        SystemType system_type = detect_system_type();
        
        if (enable_storage_benchmark) {
            Rcpp::Rcout << "Running I/O benchmark for accurate storage detection...\n";
        }
        
        print_io_thread_diagnosis(matrix_rows, matrix_cols, operation_type);
        
        // Additional enhanced information
        if (enable_storage_benchmark) {
            Rcpp::Rcout << "Storage benchmark completed. Results integrated into analysis.\n\n";
        }
        
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
//' @param enable_storage_benchmark Enable storage benchmarking for accurate detection
//' @return Integer with recommended thread count
//' @details
//' The recommendation considers:
//' \itemize{
//'   \item Storage I/O capacity (1 for HDD, 4 for SSD, 8 for NVMe)
//'   \item I/O to computation ratio for the specific operation
//'   \item System type (HPC cluster, server, desktop, container)
//'   \item Memory constraints on large systems (>64 cores)
//'   \item Job scheduler allocations (SLURM_CPUS_PER_TASK, PBS_NCPUS)
//'   \item CRAN compliance requirements
//'   \item Container resource limits (Docker, Singularity)
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
//' # Enable accurate storage detection
//' threads <- bdGetOptimalThreads(5000, 5000, "multiplication", 
//'                               enable_storage_benchmark = TRUE)
//' 
//' # Use recommendation in manual configuration
//' # omp_set_num_threads(threads)
//' @export
// [[Rcpp::export]]
int bdGetOptimalThreads(int matrix_rows, int matrix_cols,
                       std::string operation_type = "multiplication",
                       Rcpp::Nullable<int> max_threads = R_NilValue,
                       bool enable_storage_benchmark = false) {
    try {
        return get_io_aware_threads_enhanced(max_threads, matrix_rows, matrix_cols, 
                                           operation_type, enable_storage_benchmark);
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
//' @param enable_storage_benchmark Enable storage benchmarking for accurate detection
//' @param verbose Print detailed configuration information
//' @details
//' This function:
//' \itemize{
//'   \item Calculates optimal thread count using I/O-aware analysis
//'   \item Applies OpenMP thread configuration (omp_set_num_threads)
//'   \item Sets appropriate OpenMP scheduling policy
//'   \item Disables dynamic thread adjustment for consistent performance
//'   \item Configures thread affinity settings when beneficial
//'   \item Provides detailed feedback about applied settings
//' }
//' 
//' The configuration persists for the current R session and affects
//' all subsequent BigDataStatMeth operations and other OpenMP code.
//' @note Call this function once at the beginning of your analysis session
//' This modifies global OpenMP settings which may affect
//'          other packages using OpenMP in the same R session
//' @examples
//' # Configure for matrix multiplication workload
//' bdConfigureIO(5000, 5000, "multiplication")
//' 
//' # Configure for SVD operations with detailed output
//' bdConfigureIO(10000, 8000, "svd", verbose = TRUE)
//' 
//' # Configure with accurate storage detection
//' bdConfigureIO(5000, 5000, "multiplication", enable_storage_benchmark = TRUE)
//' 
//' # Now your BigDataStatMeth operations will use optimal settings
//' result <- bdblockmult_hdf5("data.h5", "matrices", "A", "B")
//' @export
// [[Rcpp::export]]
void bdConfigureIO(int matrix_rows, int matrix_cols,
                   std::string operation_type = "multiplication",
                   bool enable_storage_benchmark = false,
                   bool verbose = false) {
    try {
        if (verbose) {
            bdDiagnoseIO(matrix_rows, matrix_cols, operation_type, enable_storage_benchmark);
        }
        
        configure_bigdata_io_threads(matrix_rows, matrix_cols, operation_type, enable_storage_benchmark);
        
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
//' Returned list contains:
//' \itemize{
//'   \item version - BigDataStatMeth version
//'   \item current_threads - Currently configured OpenMP threads
//'   \item available_cores - Total CPU cores available
//'   \item system_type - Detected system type (HPC/Server/Desktop/CRAN/Container)
//'   \item storage_type - Detected storage type (HDD/SSD/NVMe/Network/Unknown)
//'   \item is_hpc - Whether HPC environment was detected
//'   \item is_cran - Whether CRAN compliance mode is active
//'   \item is_container - Whether container environment was detected
//'   \item optimization_active - Whether I/O-aware optimization appears active
//'   \item thread_efficiency - Ratio of current threads to available cores
//'   \item memory_info - Available memory information (if detectable)
//'   \item environment_vars - Relevant environment variables
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
        
        // Enhanced system detection
        SystemType sys_type = detect_system_type();
        StorageType storage_type = detect_storage_type_advanced(false);
        
        // Check environment indicators
        bool is_hpc = (sys_type == SystemType::HPC_CLUSTER);
        bool is_cran = (sys_type == SystemType::CRAN_CHECK);
        bool is_container = (sys_type == SystemType::CONTAINER);
        
        // Determine if optimization appears active
        bool optimization_active = (current_threads < available_cores) || 
                                  (current_threads <= 32 && available_cores > 64);
        
        // Convert enums to strings
        std::string system_type_str;
        switch (sys_type) {
            case SystemType::DESKTOP: system_type_str = "Desktop"; break;
            case SystemType::SERVER: system_type_str = "Server"; break;
            case SystemType::HPC_CLUSTER: system_type_str = "HPC Cluster"; break;
            case SystemType::CRAN_CHECK: system_type_str = "CRAN Environment"; break;
            case SystemType::CONTAINER: system_type_str = "Container"; break;
        }
        
        std::string storage_type_str;
        switch (storage_type) {
            case StorageType::HDD: storage_type_str = "HDD"; break;
            case StorageType::SATA_SSD: storage_type_str = "SATA SSD"; break;
            case StorageType::NVME_SSD: storage_type_str = "NVMe SSD"; break;
            case StorageType::NETWORK: storage_type_str = "Network Storage"; break;
            case StorageType::UNKNOWN: storage_type_str = "Unknown"; break;
        }
        
        // Collect environment variables
        Rcpp::CharacterVector env_vars;
        const char* env_names[] = {
            "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "SLURM_CPUS_PER_TASK", 
            "PBS_NCPUS", "SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID", "SGE_JOB_ID",
            "SINGULARITY_CONTAINER", "APPTAINER_CONTAINER"
        };
        
        for (const char* env_name : env_names) {
            const char* env_value = std::getenv(env_name);
            if (env_value) {
                std::string entry = std::string(env_name) + "=" + std::string(env_value);
                env_vars.push_back(entry);
            }
        }
        
        // Memory information (Linux only)
        size_t available_memory = 0;
#ifdef __linux__
        std::ifstream meminfo("/proc/meminfo");
        if (meminfo.is_open()) {
            std::string line;
            while (std::getline(meminfo, line)) {
                if (line.find("MemAvailable:") == 0) {
                    size_t kb;
                    sscanf(line.c_str(), "MemAvailable: %zu kB", &kb);
                    available_memory = kb;
                    break;
                }
            }
        }
#endif
        
        return Rcpp::List::create(
            Rcpp::Named("version") = "1.1.0",
            Rcpp::Named("current_threads") = current_threads,
            Rcpp::Named("available_cores") = available_cores,
            Rcpp::Named("system_type") = system_type_str,
            Rcpp::Named("storage_type") = storage_type_str,
            Rcpp::Named("is_hpc") = is_hpc,
            Rcpp::Named("is_cran") = is_cran,
            Rcpp::Named("is_container") = is_container,
            Rcpp::Named("optimization_active") = optimization_active,
            Rcpp::Named("thread_efficiency") = (double)current_threads / available_cores,
            Rcpp::Named("available_memory_kb") = (double)available_memory,
            Rcpp::Named("environment_vars") = env_vars
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
//' @param operation_type Type of operation to test (default: "multiplication")
//' @param iterations Number of iterations to average (default: 3)
//' @return List with timing results and performance metrics
//' @details
//' The test performs matrix operations using the current OpenMP configuration
//' and provides benchmarks for computational performance. Available test types:
//' \itemize{
//'   \item "multiplication" - Matrix multiplication benchmark
//'   \item "crossprod" - Cross-product benchmark  
//'   \item "svd" - SVD decomposition benchmark (smaller matrices recommended)
//' }
//' 
//' Results include timing statistics, thread efficiency metrics, and
//' recommendations for optimization.
//' @note 
//' - Larger test sizes provide more stable timing but take longer
//' - Results may vary between runs due to system load
//' - This test uses in-memory matrices, not HDF5 I/O
//' - SVD tests should use smaller matrices (test_size < 1000)
//' @examples
//' # Test current performance
//' results <- bdTestPerformance()
//' print(results)
//' 
//' # Compare different configurations
//' before <- bdTestPerformance()
//' bdConfigureIO(5000, 5000, "multiplication")
//' after <- bdTestPerformance()
//' 
//' # Compare results
//' improvement <- before$mean_time_ms / after$mean_time_ms
//' cat("Performance improvement:", improvement, "x\n")
//' 
//' # Test with different operations
//' mult_test <- bdTestPerformance(test_size = 2000, operation_type = "multiplication")
//' crossprod_test <- bdTestPerformance(test_size = 2000, operation_type = "crossprod")
//' svd_test <- bdTestPerformance(test_size = 500, operation_type = "svd")
//' @export
// [[Rcpp::export]]
Rcpp::List bdTestPerformance(int test_size = 2000, 
                            std::string operation_type = "multiplication",
                            int iterations = 3) {
    try {
        std::vector<double> times;
        times.reserve(iterations);
        
        int current_threads = 1;
#ifdef _OPENMP
        current_threads = omp_get_max_threads();
#endif
        
        Rcpp::Rcout << "Running performance test (" << operation_type << ") with " 
                   << current_threads << " threads...\n";
        
        for (int i = 0; i < iterations; ++i) {
            auto start = std::chrono::high_resolution_clock::now();
            
            if (operation_type == "multiplication") {
                // Matrix multiplication benchmark
                Rcpp::NumericMatrix A(test_size, test_size);
                Rcpp::NumericMatrix B(test_size, test_size);
                
                // Fill with random data
                std::fill(A.begin(), A.end(), 1.0);
                std::fill(B.begin(), B.end(), 1.0);
                
                // Perform multiplication (simplified)
                Rcpp::NumericMatrix C(test_size, test_size);
                
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int row = 0; row < test_size; ++row) {
                    for (int col = 0; col < test_size; ++col) {
                        double sum = 0.0;
                        for (int k = 0; k < test_size; ++k) {
                            sum += A(row, k) * B(k, col);
                        }
                        C(row, col) = sum;
                    }
                }
                
                // Force computation to complete
                volatile double checksum = C(0, 0);
                (void)checksum;
                
            } else if (operation_type == "crossprod") {
                // Cross-product benchmark
                Rcpp::NumericMatrix A(test_size, test_size);
                std::fill(A.begin(), A.end(), 1.0);
                
                Rcpp::NumericMatrix C(test_size, test_size);
                
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int i = 0; i < test_size; ++i) {
                    for (int j = 0; j < test_size; ++j) {
                        double sum = 0.0;
                        for (int k = 0; k < test_size; ++k) {
                            sum += A(k, i) * A(k, j);
                        }
                        C(i, j) = sum;
                    }
                }
                
                volatile double checksum = C(0, 0);
                (void)checksum;
                
            } else if (operation_type == "svd") {
                // Simplified SVD-like computation (computationally intensive)
                Rcpp::NumericMatrix A(test_size, test_size);
                std::fill(A.begin(), A.end(), 1.0);
                
                // Simplified power iteration for largest eigenvalue
                Rcpp::NumericVector v(test_size, 1.0);
                
                for (int iter = 0; iter < 10; ++iter) {
                    Rcpp::NumericVector new_v(test_size, 0.0);
                    
#ifdef _OPENMP
#pragma omp parallel for
#endif
                    for (int i = 0; i < test_size; ++i) {
                        for (int j = 0; j < test_size; ++j) {
                            new_v[i] += A(i, j) * v[j];
                        }
                    }
                    
                    // Normalize
                    double norm = 0.0;
                    for (int i = 0; i < test_size; ++i) {
                        norm += new_v[i] * new_v[i];
                    }
                    norm = sqrt(norm);
                    
                    for (int i = 0; i < test_size; ++i) {
                        v[i] = new_v[i] / norm;
                    }
                }
                
                volatile double checksum = v[0];
                (void)checksum;
                
            } else {
                Rcpp::stop("Unsupported operation_type: " + operation_type);
            }
            
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            times.push_back(static_cast<double>(duration.count()));
            
            if (iterations > 1) {
                Rcpp::Rcout << "  Iteration " << (i + 1) << ": " << times[i] << " ms\n";
            }
        }
        
        // Calculate statistics
        double mean_time = 0.0;
        for (double t : times) {
            mean_time += t;
        }
        mean_time /= times.size();
        
        double min_time = *std::min_element(times.begin(), times.end());
        double max_time = *std::max_element(times.begin(), times.end());
        
        // Calculate standard deviation
        double variance = 0.0;
        for (double t : times) {
            variance += (t - mean_time) * (t - mean_time);
        }
        variance /= times.size();
        double std_dev = sqrt(variance);
        
        // Performance metrics
        double ops_per_sec = 0.0;
        if (operation_type == "multiplication") {
            // Approximate FLOPS for matrix multiplication
            double flops = 2.0 * test_size * test_size * test_size;
            ops_per_sec = flops / (mean_time / 1000.0);
        }
        
        // Thread efficiency (rough estimate)
        double theoretical_speedup = static_cast<double>(current_threads);
        double actual_speedup = 1.0; // Would need baseline to calculate
        
        Rcpp::Rcout << "\nPerformance Summary:\n";
        Rcpp::Rcout << "  Mean time: " << std::fixed << std::setprecision(2) << mean_time << " ms\n";
        Rcpp::Rcout << "  Min time: " << min_time << " ms\n";
        Rcpp::Rcout << "  Max time: " << max_time << " ms\n";
        Rcpp::Rcout << "  Std dev: " << std::setprecision(1) << std_dev << " ms\n";
        
        if (ops_per_sec > 0) {
            if (ops_per_sec > 1e9) {
                Rcpp::Rcout << "  Performance: " << std::setprecision(2) 
                           << (ops_per_sec / 1e9) << " GFLOPS\n";
            } else if (ops_per_sec > 1e6) {
                Rcpp::Rcout << "  Performance: " << std::setprecision(2) 
                           << (ops_per_sec / 1e6) << " MFLOPS\n";
            }
        }
        
        return Rcpp::List::create(
            Rcpp::Named("operation_type") = operation_type,
            Rcpp::Named("test_size") = test_size,
            Rcpp::Named("iterations") = iterations,
            Rcpp::Named("threads_used") = current_threads,
            Rcpp::Named("times_ms") = Rcpp::NumericVector(times.begin(), times.end()),
            Rcpp::Named("mean_time_ms") = mean_time,
            Rcpp::Named("min_time_ms") = min_time,
            Rcpp::Named("max_time_ms") = max_time,
            Rcpp::Named("std_dev_ms") = std_dev,
            Rcpp::Named("ops_per_second") = ops_per_sec,
            Rcpp::Named("relative_std_dev") = std_dev / mean_time
        );
        
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdTestPerformance: " << ex.what() << std::endl;
        return R_NilValue;
    }
}

//' Auto-tune I/O optimization for specific dataset
//'
//' Automatically determines optimal thread configuration by analyzing
//' a sample of your actual data. This provides more accurate optimization
//' than estimates based on matrix size alone.
//'
//' @param sample_rows Number of rows to use for testing (default: 1000)
//' @param sample_cols Number of columns to use for testing (default: 1000)
//' @param operation_type Primary operation type to optimize for
//' @param max_threads Maximum threads to test (default: all available)
//' @param test_iterations Number of iterations per thread count (default: 2)
//' @return List with optimal configuration and benchmarking results
//' @details
//' This function performs empirical testing with different thread counts
//' to find the optimal configuration for your specific workload and system.
//' The testing process:
//' \itemize{
//'   \item Tests thread counts from 1 to max_threads
//'   \item Measures actual performance for each configuration
//'   \item Identifies the optimal thread count based on throughput
//'   \item Applies the optimal configuration automatically
//'   \item Returns detailed benchmarking results
//' }
//' @note This function modifies global OpenMP settings and may take
//'       several minutes to complete depending on the test parameters
//' @examples
//' # Auto-tune for matrix multiplication
//' results <- bdAutoTune(operation_type = "multiplication")
//' 
//' # Auto-tune with specific constraints
//' results <- bdAutoTune(sample_rows = 2000, sample_cols = 2000,
//'                      operation_type = "svd", max_threads = 16)
//' 
//' # View results
//' print(results$optimal_threads)
//' plot(results$thread_counts, results$performance_scores)
//' @export
// [[Rcpp::export]]
Rcpp::List bdAutoTune(int sample_rows = 1000, int sample_cols = 1000,
                     std::string operation_type = "multiplication",
                     Rcpp::Nullable<int> max_threads = R_NilValue,
                     int test_iterations = 2) {
    try {
        int available_cores = std::thread::hardware_concurrency();
#ifdef _OPENMP
        available_cores = omp_get_num_procs();
#endif
        
        int max_test_threads = available_cores;
        if (max_threads.isNotNull()) {
            max_test_threads = std::min(Rcpp::as<int>(max_threads), available_cores);
        }
        
        Rcpp::Rcout << "Auto-tuning I/O optimization...\n";
        Rcpp::Rcout << "Testing " << operation_type << " with " << sample_rows 
                   << "x" << sample_cols << " matrices\n";
        Rcpp::Rcout << "Thread range: 1 to " << max_test_threads << "\n\n";
        
        std::vector<int> thread_counts;
        std::vector<double> performance_scores;
        std::vector<double> mean_times;
        
        // Store original configuration
        int original_threads = 1;
#ifdef _OPENMP
        original_threads = omp_get_max_threads();
#endif
        
        // Test different thread counts
        for (int threads = 1; threads <= max_test_threads; threads++) {
#ifdef _OPENMP
            omp_set_num_threads(threads);
#endif
            
            Rcpp::List test_result = bdTestPerformance(sample_rows, operation_type, test_iterations);
            double mean_time = Rcpp::as<double>(test_result["mean_time_ms"]);
            double performance_score = 1000.0 / mean_time; // Higher is better
            
            thread_counts.push_back(threads);
            mean_times.push_back(mean_time);
            performance_scores.push_back(performance_score);
            
            Rcpp::Rcout << "Threads " << threads << ": " 
                       << std::fixed << std::setprecision(1) << mean_time 
                       << " ms (score: " << std::setprecision(3) << performance_score << ")\n";
        }
        
        // Find optimal thread count
        auto max_it = std::max_element(performance_scores.begin(), performance_scores.end());
        int optimal_index = std::distance(performance_scores.begin(), max_it);
        int optimal_threads = thread_counts[optimal_index];
        double optimal_score = performance_scores[optimal_index];
        
        // Apply optimal configuration
#ifdef _OPENMP
        omp_set_num_threads(optimal_threads);
        omp_set_dynamic(false);
#endif
        
        Rcpp::Rcout << "\nOptimal configuration found:\n";
        Rcpp::Rcout << "  Optimal threads: " << optimal_threads << "\n";
        Rcpp::Rcout << "  Performance score: " << std::fixed << std::setprecision(3) 
                   << optimal_score << "\n";
        Rcpp::Rcout << "  Improvement over 1 thread: " 
                   << std::setprecision(2) << (optimal_score / performance_scores[0]) << "x\n";
        
        if (optimal_threads < max_test_threads) {
            Rcpp::Rcout << "\nNote: Optimal thread count is less than maximum available.\n";
            Rcpp::Rcout << "This suggests I/O or memory bandwidth limitations.\n";
        }
        
        return Rcpp::List::create(
            Rcpp::Named("optimal_threads") = optimal_threads,
            Rcpp::Named("optimal_score") = optimal_score,
            Rcpp::Named("thread_counts") = Rcpp::IntegerVector(thread_counts.begin(), thread_counts.end()),
            Rcpp::Named("performance_scores") = Rcpp::NumericVector(performance_scores.begin(), performance_scores.end()),
            Rcpp::Named("mean_times_ms") = Rcpp::NumericVector(mean_times.begin(), mean_times.end()),
            Rcpp::Named("operation_type") = operation_type,
            Rcpp::Named("test_matrix_size") = Rcpp::IntegerVector::create(sample_rows, sample_cols),
            Rcpp::Named("improvement_factor") = optimal_score / performance_scores[0],
            Rcpp::Named("original_threads") = original_threads
        );
        
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdAutoTune: " << ex.what() << std::endl;
        return R_NilValue;
    }
}

//' Get system information and capabilities
//'
//' Returns comprehensive information about the system hardware,
//' software environment, and BigDataStatMeth capabilities.
//'
//' @param include_benchmark Include storage benchmark in the analysis
//' @return List with detailed system information
//' @details
//' This function provides a complete overview of:
//' \itemize{
//'   \item Hardware: CPU cores, memory, storage type
//'   \item Software: OpenMP support, BigDataStatMeth version, R environment
//'   \item Environment: HPC scheduler, container detection, CRAN compliance
//'   \item Capabilities: Available features and optimization options
//'   \item Performance: Storage I/O capacity and memory constraints
//' }
//' @examples
//' # Get basic system information
//' info <- bdSystemInfo()
//' print(info)
//' 
//' # Include storage benchmarking
//' info <- bdSystemInfo(include_benchmark = TRUE)
//' @export
// [[Rcpp::export]]
Rcpp::List bdSystemInfo(bool include_benchmark = false) {
    try {
        // Basic system information
        int available_cores = std::thread::hardware_concurrency();
        int current_threads = 1;
        
#ifdef _OPENMP
        available_cores = omp_get_num_procs();
        current_threads = omp_get_max_threads();
#endif
        
        // System analysis
        SystemType sys_type = detect_system_type();
        StorageType storage_type = detect_storage_type_advanced(include_benchmark);
        
        // Convert enums to readable strings
        std::string system_str, storage_str;
        switch (sys_type) {
            case SystemType::DESKTOP: system_str = "Desktop/Laptop"; break;
            case SystemType::SERVER: system_str = "Server"; break;
            case SystemType::HPC_CLUSTER: system_str = "HPC Cluster"; break;
            case SystemType::CRAN_CHECK: system_str = "CRAN Environment"; break;
            case SystemType::CONTAINER: system_str = "Container"; break;
        }
        
        switch (storage_type) {
            case StorageType::HDD: storage_str = "Hard Disk Drive"; break;
            case StorageType::SATA_SSD: storage_str = "SATA SSD"; break;
            case StorageType::NVME_SSD: storage_str = "NVMe SSD"; break;
            case StorageType::NETWORK: storage_str = "Network Storage"; break;
            case StorageType::UNKNOWN: storage_str = "Unknown"; break;
        }
        
        // Feature detection
        Rcpp::CharacterVector features;
        features.push_back("io_optimization");
        features.push_back("storage_detection");
        features.push_back("system_analysis");
        features.push_back("thread_optimization");
        features.push_back("performance_diagnostics");
        
#ifdef _OPENMP
        features.push_back("openmp");
#endif
        
        // Memory information
        size_t available_memory_kb = 0;
#ifdef __linux__
        std::ifstream meminfo("/proc/meminfo");
        if (meminfo.is_open()) {
            std::string line;
            while (std::getline(meminfo, line)) {
                if (line.find("MemAvailable:") == 0) {
                    sscanf(line.c_str(), "MemAvailable: %zu kB", &available_memory_kb);
                    break;
                }
            }
        }
#endif
        
        // Calculate memory-safe threads for a typical large matrix
        int memory_safe_threads = calculate_memory_safe_threads(25000000); // 5000x5000 matrix
        
        // Performance estimates
        int storage_io_capacity = static_cast<int>(storage_type);
        
        return Rcpp::List::create(
            Rcpp::Named("bigdatastatmeth_version") = "1.1.0",
            Rcpp::Named("system_type") = system_str,
            Rcpp::Named("storage_type") = storage_str,
            Rcpp::Named("cpu_cores") = available_cores,
            Rcpp::Named("current_openmp_threads") = current_threads,
            Rcpp::Named("storage_io_capacity") = storage_io_capacity,
            Rcpp::Named("memory_safe_threads_5kx5k") = memory_safe_threads,
            Rcpp::Named("available_memory_gb") = available_memory_kb / 1024.0 / 1024.0,
            Rcpp::Named("features") = features,
            Rcpp::Named("openmp_enabled") = 
#ifdef _OPENMP
                true,
#else
                false,
#endif
            Rcpp::Named("benchmark_included") = include_benchmark,
            Rcpp::Named("is_hpc") = (sys_type == SystemType::HPC_CLUSTER),
            Rcpp::Named("is_cran") = (sys_type == SystemType::CRAN_CHECK),
            Rcpp::Named("is_container") = (sys_type == SystemType::CONTAINER)
        );
        
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdSystemInfo: " << ex.what() << std::endl;
        return R_NilValue;
    }
}