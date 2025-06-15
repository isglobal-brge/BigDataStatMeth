/**
 * @file performance_integration.hpp
 * @brief Integration helpers and performance wrappers for BigDataStatMeth
 * 
 * Header-only integration functions that provide seamless integration
 * of I/O-aware optimization into existing BigDataStatMeth operations
 * and external C++ code.
 */

#ifndef BIGDATASTATMETH_PERFORMANCE_INTEGRATION_HPP
#define BIGDATASTATMETH_PERFORMANCE_INTEGRATION_HPP

namespace BigDataStatMeth {

/**
 * @brief Auto-configure based on operation parameters
 * 
 * Smart configuration function that automatically detects optimal
 * thread count and applies configuration. Designed for internal
 * use by BigDataStatMeth operations.
 * 
 * @param matrix_rows Number of matrix rows
 * @param matrix_cols Number of matrix columns
 * @param operation_type Operation type being performed
 * @param user_threads User-specified threads (if any)
 * @return Configured thread count for immediate use
 */
inline int auto_configure_threads(int matrix_rows, int matrix_cols,
                                 const std::string& operation_type,
                                 int user_threads) {
    
    // Fast path: if user specified threads, use them directly
    if (user_threads > 0) {
        int available_cores = std::thread::hardware_concurrency();
#ifdef _OPENMP
        available_cores = omp_get_num_procs();
        omp_set_num_threads(std::min(user_threads, available_cores));
#endif
        return std::min(user_threads, available_cores);
    }
    
    // Get optimal threads using enhanced analysis
    int optimal_threads = get_io_aware_threads_enhanced(R_NilValue, matrix_rows, matrix_cols, 
                                                       operation_type, false);
    
#ifdef _OPENMP
    // Apply the configuration
    omp_set_num_threads(optimal_threads);
    omp_set_dynamic(false);
    
    // Quick scheduling decision
    if (optimal_threads <= 4) {
        omp_set_schedule(omp_sched_static, 0);
    } else {
        omp_set_schedule(omp_sched_dynamic, 1000);
    }
#endif
    
    return optimal_threads;
}

// /**
//  * @brief Integration helper for BigDataStatMeth matrix operations
//  * 
//  * Provides standardized integration of I/O-aware optimization
//  * into existing BigDataStatMeth functions with minimal code changes.
//  * 
//  * Usage pattern:
//  * @code
//  * int num_threads = integrate_io_optimization(nrows, ncols, "multiplication");
//  * #pragma omp parallel num_threads(num_threads)
//  * {
//  *     // Your matrix operation code here
//  * }
//  * @endcode
//  * 
//  * @param matrix_rows Number of matrix rows
//  * @param matrix_cols Number of matrix columns
//  * @param operation_type Type of operation being performed
//  * @param user_threads User-specified thread count (optional)
//  * @param verbose Whether to print optimization information
//  * @return Optimal thread count for the operation
//  */
// inline int integrate_io_optimization(int matrix_rows, int matrix_cols,
//                                     const std::string& operation_type,
//                                     int user_threads, bool verbose) {
//     
//     if (verbose) {
//         diagnose_io_cpp(matrix_rows, matrix_cols, operation_type, false, true);
//     }
//     
//     return auto_configure_threads(matrix_rows, matrix_cols, operation_type, user_threads);
// }

/**
 * @brief Performance monitoring wrapper for matrix operations
 * 
 * Wraps any matrix operation with performance monitoring and
 * automatic I/O-aware optimization. Template-based for maximum
 * flexibility with different function signatures.
 * 
 * @tparam Func Function type to wrap
 * @tparam Args Argument types
 * @param operation_name Name of the operation for logging
 * @param matrix_rows Number of matrix rows
 * @param matrix_cols Number of matrix columns
 * @param func Function to execute with optimization
 * @param args Arguments to pass to the function
 * @return Result of the wrapped function
 */
template<typename Func, typename... Args>
inline auto performance_wrapper(const std::string& operation_name,
                               int matrix_rows, int matrix_cols,
                               Func&& func, Args&&... args) -> decltype(func(args...)) {
    
    // Apply I/O-aware optimization
    auto start_optimization = std::chrono::high_resolution_clock::now();
    int threads = integrate_io_optimization(matrix_rows, matrix_cols, operation_name);
    auto end_optimization = std::chrono::high_resolution_clock::now();
    
    // Execute the operation
    auto start_execution = std::chrono::high_resolution_clock::now();
    auto result = func(std::forward<Args>(args)...);
    auto end_execution = std::chrono::high_resolution_clock::now();
    
    // Calculate timings
    auto optimization_time = std::chrono::duration_cast<std::chrono::microseconds>(
        end_optimization - start_optimization).count();
    auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_execution - start_execution).count();
    
    // Report performance (optional, can be controlled by environment variable)
    const char* verbose_env = std::getenv("BIGDATASTATMETH_VERBOSE");
    if (verbose_env && std::string(verbose_env) == "1") {
        std::cout << "\nPerformance Report for " << operation_name << ":\n";
        std::cout << "  Matrix size: " << matrix_rows << "x" << matrix_cols << "\n";
        std::cout << "  Threads used: " << threads << "\n";
        std::cout << "  Optimization time: " << optimization_time << " μs\n";
        std::cout << "  Execution time: " << execution_time << " ms\n";
        
        if (execution_time > 1000) {
            double seconds = execution_time / 1000.0;
            std::cout << "  Total time: " << std::fixed << std::setprecision(2) 
                      << seconds << " seconds\n";
        }
        std::cout << "\n";
    }
    
    return result;
}

/**
 * @brief Batch optimization for multiple operations
 * 
 * Analyzes and optimizes thread configuration for a sequence of operations.
 * Useful for workflows with multiple different matrix operations.
 * 
 * @param operations Vector of operation specifications
 * @param global_optimize Use global optimization across all operations
 * @return Pair of (global_threads, individual_threads_vector)
 */
inline std::pair<int, std::vector<int>> batch_optimize_cpp(
    const std::vector<std::tuple<std::string, int, int, double>>& operations,
    bool global_optimize = true) {
    
    std::vector<int> individual_threads;
    std::vector<double> io_intensities;
    
    double total_weight = 0.0;
    double weighted_io_intensity = 0.0;
    
    // Analyze each operation
    for (const auto& op : operations) {
        std::string name = std::get<0>(op);
        int rows = std::get<1>(op);
        int cols = std::get<2>(op);
        double weight = std::get<3>(op);
        
        // Calculate individual optimization
        int threads = get_io_aware_threads(R_NilValue, rows, cols, name);
        size_t elements = static_cast<size_t>(rows) * cols;
        double io_intensity = calculate_io_intensity_ratio(elements, name);
        
        individual_threads.push_back(threads);
        io_intensities.push_back(io_intensity);
        
        total_weight += weight;
        weighted_io_intensity += io_intensity * weight;
    }
    
    // Calculate global optimization if requested
    int global_threads = 1;
    if (global_optimize && !operations.empty()) {
        // Use weighted average I/O intensity to determine global strategy
        double avg_io_intensity = weighted_io_intensity / total_weight;
        
        // Find typical matrix size (weighted average)
        double avg_rows = 0.0, avg_cols = 0.0;
        for (const auto& op : operations) {
            avg_rows += std::get<1>(op) * std::get<3>(op);
            avg_cols += std::get<2>(op) * std::get<3>(op);
        }
        avg_rows /= total_weight;
        avg_cols /= total_weight;
        
        global_threads = get_io_aware_threads(R_NilValue, 
                                             static_cast<int>(avg_rows), 
                                             static_cast<int>(avg_cols), 
                                             "mixed");
    }
    
    return std::make_pair(global_threads, individual_threads);
}

/**
 * @brief Enhanced matrix multiplication with I/O optimization
 * 
 * Template-based matrix multiplication that works with different matrix types
 * and automatically applies I/O-aware optimization.
 * 
 * @tparam MatrixType Type of input matrices
 * @param A First matrix
 * @param B Second matrix
 * @param user_threads Optional user-specified thread count
 * @param verbose Print optimization information
 * @return Product matrix A * B
 */
template<typename MatrixType>
inline MatrixType multiply_optimized(const MatrixType& A, const MatrixType& B,
                                    int user_threads = 0, bool verbose = false) {
    
    // Validate dimensions (assuming .rows() and .cols() methods)
    if (A.cols() != B.rows()) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }
    
    int rows_A = A.rows();
    int cols_A = A.cols();
    int cols_B = B.cols();
    
    // Apply I/O-aware optimization
    int optimal_threads = integrate_io_optimization(rows_A, cols_B, "multiplication", 
                                                   user_threads, verbose);
    
    // Create result matrix
    MatrixType C(rows_A, cols_B);
    
    // Perform optimized matrix multiplication
    auto start_time = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(optimal_threads) schedule(dynamic)
#endif
    for (int i = 0; i < rows_A; ++i) {
        for (int j = 0; j < cols_B; ++j) {
            typename MatrixType::Scalar sum = 0;
            for (int k = 0; k < cols_A; ++k) {
                sum += A(i, k) * B(k, j);
            }
            C(i, j) = sum;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    if (verbose) {
        std::cout << "Matrix multiplication completed in " << duration.count() 
                  << " ms using " << optimal_threads << " threads\n";
    }
    
    return C;
}

/**
 * @brief Enhanced cross-product with I/O optimization
 * 
 * Template-based cross-product computation (t(A) * A) with automatic optimization.
 * 
 * @tparam MatrixType Type of input matrix
 * @param A Input matrix
 * @param user_threads Optional user-specified thread count
 * @param verbose Print optimization information
 * @return Cross-product matrix t(A) * A
 */
template<typename MatrixType>
inline MatrixType crossprod_optimized(const MatrixType& A, int user_threads = 0, bool verbose = false) {
    
    int rows = A.rows();
    int cols = A.cols();
    
    // Cross-product has better cache behavior
    int optimal_threads = integrate_io_optimization(rows, cols, "crossprod", 
                                                   user_threads, verbose);
    
    MatrixType AtA(cols, cols);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(optimal_threads) schedule(static)
#endif
    for (int i = 0; i < cols; ++i) {
        for (int j = i; j < cols; ++j) {  // Only compute upper triangle
            typename MatrixType::Scalar sum = 0;
            for (int k = 0; k < rows; ++k) {
                sum += A(k, i) * A(k, j);
            }
            AtA(i, j) = sum;
            if (i != j) {
                AtA(j, i) = sum;  // Symmetric matrix
            }
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    if (verbose) {
        std::cout << "Cross-product completed in " << duration.count() 
                  << " ms using " << optimal_threads << " threads\n";
    }
    
    return AtA;
}

/**
 * @brief HDF5 operation optimization helper
 * 
 * Specialized optimization for HDF5-based operations that considers
 * file I/O characteristics and storage type.
 * 
 * @param file_path Path to HDF5 file
 * @param matrix_rows Estimated matrix rows
 * @param matrix_cols Estimated matrix columns
 * @param operation_type Type of HDF5 operation
 * @param user_threads User-specified threads (optional)
 * @return Optimal thread count for HDF5 operation
 */
inline int optimize_hdf5_operation(const std::string& file_path,
                                  int matrix_rows, int matrix_cols,
                                  const std::string& operation_type,
                                  int user_threads) {
    
    // Enhanced optimization for HDF5 operations
    StorageType storage_type = detect_storage_type_advanced(false);
    
    // Network storage is common for shared HDF5 files
    if (file_path.find("nfs") != std::string::npos || 
        file_path.find("gpfs") != std::string::npos ||
        file_path.find("lustre") != std::string::npos ||
        file_path.find("beegfs") != std::string::npos) {
        storage_type = StorageType::NETWORK;
    }
    
    // Apply optimization considering HDF5 I/O characteristics
    int optimal_threads = get_io_aware_threads_enhanced(
        user_threads > 0 ? Rcpp::wrap(user_threads) : R_NilValue, 
        matrix_rows, matrix_cols, operation_type, false);
    
    // Adjust for HDF5-specific considerations
    if (storage_type == StorageType::NETWORK) {
        // Network storage with HDF5 benefits from fewer threads
        optimal_threads = std::min(optimal_threads, 4);
    } else if (storage_type == StorageType::HDD) {
        // HDD with HDF5 should be very conservative
        optimal_threads = std::min(optimal_threads, 2);
    }
    
#ifdef _OPENMP
    omp_set_num_threads(optimal_threads);
    omp_set_dynamic(false);
    
    // HDF5 operations benefit from static scheduling
    omp_set_schedule(omp_sched_static, 0);
#endif
    
    return optimal_threads;
}

/**
 * @brief Memory monitoring for large operations
 * 
 * Monitors memory usage during large matrix operations and adjusts
 * thread count if memory pressure is detected.
 * 
 * @param matrix_elements Total number of matrix elements
 * @param current_threads Current thread count
 * @return Adjusted thread count based on memory pressure
 */
inline int monitor_memory_pressure(size_t matrix_elements, int current_threads) {
    
    // Get current memory usage
    size_t available_memory = 0;
    
#ifdef __linux__
    std::ifstream meminfo("/proc/meminfo");
    if (meminfo.is_open()) {
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.find("MemAvailable:") == 0) {
                size_t kb;
                sscanf(line.c_str(), "MemAvailable: %zu kB", &kb);
                available_memory = kb * 1024;  // Convert to bytes
                break;
            }
        }
    }
#endif
    
    if (available_memory == 0) {
        return current_threads;  // Can't monitor, keep current
    }
    
    // Estimate memory requirements
    size_t matrix_memory = matrix_elements * 8;  // Assume double precision
    size_t total_memory_needed = matrix_memory * 3;  // Input + intermediate + output
    size_t per_thread_memory = 64 * 1024 * 1024;  // 64MB per thread
    
    // Check if current configuration fits in available memory
    size_t total_estimated = total_memory_needed + (current_threads * per_thread_memory);
    
    if (total_estimated > available_memory * 0.8) {  // 80% memory threshold
        // Reduce threads to fit memory constraints
        size_t safe_memory = static_cast<size_t>(available_memory * 0.8);
        if (total_memory_needed < safe_memory) {
            size_t remaining = safe_memory - total_memory_needed;
            int safe_threads = static_cast<int>(remaining / per_thread_memory);
            return std::max(1, std::min(safe_threads, current_threads));
        } else {
            return 1;  // Even basic operation might be tight
        }
    }
    
    return current_threads;  // Memory looks good
}

/**
 * @brief Adaptive thread adjustment during operation
 * 
 * Dynamically adjusts thread count during long-running operations
 * based on system load and performance feedback.
 * 
 * @param operation_start_time Start time of the operation
 * @param initial_threads Initial thread count
 * @param target_efficiency Target efficiency threshold
 * @return Adjusted thread count
 */
inline int adaptive_thread_adjustment(std::chrono::high_resolution_clock::time_point operation_start_time,
                                     int initial_threads, double target_efficiency) {
    
    auto current_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - operation_start_time);
    
    // Only adjust for operations running longer than 10 seconds
    if (elapsed.count() < 10) {
        return initial_threads;
    }
    
    // Check system load (Linux only)
#ifdef __linux__
    std::ifstream loadavg("/proc/loadavg");
    if (loadavg.is_open()) {
        double load1, load5, load15;
        loadavg >> load1 >> load5 >> load15;
        
        int available_cores = std::thread::hardware_concurrency();
        double load_per_core = load1 / available_cores;
        
        // If system is heavily loaded, reduce threads
        if (load_per_core > 1.5) {
            return std::max(1, initial_threads / 2);
        }
        // If system is lightly loaded, could potentially increase
        else if (load_per_core < 0.5 && initial_threads < available_cores / 2) {
            return std::min(available_cores / 2, initial_threads + 1);
        }
    }
#endif
    
    return initial_threads;  // No adjustment needed
}

/**
 * @brief Integration macro for existing BigDataStatMeth functions
 * 
 * Provides a simple macro that can be added to existing functions
 * with minimal code changes.
 */
#define BIGDATASTATMETH_OPTIMIZE(rows, cols, operation) \
    int _bd_optimal_threads = BigDataStatMeth::integrate_io_optimization(rows, cols, operation); \
    _Pragma("omp parallel num_threads(_bd_optimal_threads)")

/**
 * @brief RAII class for automatic thread management
 * 
 * Automatically applies and restores thread configuration using RAII pattern.
 * Useful for temporary optimizations within specific scopes.
 */
class ScopedThreadOptimizer {
private:
    int original_threads_;
    bool openmp_available_;
    
public:
    /**
     * @brief Constructor that applies optimization
     * 
     * @param matrix_rows Number of matrix rows
     * @param matrix_cols Number of matrix columns
     * @param operation_type Type of operation
     */
    ScopedThreadOptimizer(int matrix_rows, int matrix_cols, const std::string& operation_type) 
        : original_threads_(1), openmp_available_(false) {
        
#ifdef _OPENMP
        openmp_available_ = true;
        original_threads_ = omp_get_max_threads();
        
        int optimal_threads = integrate_io_optimization(matrix_rows, matrix_cols, operation_type);
        omp_set_num_threads(optimal_threads);
#endif
    }
    
    /**
     * @brief Destructor that restores original configuration
     */
    ~ScopedThreadOptimizer() {
#ifdef _OPENMP
        if (openmp_available_) {
            omp_set_num_threads(original_threads_);
        }
#endif
    }
    
    // Non-copyable, non-movable
    ScopedThreadOptimizer(const ScopedThreadOptimizer&) = delete;
    ScopedThreadOptimizer& operator=(const ScopedThreadOptimizer&) = delete;
    ScopedThreadOptimizer(ScopedThreadOptimizer&&) = delete;
    ScopedThreadOptimizer& operator=(ScopedThreadOptimizer&&) = delete;
};

/**
 * @brief Performance logger for operations
 * 
 * Logs performance metrics to file for analysis and optimization tracking.
 */
class PerformanceLogger {
private:
    std::string log_file_;
    bool enabled_;
    
public:
    /**
     * @brief Constructor
     * 
     * @param log_file Path to log file
     * @param enabled Whether logging is enabled
     */
    PerformanceLogger(const std::string& log_file = "bigdatastatmeth_performance.log", bool enabled = true)
        : log_file_(log_file), enabled_(enabled) {
        
        if (enabled_ && !std::ifstream(log_file_).good()) {
            // Create log file with header
            std::ofstream ofs(log_file_);
            if (ofs.is_open()) {
                ofs << "timestamp,operation,matrix_rows,matrix_cols,threads_used,execution_time_ms,system_load\n";
            }
        }
    }
    
    /**
     * @brief Log performance data
     * 
     * @param operation Operation name
     * @param matrix_rows Number of matrix rows
     * @param matrix_cols Number of matrix columns
     * @param threads_used Number of threads used
     * @param execution_time_ms Execution time in milliseconds
     */
    void log(const std::string& operation, int matrix_rows, int matrix_cols, 
             int threads_used, double execution_time_ms) {
        
        if (!enabled_) return;
        
        std::ofstream ofs(log_file_, std::ios::app);
        if (!ofs.is_open()) return;
        
        // Get current timestamp
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        
        // Get system load (if available)
        double system_load = 0.0;
#ifdef __linux__
        std::ifstream loadavg("/proc/loadavg");
        if (loadavg.is_open()) {
            loadavg >> system_load;
        }
#endif
        
        // Write log entry
        ofs << time_t << "," << operation << "," << matrix_rows << "," << matrix_cols 
            << "," << threads_used << "," << execution_time_ms << "," << system_load << "\n";
    }
    
    /**
     * @brief Enable or disable logging
     */
    void set_enabled(bool enabled) { enabled_ = enabled; }
    
    /**
     * @brief Check if logging is enabled
     */
    bool is_enabled() const { return enabled_; }
};

// Global performance logger instance
static PerformanceLogger global_performance_logger;

/**
 * @brief Enable global performance logging
 * 
 * @param log_file Path to log file
 */
inline void enable_performance_logging(const std::string& log_file) {
    global_performance_logger = PerformanceLogger(log_file, true);
}

/**
 * @brief Disable global performance logging
 */
inline void disable_performance_logging() {
    global_performance_logger.set_enabled(false);
}

/**
 * @brief Complete integration example for matrix multiplication
 * 
 * Shows how to integrate I/O-aware optimization into a complete
 * matrix operation with all features enabled.
 * 
 * @tparam MatrixType Type of matrices
 * @param A First matrix
 * @param B Second matrix
 * @param operation_name Name for logging
 * @param user_threads Optional user-specified threads
 * @param enable_logging Enable performance logging
 * @param enable_monitoring Enable memory monitoring
 * @return Result matrix
 */
template<typename MatrixType>
inline MatrixType complete_matrix_multiply(const MatrixType& A, const MatrixType& B,
                                          const std::string& operation_name = "matrix_multiply",
                                          int user_threads = 0,
                                          bool enable_logging = false,
                                          bool enable_monitoring = true) {
    
    // Validate dimensions
    if (A.cols() != B.rows()) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }
    
    int rows_A = A.rows();
    int cols_A = A.cols();
    int cols_B = B.cols();
    size_t total_elements = static_cast<size_t>(rows_A) * cols_B;
    
    // Apply I/O-aware optimization with RAII
    ScopedThreadOptimizer optimizer(rows_A, cols_B, "multiplication");
    
    // Get actual thread count being used
    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = omp_get_max_threads();
#endif
    
    // Memory monitoring
    if (enable_monitoring) {
        actual_threads = monitor_memory_pressure(total_elements, actual_threads);
#ifdef _OPENMP
        omp_set_num_threads(actual_threads);
#endif
    }
    
    // Create result matrix
    MatrixType C(rows_A, cols_B);
    
    // Perform multiplication with timing
    auto start_time = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(actual_threads) schedule(dynamic)
#endif
    for (int i = 0; i < rows_A; ++i) {
        for (int j = 0; j < cols_B; ++j) {
            typename MatrixType::Scalar sum = 0;
            for (int k = 0; k < cols_A; ++k) {
                sum += A(i, k) * B(k, j);
            }
            C(i, j) = sum;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Performance logging
    if (enable_logging) {
        global_performance_logger.log(operation_name, rows_A, cols_B, 
                                     actual_threads, static_cast<double>(duration.count()));
    }
    
    return C;
}

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_PERFORMANCE_INTEGRATION_HPP