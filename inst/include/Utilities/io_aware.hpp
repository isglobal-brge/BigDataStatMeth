/**
 * @file io-aware-openmp.hpp
 * @brief I/O-aware OpenMP optimization for BigDataStatMeth
 * 
 * Public API header for external libraries using BigDataStatMeth.
 * Located in inst/include/ for Rcpp LinkingTo compatibility.
 * 
 * This header provides I/O-aware thread optimization that prevents
 * libgomp memory overflow on large systems while maximizing performance
 * based on storage characteristics and workload analysis.
 * 
 * @author BigDataStatMeth Team
 * @version 1.0
 */

#ifndef BIGDATASTATMETH_IO_AWARE_OPENMP_HPP
#define BIGDATASTATMETH_IO_AWARE_OPENMP_HPP

// #include <Rcpp.h>
// #include <algorithm>
// #include <string>
// #include <fstream>
// 
// #ifdef _OPENMP
// #include <omp.h>
// #endif

namespace BigDataStatMeth {

/**
 * @brief Detect storage type and optimal I/O thread count
 * 
 * Analyzes the filesystem and storage characteristics to determine
 * the optimal number of threads for I/O operations. Uses heuristics
 * based on mount information and environment variables.
 * 
 * @return Optimal threads for I/O operations on detected storage
 * @retval 1 Conservative default (HDD or unknown storage)
 * @retval 4 SATA SSD detected
 * @retval 6 Network storage detected  
 * @retval 8 NVMe SSD detected
 * 
 * @note On non-Linux systems, returns conservative default of 4
 * @since BigDataStatMeth 1.0
 */
inline int detect_storage_io_threads() {
    
#ifdef __linux__
    // Detect storage type from filesystem information
    std::ifstream mounts("/proc/mounts");
    std::string line;
    
    // Check for high-performance storage first
    while (std::getline(mounts, line)) {
        if (line.find("nvme") != std::string::npos) {
            return 8;  // NVMe: Can handle aggressive parallel I/O
        } else if (line.find("ssd") != std::string::npos || 
                   line.find("ext4") != std::string::npos) {
            return 4;  // SSD: Moderate parallel I/O capability
        }
    }
    
    // Check for network/distributed storage
    const char* pwd = std::getenv("PWD");
    if (pwd) {
        std::string current_path(pwd);
        if (current_path.find("nfs") != std::string::npos ||
            current_path.find("gpfs") != std::string::npos ||
            current_path.find("lustre") != std::string::npos ||
            current_path.find("beegfs") != std::string::npos) {
            return 6;  // Network storage: Variable but generally good
        }
    }
    
    return 1;  // Conservative default (assume HDD)
#else
    return 4;  // Conservative default for non-Linux systems
#endif
}

/**
 * @brief Calculate I/O to computation ratio for BigDataStatMeth operations
 * 
 * Estimates the ratio of time spent on I/O operations versus computation
 * for different BigDataStatMeth operation types. Higher ratios indicate
 * I/O-bound operations that benefit from fewer threads to reduce contention.
 * 
 * @param matrix_elements Total number of matrix elements (rows × cols)
 * @param operation_type Type of BigDataStatMeth operation
 * @return Ratio indicating I/O intensity
 * @retval >4.0 Heavily I/O bound - minimize threads
 * @retval 2.0-4.0 Moderately I/O bound - balance threads
 * @retval 0.8-2.0 Compute bound - more threads acceptable  
 * @retval <0.8 Heavily compute bound - aggressive parallelization
 * 
 * @details
 * Supported operation types:
 * - "multiplication", "bdblockmult": Matrix A × B operations
 * - "svd", "bdSVD": Singular Value Decomposition
 * - "crossprod", "bdCrossprod": t(A) × A cross-product
 * - "tcrossprod", "bdtCrossprod": A × t(A) transpose cross-product
 * - "pca", "bdPCA": Principal Component Analysis
 * 
 * @note The ratio is adjusted based on matrix size, with larger matrices
 *       being more I/O intensive due to memory hierarchy effects
 * @since BigDataStatMeth 1.0
 */
inline double calculate_io_intensity_ratio(size_t matrix_elements, const std::string& operation_type) {
    
    // Base I/O factor for different BigDataStatMeth operations
    double base_io_factor;
    
    if (operation_type == "multiplication" || operation_type == "bdblockmult") {
        base_io_factor = 3.0;  // Read matrix A + Read matrix B + Write result C
    } else if (operation_type == "svd" || operation_type == "bdSVD") {
        base_io_factor = 4.5;  // Read matrix + Write U,S,V + Intermediate operations
    } else if (operation_type == "crossprod" || operation_type == "bdCrossprod") {
        base_io_factor = 1.8;  // Read matrix + Write t(A) %*% A result
    } else if (operation_type == "tcrossprod" || operation_type == "bdtCrossprod") {
        base_io_factor = 1.8;  // Read matrix + Write A %*% t(A) result
    } else if (operation_type == "pca" || operation_type == "bdPCA") {
        base_io_factor = 3.5;  // Read data + Write loadings/scores + Intermediate
    } else {
        base_io_factor = 2.5;  // Generic BigDataStatMeth operation
    }
    
    // Adjust based on matrix size (larger matrices are more I/O intensive)
    double size_multiplier;
    if (matrix_elements < 1000000) {        // < 1M elements (small matrices)
        size_multiplier = 0.3;              // Compute dominates
    } else if (matrix_elements < 25000000) { // < 25M elements (medium matrices)
        size_multiplier = 1.0;              // Balanced
    } else if (matrix_elements < 100000000) { // < 100M elements (large matrices)
        size_multiplier = 2.0;              // I/O becomes significant
    } else {
        size_multiplier = 3.5;              // Very large matrices: heavily I/O bound
    }
    
    return base_io_factor * size_multiplier;
}

/**
 * @brief Main I/O-aware thread calculation function
 * 
 * This function replaces get_optimal_threads() in BigDataStatMeth operations
 * with intelligent analysis of I/O bottlenecks, storage characteristics,
 * and system constraints to prevent performance issues on large systems.
 * 
 * The optimization strategy considers:
 * - Storage type (HDD/SSD/NVMe/Network) and I/O capacity
 * - Matrix size and operation type to estimate I/O intensity
 * - System characteristics (HPC cluster, large server, desktop)
 * - Environment constraints (CRAN, job schedulers, memory limits)
 * 
 * @param user_threads User-specified thread count (optional)
 * @param matrix_rows Number of matrix rows
 * @param matrix_cols Number of matrix columns  
 * @param operation_type BigDataStatMeth operation type (default: "multiplication")
 * @return Optimal number of threads considering all factors
 * 
 * @details
 * Optimization strategies:
 * - **I/O BOUND** (ratio >4.0): 1-2 threads to minimize contention
 * - **BALANCED** (ratio 2.0-4.0): Storage capacity × 2, limited by cores/4
 * - **COMPUTE BOUND** (ratio 0.8-2.0): cores/2, max 32 threads
 * - **CPU INTENSIVE** (ratio <0.8): aggressive parallelization, max 64 threads
 * 
 * @note
 * - Automatically detects and respects HPC job scheduler allocations
 * - Prevents libgomp memory overflow on systems with >64 cores
 * - Maintains CRAN compliance when OMP_THREAD_LIMIT ≤ 2
 * - Returns minimum 1 thread, maximum available_cores
 * 
 * @see detect_storage_io_threads() for storage detection details
 * @see calculate_io_intensity_ratio() for I/O analysis details
 * @since BigDataStatMeth 1.0
 */
inline int get_io_aware_threads(Rcpp::Nullable<int> user_threads,
                               int matrix_rows, int matrix_cols,
                               const std::string& operation_type = "multiplication") {
    
    // Get system information
    int available_cores = 1;
#ifdef _OPENMP
    available_cores = omp_get_num_procs();
#endif
    
    // Honor user specification if provided
    if (user_threads.isNotNull()) {
        int requested = Rcpp::as<int>(user_threads);
        return std::min(requested, available_cores);
    }
    
    // Calculate workload characteristics
    size_t matrix_elements = (size_t)matrix_rows * matrix_cols;
    double io_intensity = calculate_io_intensity_ratio(matrix_elements, operation_type);
    int storage_io_capacity = detect_storage_io_threads();
    
    // Core optimization logic: Choose strategy based on I/O intensity
    int optimal_threads;
    
    if (io_intensity > 4.0) {
        // HEAVILY I/O BOUND: Minimize thread contention on I/O critical sections
        optimal_threads = std::min(2, storage_io_capacity);
        
    } else if (io_intensity > 2.0) {
        // MODERATELY I/O BOUND: Balance I/O capacity with computation needs
        optimal_threads = std::min(storage_io_capacity * 2, available_cores / 4);
        
    } else if (io_intensity > 0.8) {
        // COMPUTE BOUND: More threads acceptable, but respect I/O constraints
        optimal_threads = std::min(available_cores / 2, 32);
        
    } else {
        // HEAVILY COMPUTE BOUND: Aggressive parallelization with I/O safety limit
        optimal_threads = std::min(available_cores, 64);  // Hard limit for libgomp safety
    }
    
    // Apply HPC environment constraints
    bool is_hpc = (std::getenv("SLURM_JOB_ID") != nullptr) ||
                  (std::getenv("PBS_JOBID") != nullptr) ||
                  (std::getenv("LSB_JOBID") != nullptr) ||
                  (std::getenv("SGE_JOB_ID") != nullptr) ||
                  (std::getenv("SLURM_CPUS_PER_TASK") != nullptr) ||
                  (std::getenv("PBS_NCPUS") != nullptr);
    
    if (is_hpc) {
        // Respect job scheduler allocations
        const char* slurm_cpus = std::getenv("SLURM_CPUS_PER_TASK");
        const char* pbs_cpus = std::getenv("PBS_NCPUS");
        const char* omp_threads = std::getenv("OMP_NUM_THREADS");
        
        if (slurm_cpus) {
            optimal_threads = std::min(optimal_threads, std::atoi(slurm_cpus));
        } else if (pbs_cpus) {
            optimal_threads = std::min(optimal_threads, std::atoi(pbs_cpus));
        } else if (omp_threads) {
            optimal_threads = std::min(optimal_threads, std::atoi(omp_threads));
        }
    }
    
    // Large system safety: Prevent libgomp memory overflow
    if (available_cores > 64) {
        optimal_threads = std::min(optimal_threads, 32);
    } else if (available_cores > 32) {
        optimal_threads = std::min(optimal_threads, 24);
    }
    
    // CRAN compliance check
    const char* thread_limit = std::getenv("OMP_THREAD_LIMIT");
    if (thread_limit && std::atoi(thread_limit) <= 2) {
        optimal_threads = 2;  // CRAN environment detected
    }
    
    // Final safety bounds
    return std::max(1, std::min(optimal_threads, available_cores));
}

/**
 * @brief Diagnostic function for I/O-aware thread optimization
 * 
 * Provides detailed analysis of the thread optimization decision process,
 * including system detection, storage characteristics, I/O intensity analysis,
 * and the chosen optimization strategy. Useful for debugging performance 
 * issues and understanding system behavior.
 * 
 * @param matrix_rows Number of matrix rows for analysis
 * @param matrix_cols Number of matrix columns for analysis
 * @param operation_type BigDataStatMeth operation type (default: "multiplication")
 * 
 * @details
 * Output includes:
 * - System type detection (HPC/Server/Desktop/CRAN)
 * - Storage analysis and I/O thread capacity
 * - I/O intensity ratio and optimization strategy
 * - Thread recommendations and reasoning
 * - Performance expectations for large systems
 * 
 * @note This function only prints diagnostics and does not modify
 *       any OpenMP or system settings
 * @since BigDataStatMeth 1.0
 */
inline void print_io_thread_diagnosis(int matrix_rows, int matrix_cols, 
                                     const std::string& operation_type = "multiplication") {
    
    size_t elements = (size_t)matrix_rows * matrix_cols;
    double io_intensity = calculate_io_intensity_ratio(elements, operation_type);
    int storage_capacity = detect_storage_io_threads();
    int recommended = get_io_aware_threads(R_NilValue, matrix_rows, matrix_cols, operation_type);
    
    int available_cores = 1;
#ifdef _OPENMP
    available_cores = omp_get_num_procs();
#endif
    
    // System environment detection
    bool is_hpc = (std::getenv("SLURM_JOB_ID") != nullptr) ||
                  (std::getenv("PBS_JOBID") != nullptr);
    bool is_cran = (std::getenv("OMP_THREAD_LIMIT") != nullptr) &&
                   (std::atoi(std::getenv("OMP_THREAD_LIMIT")) <= 2);
    
    Rcpp::Rcout << "\n=== BigDataStatMeth I/O-Aware Thread Analysis ===\n";
    Rcpp::Rcout << "Matrix: " << matrix_rows << "x" << matrix_cols 
               << " (" << elements << " elements)\n";
    Rcpp::Rcout << "Operation: " << operation_type << "\n";
    Rcpp::Rcout << "System: " << available_cores << " cores available";
    
    if (is_hpc) Rcpp::Rcout << " (HPC environment)";
    if (is_cran) Rcpp::Rcout << " (CRAN environment)";
    Rcpp::Rcout << "\n";
    
    Rcpp::Rcout << "\nAnalysis:\n";
    Rcpp::Rcout << "  I/O intensity ratio: " << std::fixed << std::setprecision(2) 
               << io_intensity << "\n";
    Rcpp::Rcout << "  Storage I/O capacity: " << storage_capacity << " threads\n";
    Rcpp::Rcout << "  Recommended threads: " << recommended << "\n";
    
    // Explain the strategy
    Rcpp::Rcout << "\nOptimization strategy: ";
    if (io_intensity > 4.0) {
        Rcpp::Rcout << "I/O BOUND - Minimize thread contention\n";
        Rcpp::Rcout << "  → Using minimal threads to avoid I/O critical section bottlenecks\n";
    } else if (io_intensity > 2.0) {
        Rcpp::Rcout << "BALANCED - I/O and compute optimization\n";
        Rcpp::Rcout << "  → Balancing I/O capacity with computational parallelism\n";
    } else if (io_intensity > 0.8) {
        Rcpp::Rcout << "COMPUTE BOUND - CPU-focused parallelization\n";
        Rcpp::Rcout << "  → Using more threads while respecting I/O constraints\n";
    } else {
        Rcpp::Rcout << "CPU INTENSIVE - Aggressive parallelization\n";
        Rcpp::Rcout << "  → Maximum parallelization with libgomp safety limits\n";
    }
    
    // Performance expectations
    if (available_cores > 64 && recommended < available_cores/3) {
        Rcpp::Rcout << "\nNote: Large system detected (" << available_cores 
                   << " cores) - threads limited to prevent:\n";
        Rcpp::Rcout << "  • libgomp memory overflow\n";
        Rcpp::Rcout << "  • I/O critical section contention\n";
        Rcpp::Rcout << "  • Memory bandwidth saturation\n";
        
        double improvement = (double)available_cores / recommended * 1.5;
        Rcpp::Rcout << "  Expected performance improvement: " 
                   << std::fixed << std::setprecision(1) << improvement << "x faster\n";
    }
    
    Rcpp::Rcout << "================================================\n\n";
}

/**
 * @brief Configure BigDataStatMeth with I/O-aware optimization
 * 
 * Applies I/O-aware thread configuration globally for the current session
 * using OpenMP settings. This function calculates optimal thread count
 * and applies appropriate OpenMP scheduling and configuration.
 * 
 * @param matrix_rows Typical matrix rows for your workload
 * @param matrix_cols Typical matrix columns for your workload
 * @param operation_type Primary operation type (default: "multiplication")
 * 
 * @details
 * Applied OpenMP settings:
 * - Thread count: Based on I/O-aware analysis
 * - Dynamic threads: Disabled for consistent performance
 * - Nested parallelism: Disabled to avoid overhead
 * - Scheduling: Static for HPC/few threads, dynamic for desktop
 * 
 * @note
 * - Configuration persists for the current R session
 * - Call this function once at the beginning of your analysis
 * - Subsequent BigDataStatMeth operations will use these settings
 * 
 * @warning This function modifies global OpenMP settings which may
 *          affect other OpenMP-enabled code in the same session
 * @since BigDataStatMeth 1.0
 */
inline void configure_bigdata_io_threads(int matrix_rows, int matrix_cols,
                                        const std::string& operation_type = "multiplication") {
    
    int optimal_threads = get_io_aware_threads(R_NilValue, matrix_rows, matrix_cols, operation_type);
    
#ifdef _OPENMP
    omp_set_num_threads(optimal_threads);
    
    // Apply additional OpenMP optimizations
    omp_set_dynamic(false);  // Disable dynamic thread adjustment
    omp_set_nested(false);   // Disable nested parallelism
    
    // Set appropriate scheduling
    bool is_hpc = (std::getenv("SLURM_JOB_ID") != nullptr);
    if (is_hpc || optimal_threads <= 8) {
        omp_set_schedule(omp_sched_static, 0);  // Static for HPC or few threads
    } else {
        omp_set_schedule(omp_sched_dynamic, 1000);  // Dynamic for desktop
    }
#endif
    
    Rcpp::Rcout << "BigDataStatMeth I/O-aware configuration applied:\n";
    Rcpp::Rcout << "  Threads configured: " << optimal_threads << "\n";
    Rcpp::Rcout << "  Operation optimized for: " << operation_type << "\n";
    Rcpp::Rcout << "  Configuration active for this R session.\n\n";
}

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_IO_AWARE_OPENMP_HPP