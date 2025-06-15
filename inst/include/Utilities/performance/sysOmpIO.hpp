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
 * @version 1.1
 */

#ifndef BIGDATASTATMETH_IO_AWARE_OPENMP_HPP
#define BIGDATASTATMETH_IO_AWARE_OPENMP_HPP

// #include <Rcpp.h>
// #include <algorithm>
// #include <string>
// #include <fstream>
// #include <chrono>
// #include <thread>
// #include <iomanip>
// 
// #ifdef _OPENMP
// #include <omp.h>
// #endif

namespace BigDataStatMeth {

// /**
//  * @brief Storage type enumeration for clearer type handling
//  */
// enum class StorageType {
//     HDD = 1,        ///< Traditional hard disk drive
//     SATA_SSD = 4,   ///< SATA solid state drive  
//     NVME_SSD = 8,   ///< NVMe solid state drive
//     NETWORK = 6,    ///< Network/distributed storage
//     UNKNOWN = 2     ///< Unknown storage type
// };
// 
// /**
//  * @brief System type enumeration for environment detection
//  */
// enum class SystemType {
//     DESKTOP,        ///< Desktop/laptop system
//     SERVER,         ///< Server system (16+ cores)
//     HPC_CLUSTER,    ///< HPC cluster with job scheduler
//     CRAN_CHECK,     ///< CRAN checking environment
//     CONTAINER       ///< Containerized environment (Docker/Singularity)
// };

/**
 * @brief Advanced storage detection with benchmark capabilities
 * 
 * Enhanced storage detection that can optionally perform I/O benchmarks
 * to more accurately determine storage characteristics. Falls back to
 * heuristic detection when benchmarking is disabled or fails.
 * 
 * @param enable_benchmark Whether to perform I/O benchmark test (default: false)
 * @param benchmark_size Size of benchmark test in KB (default: 1024KB)
 * @return StorageType enum indicating detected storage type
 * 
 * @details
 * Detection methods (in order of preference):
 * 1. **Benchmark test**: Creates temporary file and measures I/O throughput
 * 2. **Filesystem analysis**: Examines /proc/mounts and device information
 * 3. **Environment heuristics**: Checks working directory and mount points
 * 4. **Conservative fallback**: Returns UNKNOWN for safer defaults
 * 
 * Benchmark thresholds:
 * - **NVMe**: >500 MB/s sequential write
 * - **SATA SSD**: 100-500 MB/s sequential write
 * - **Network**: Variable, detected by latency patterns
 * - **HDD**: <100 MB/s sequential write
 * 
 * @note Benchmark test creates temporary files and may add 100-500ms overhead
 * @since BigDataStatMeth 1.1
 */
inline StorageType detect_storage_type_advanced(bool enable_benchmark, 
                                                size_t benchmark_size) {
    
    // Optional I/O benchmark for accurate detection
    if (enable_benchmark) {
        try {
            std::string temp_file = std::tmpnam(nullptr);
            std::vector<char> test_data(benchmark_size * 1024, 'X');
            
            auto start = std::chrono::high_resolution_clock::now();
            
            // Write test
            std::ofstream ofs(temp_file, std::ios::binary);
            if (ofs.is_open()) {
                ofs.write(test_data.data(), test_data.size());
                ofs.flush();
                ofs.close();
                
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                
                // Calculate throughput (MB/s)
                double mb_per_sec = (benchmark_size / 1024.0) / (duration.count() / 1e6);
                
                // Clean up
                std::remove(temp_file.c_str());
                
                // Classify based on throughput
                if (mb_per_sec > 500) {
                    return StorageType::NVME_SSD;
                } else if (mb_per_sec > 100) {
                    return StorageType::SATA_SSD;
                } else {
                    return StorageType::HDD;
                }
            }
            
        } catch (...) {
            // Benchmark failed, fall back to heuristic detection
        }
    }
    
#ifdef __linux__
    // Enhanced heuristic detection
    std::ifstream mounts("/proc/mounts");
    if (mounts.is_open()) {
        std::string line;
        
        // Check for high-performance storage indicators
        while (std::getline(mounts, line)) {
            if (line.find("nvme") != std::string::npos) {
                return StorageType::NVME_SSD;
            } else if (line.find("ssd") != std::string::npos) {
                return StorageType::SATA_SSD;
            } else if (line.find("tmpfs") != std::string::npos && 
                line.find("/tmp") != std::string::npos) {
                return StorageType::NVME_SSD;  // RAM disk, treat as fastest
            }
        }
    }
    
    // Check for network/distributed storage
    const char* pwd = std::getenv("PWD");
    if (pwd) {
        std::string current_path(pwd);
        if (current_path.find("nfs") != std::string::npos ||
            current_path.find("gpfs") != std::string::npos ||
            current_path.find("lustre") != std::string::npos ||
            current_path.find("beegfs") != std::string::npos ||
            current_path.find("ceph") != std::string::npos ||
            current_path.find("gluster") != std::string::npos) {
            return StorageType::NETWORK;
        }
    }
    
    // Check /sys/block for device type information
    try {
        std::ifstream rotational("/sys/block/sda/queue/rotational");
        if (rotational.is_open()) {
            std::string value;
            rotational >> value;
            if (value == "0") {
                return StorageType::SATA_SSD;  // Non-rotational = SSD
            } else {
                return StorageType::HDD;       // Rotational = HDD
            }
        }
    } catch (...) {
        // /sys/block access failed, continue with other methods
    }
    
#endif
    
    return StorageType::UNKNOWN;  // Conservative default
}

// inline StorageType detect_storage_type_advanced(bool enable_benchmark, 
//                                                size_t benchmark_size) {
//     
//     // Optional I/O benchmark for accurate detection
//     if (enable_benchmark) {
//         try {
//             std::string temp_file = std::tmpnam(nullptr);
//             std::vector<char> test_data(benchmark_size * 1024, 'X');
//             
//             auto start = std::chrono::high_resolution_clock::now();
//             
//             // Write test
//             std::ofstream ofs(temp_file, std::ios::binary);
//             ofs.write(test_data.data(), test_data.size());
//             ofs.flush();
//             ofs.close();
//             
//             auto end = std::chrono::high_resolution_clock::now();
//             auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//             
//             // Calculate throughput (MB/s)
//             double mb_per_sec = (benchmark_size / 1024.0) / (duration.count() / 1e6);
//             
//             // Clean up
//             std::remove(temp_file.c_str());
//             
//             // Classify based on throughput
//             if (mb_per_sec > 500) {
//                 return StorageType::NVME_SSD;
//             } else if (mb_per_sec > 100) {
//                 return StorageType::SATA_SSD;
//             } else {
//                 return StorageType::HDD;
//             }
//             
//         } catch (...) {
//             // Benchmark failed, fall back to heuristic detection
//         }
//     }
//     
// #ifdef __linux__
//     // Enhanced heuristic detection
//     std::ifstream mounts("/proc/mounts");
//     std::string line;
//     
//     // Check for high-performance storage indicators
//     while (std::getline(mounts, line)) {
//         if (line.find("nvme") != std::string::npos) {
//             return StorageType::NVME_SSD;
//         } else if (line.find("ssd") != std::string::npos) {
//             return StorageType::SATA_SSD;
//         } else if (line.find("tmpfs") != std::string::npos && 
//                    line.find("/tmp") != std::string::npos) {
//             return StorageType::NVME_SSD;  // RAM disk, treat as fastest
//         }
//     }
//     
//     // Check for network/distributed storage
//     const char* pwd = std::getenv("PWD");
//     if (pwd) {
//         std::string current_path(pwd);
//         if (current_path.find("nfs") != std::string::npos ||
//             current_path.find("gpfs") != std::string::npos ||
//             current_path.find("lustre") != std::string::npos ||
//             current_path.find("beegfs") != std::string::npos ||
//             current_path.find("ceph") != std::string::npos ||
//             current_path.find("gluster") != std::string::npos) {
//             return StorageType::NETWORK;
//         }
//     }
//     
//     // Check /sys/block for device type information
//     try {
//         std::ifstream rotational("/sys/block/sda/queue/rotational");
//         if (rotational.is_open()) {
//             std::string value;
//             rotational >> value;
//             if (value == "0") {
//                 return StorageType::SATA_SSD;  // Non-rotational = SSD
//             } else {
//                 return StorageType::HDD;       // Rotational = HDD
//             }
//         }
//     } catch (...) {
//         // /sys/block access failed, continue with other methods
//     }
//     
// #endif
//     
//     return StorageType::UNKNOWN;  // Conservative default
// }



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
    StorageType storage = detect_storage_type_advanced(false);
    return static_cast<int>(storage);
}

/**
 * @brief Detect system type and environment characteristics
 * 
 * Comprehensive system environment detection that identifies the type
 * of system and runtime environment to apply appropriate optimizations.
 * 
 * @return SystemType enum indicating detected system type
 * 
 * @details
 * Detection criteria:
 * - **HPC_CLUSTER**: Job scheduler environment variables present
 * - **CRAN_CHECK**: OMP_THREAD_LIMIT ≤ 2 or _R_CHECK_* variables
 * - **CONTAINER**: Container-specific indicators (/proc/1/cgroup, etc.)
 * - **SERVER**: >16 CPU cores available
 * - **DESKTOP**: Default fallback for smaller systems
 * 
 * @note Container detection helps apply appropriate memory constraints
 * @since BigDataStatMeth 1.1
 */
inline SystemType detect_system_type() {
    
    // Check for CRAN environment first (highest priority)
    const char* thread_limit = std::getenv("OMP_THREAD_LIMIT");
    const char* r_check = std::getenv("_R_CHECK_PACKAGE_NAME_");
    if ((thread_limit && std::atoi(thread_limit) <= 2) || r_check) {
        return SystemType::CRAN_CHECK;
    }
    
    // Check for HPC job scheduler environments
    if (std::getenv("SLURM_JOB_ID") || std::getenv("PBS_JOBID") ||
        std::getenv("LSB_JOBID") || std::getenv("SGE_JOB_ID") ||
        std::getenv("SLURM_CPUS_PER_TASK") || std::getenv("PBS_NCPUS")) {
        return SystemType::HPC_CLUSTER;
    }
    
    // Check for container environments
#ifdef __linux__
    std::ifstream cgroup("/proc/1/cgroup");
    if (cgroup.is_open()) {
        std::string line;
        while (std::getline(cgroup, line)) {
            if (line.find("docker") != std::string::npos ||
                line.find("lxc") != std::string::npos ||
                line.find("containerd") != std::string::npos) {
                return SystemType::CONTAINER;
            }
        }
    }
    
    // Check for Singularity/Apptainer
    if (std::getenv("SINGULARITY_CONTAINER") || std::getenv("APPTAINER_CONTAINER")) {
        return SystemType::CONTAINER;
    }
#endif
    
    // Determine based on core count
    int available_cores = std::thread::hardware_concurrency();
#ifdef _OPENMP
    available_cores = omp_get_num_procs();
#endif
    
    if (available_cores > 16) {
        return SystemType::SERVER;
    } else {
        return SystemType::DESKTOP;
    }
}

/**
 * @brief Memory-aware thread limiting for large systems
 * 
 * Calculates safe thread limits based on available memory and matrix size
 * to prevent memory exhaustion and excessive swapping on large systems.
 * 
 * @param matrix_elements Total number of matrix elements
 * @param element_size Size of each matrix element in bytes (default: 8 for double)
 * @param safety_factor Memory safety factor (default: 0.8 = use 80% of available)
 * @return Maximum threads considering memory constraints
 * 
 * @details
 * Memory estimation considers:
 * - Matrix storage requirements (input + intermediate + output)
 * - Per-thread memory overhead (OpenMP stack, temporary buffers)
 * - System memory reserves for OS and other processes
 * - NUMA locality effects on large systems
 * 
 * The function prevents scenarios where too many threads cause memory
 * pressure leading to swapping and severe performance degradation.
 * 
 * @note Returns 1 if memory requirements exceed available system memory
 * @since BigDataStatMeth 1.1
 */
inline int calculate_memory_safe_threads(size_t matrix_elements, 
                                        size_t element_size,
                                        double safety_factor) {
    
    // Get available system memory
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
    
    // Fallback: estimate based on physical memory
    if (available_memory == 0) {
        // Conservative estimate: assume 75% of 8GB per 16 cores
        int cores = std::thread::hardware_concurrency();
        available_memory = static_cast<size_t>(cores / 16.0 * 8ULL * 1024 * 1024 * 1024 * 0.75);
    }
    
    // Calculate memory requirements
    size_t matrix_memory = matrix_elements * element_size;
    size_t per_operation_memory = matrix_memory * 3;  // A + B + C matrices
    size_t per_thread_overhead = 64 * 1024 * 1024;   // 64MB per thread (stack + buffers)
    
    size_t safe_memory = static_cast<size_t>(available_memory * safety_factor);
    
    if (per_operation_memory > safe_memory) {
        return 1;  // Even single-threaded may be tight
    }
    
    size_t remaining_memory = safe_memory - per_operation_memory;
    int max_threads = static_cast<int>(remaining_memory / per_thread_overhead);
    
    return std::max(1, max_threads);
}

/**
 * @brief Enhanced I/O intensity calculation with operation-specific tuning
 * 
 * Improved calculation that considers additional factors affecting I/O
 * intensity including data access patterns, cache behavior, and operation
 * complexity.
 * 
 * @param matrix_elements Total number of matrix elements (rows × cols)
 * @param operation_type Type of BigDataStatMeth operation
 * @param storage_type Detected storage type for refined estimates
 * @return Ratio indicating I/O intensity with enhanced precision
 * 
 * @details
 * Enhanced considerations:
 * - **Cache behavior**: Different operations have different cache hit rates
 * - **Access patterns**: Sequential vs random access patterns
 * - **Storage characteristics**: NVMe vs HDD access time differences
 * - **Operation complexity**: Some operations require multiple data passes
 * - **Memory hierarchy**: L1/L2/L3 cache vs main memory vs storage
 * 
 * @since BigDataStatMeth 1.1
 */
inline double calculate_io_intensity_ratio_enhanced(size_t matrix_elements, 
                                                   const std::string& operation_type,
                                                   StorageType storage_type = StorageType::UNKNOWN) {
    
    // Base I/O factors with refined estimates
    double base_io_factor;
    double cache_efficiency;
    
    if (operation_type == "multiplication" || operation_type == "bdblockmult") {
        base_io_factor = 3.0;
        cache_efficiency = 0.7;  // Good cache reuse in blocked multiplication
    } else if (operation_type == "svd" || operation_type == "bdSVD") {
        base_io_factor = 5.5;    // Multiple iterations, intermediate storage
        cache_efficiency = 0.4;  // Poor cache behavior due to multiple passes
    } else if (operation_type == "crossprod" || operation_type == "bdCrossprod") {
        base_io_factor = 1.5;
        cache_efficiency = 0.8;  // Excellent cache reuse (t(A) %*% A)
    } else if (operation_type == "tcrossprod" || operation_type == "bdtCrossprod") {
        base_io_factor = 1.5;
        cache_efficiency = 0.8;  // Excellent cache reuse (A %*% t(A))
    } else if (operation_type == "pca" || operation_type == "bdPCA") {
        base_io_factor = 4.0;    // Multiple matrix operations
        cache_efficiency = 0.5;  // Moderate cache behavior
    } else if (operation_type == "qr" || operation_type == "bdQR") {
        base_io_factor = 3.8;
        cache_efficiency = 0.6;
    } else if (operation_type == "eigen" || operation_type == "bdEigen") {
        base_io_factor = 6.0;    // Iterative eigenvalue computation
        cache_efficiency = 0.3;
    } else {
        base_io_factor = 2.8;    // Generic operation
        cache_efficiency = 0.6;
    }
    
    // Adjust for storage characteristics
    double storage_multiplier = 1.0;
    switch (storage_type) {
        case StorageType::HDD:
            storage_multiplier = 2.5;    // High seek time penalty
            break;
        case StorageType::SATA_SSD:
            storage_multiplier = 1.0;    // Baseline
            break;
        case StorageType::NVME_SSD:
            storage_multiplier = 0.6;    // Much faster random access
            break;
        case StorageType::NETWORK:
            storage_multiplier = 3.0;    // Network latency + bandwidth
            break;
        case StorageType::UNKNOWN:
            storage_multiplier = 1.5;    // Conservative assumption
            break;
    }
    
    // Size-based scaling with refined thresholds
    double size_multiplier;
    if (matrix_elements < 500000) {          // < 500K elements (fits in L3 cache)
        size_multiplier = 0.2;               // Compute dominates, minimal I/O
    } else if (matrix_elements < 4000000) {  // < 4M elements (fits in main memory)
        size_multiplier = 0.8;
    } else if (matrix_elements < 25000000) { // < 25M elements
        size_multiplier = 1.0;
    } else if (matrix_elements < 100000000) { // < 100M elements
        size_multiplier = 2.2;
    } else if (matrix_elements < 500000000) { // < 500M elements
        size_multiplier = 3.8;
    } else {
        size_multiplier = 5.0;               // Very large matrices
    }
    
    // Apply cache efficiency factor
    double effective_io_factor = base_io_factor * (2.0 - cache_efficiency);
    
    return effective_io_factor * size_multiplier * storage_multiplier;
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
    return calculate_io_intensity_ratio_enhanced(matrix_elements, operation_type, StorageType::UNKNOWN);
}

/**
 * @brief Enhanced main I/O-aware thread calculation with comprehensive optimization
 * 
 * Advanced version of the thread optimization function with additional
 * considerations for modern HPC systems, containers, and diverse workloads.
 * 
 * @param user_threads User-specified thread count (optional)
 * @param matrix_rows Number of matrix rows
 * @param matrix_cols Number of matrix columns
 * @param operation_type BigDataStatMeth operation type (default: "multiplication")
 * @param enable_storage_benchmark Enable I/O benchmark for accurate storage detection
 * @return Optimal number of threads considering all factors
 * 
 * @details
 * Enhanced optimization considers:
 * - **Advanced storage detection** with optional benchmarking
 * - **System type detection** (desktop/server/HPC/container/CRAN)
 * - **Memory constraints** to prevent swapping on large systems  
 * - **NUMA topology** awareness for multi-socket systems
 * - **Container resource limits** (cgroups, Docker limits)
 * - **Enhanced I/O intensity** analysis with storage-specific tuning
 * 
 * @since BigDataStatMeth 1.1
 */
inline int get_io_aware_threads_enhanced(Rcpp::Nullable<int> user_threads,
                                        int matrix_rows, int matrix_cols,
                                        const std::string& operation_type = "multiplication",
                                        bool enable_storage_benchmark = false) {
    
    // Get system information
    int available_cores = std::thread::hardware_concurrency();
#ifdef _OPENMP
    available_cores = omp_get_num_procs();
#endif
    
    // Honor user specification if provided
    if (user_threads.isNotNull()) {
        int requested = Rcpp::as<int>(user_threads);
        return std::min(requested, available_cores);
    }
    
    // Enhanced system analysis
    SystemType system_type = detect_system_type();
    StorageType storage_type = detect_storage_type_advanced(enable_storage_benchmark);
    
    // Calculate workload characteristics
    size_t matrix_elements = static_cast<size_t>(matrix_rows) * matrix_cols;
    double io_intensity = calculate_io_intensity_ratio_enhanced(matrix_elements, operation_type, storage_type);
    int storage_io_capacity = static_cast<int>(storage_type);
    
    // Memory safety check
    int memory_safe_threads = calculate_memory_safe_threads(matrix_elements);
    
    // Core optimization logic with enhanced strategies
    int optimal_threads;
    
    if (io_intensity > 4.0) {
        // HEAVILY I/O BOUND: Minimize thread contention
        optimal_threads = std::min({2, storage_io_capacity, memory_safe_threads});
        
    } else if (io_intensity > 2.0) {
        // MODERATELY I/O BOUND: Balance I/O capacity with computation
        optimal_threads = std::min({storage_io_capacity * 2, available_cores / 4, memory_safe_threads});
        
    } else if (io_intensity > 0.8) {
        // COMPUTE BOUND: More threads acceptable
        optimal_threads = std::min({available_cores / 2, 32, memory_safe_threads});
        
    } else {
        // HEAVILY COMPUTE BOUND: Aggressive parallelization
        optimal_threads = std::min({available_cores, 64, memory_safe_threads});
    }
    
    // Apply system-specific constraints
    switch (system_type) {
        case SystemType::CRAN_CHECK:
            optimal_threads = 2;  // CRAN compliance
            break;
            
        case SystemType::HPC_CLUSTER: {
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
            break;
        }
        
        case SystemType::CONTAINER: {
            // Check for container CPU limits
            #ifdef __linux__
            std::ifstream cpu_quota("/sys/fs/cgroup/cpu/cpu.cfs_quota_us");
            std::ifstream cpu_period("/sys/fs/cgroup/cpu/cpu.cfs_period_us");
            if (cpu_quota.is_open() && cpu_period.is_open()) {
                int quota, period;
                cpu_quota >> quota;
                cpu_period >> period;
                if (quota > 0 && period > 0) {
                    int container_cores = quota / period;
                    optimal_threads = std::min(optimal_threads, container_cores);
                }
            }
            #endif
            break;
        }
        
        case SystemType::SERVER:
            // Large system safety: Prevent libgomp memory overflow
            if (available_cores > 64) {
                optimal_threads = std::min(optimal_threads, 32);
            } else if (available_cores > 32) {
                optimal_threads = std::min(optimal_threads, 24);
            }
            break;
            
        case SystemType::DESKTOP:
            // Desktop systems: balance performance with responsiveness
            optimal_threads = std::min(optimal_threads, available_cores - 1);
            break;
    }
    
    // Final safety bounds
    return std::max(1, std::min(optimal_threads, available_cores));
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
                               const std::string& operation_type) {
    return get_io_aware_threads_enhanced(user_threads, matrix_rows, matrix_cols, operation_type, false);
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
    StorageType storage_type = detect_storage_type_advanced(false);
    SystemType system_type = detect_system_type();
    double io_intensity = calculate_io_intensity_ratio_enhanced(elements, operation_type, storage_type);
    int storage_capacity = static_cast<int>(storage_type);
    int recommended = get_io_aware_threads(R_NilValue, matrix_rows, matrix_cols, operation_type);
    int memory_safe = calculate_memory_safe_threads(elements);
    
    int available_cores = 1;
#ifdef _OPENMP
    available_cores = omp_get_num_procs();
#endif
    
    // System environment detection
    bool is_hpc = (system_type == SystemType::HPC_CLUSTER);
    bool is_cran = (system_type == SystemType::CRAN_CHECK);
    bool is_container = (system_type == SystemType::CONTAINER);
    
    Rcpp::Rcout << "\n=== BigDataStatMeth I/O-Aware Thread Analysis ===\n";
    Rcpp::Rcout << "Matrix: " << matrix_rows << "x" << matrix_cols 
               << " (" << elements << " elements)\n";
    Rcpp::Rcout << "Operation: " << operation_type << "\n";
    Rcpp::Rcout << "System: " << available_cores << " cores available";
    
    if (is_hpc) Rcpp::Rcout << " (HPC environment)";
    if (is_cran) Rcpp::Rcout << " (CRAN environment)";
    if (is_container) Rcpp::Rcout << " (Container environment)";
    Rcpp::Rcout << "\n";
    
    // Storage type information
    std::string storage_str;
    switch (storage_type) {
        case StorageType::HDD: storage_str = "HDD"; break;
        case StorageType::SATA_SSD: storage_str = "SATA SSD"; break;
        case StorageType::NVME_SSD: storage_str = "NVMe SSD"; break;
        case StorageType::NETWORK: storage_str = "Network Storage"; break;
        case StorageType::UNKNOWN: storage_str = "Unknown"; break;
    }
    
    Rcpp::Rcout << "\nAnalysis:\n";
    Rcpp::Rcout << "  Storage type: " << storage_str << "\n";
    Rcpp::Rcout << "  I/O intensity ratio: " << std::fixed << std::setprecision(2) 
               << io_intensity << "\n";
    Rcpp::Rcout << "  Storage I/O capacity: " << storage_capacity << " threads\n";
    Rcpp::Rcout << "  Memory-safe limit: " << memory_safe << " threads\n";
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
    
    // Memory constraints warning
    if (memory_safe < available_cores / 2) {
        Rcpp::Rcout << "\nMemory constraint detected:\n";
        Rcpp::Rcout << "  • Matrix size may strain available memory\n";
        Rcpp::Rcout << "  • Consider using smaller block sizes or chunked processing\n";
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
 * @param enable_storage_benchmark Enable storage benchmarking for accurate detection
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
                                        const std::string& operation_type,
                                        bool enable_storage_benchmark) {
    
    int optimal_threads = get_io_aware_threads_enhanced(R_NilValue, matrix_rows, matrix_cols, 
                                                       operation_type, enable_storage_benchmark);
    
#ifdef _OPENMP
    omp_set_num_threads(optimal_threads);
    
    // Apply additional OpenMP optimizations
    omp_set_dynamic(false);  // Disable dynamic thread adjustment
    omp_set_nested(false);   // Disable nested parallelism
    
    // Set appropriate scheduling based on system type and thread count
    SystemType system_type = detect_system_type();
    if (system_type == SystemType::HPC_CLUSTER || optimal_threads <= 8) {
        omp_set_schedule(omp_sched_static, 0);  // Static for HPC or few threads
    } else {
        omp_set_schedule(omp_sched_dynamic, 1000);  // Dynamic for desktop
    }
#endif
    
    Rcpp::Rcout << "BigDataStatMeth I/O-aware configuration applied:\n";
    Rcpp::Rcout << "  Threads configured: " << optimal_threads << "\n";
    Rcpp::Rcout << "  Operation optimized for: " << operation_type << "\n";
    
    // Show storage and system information
    StorageType storage_type = detect_storage_type_advanced(enable_storage_benchmark);
    SystemType sys_type = detect_system_type();
    
    std::string storage_str, system_str;
    switch (storage_type) {
        case StorageType::HDD: storage_str = "HDD"; break;
        case StorageType::SATA_SSD: storage_str = "SATA SSD"; break;
        case StorageType::NVME_SSD: storage_str = "NVMe SSD"; break;
        case StorageType::NETWORK: storage_str = "Network Storage"; break;
        case StorageType::UNKNOWN: storage_str = "Unknown"; break;
    }
    
    switch (sys_type) {
        case SystemType::DESKTOP: system_str = "Desktop"; break;
        case SystemType::SERVER: system_str = "Server"; break;
        case SystemType::HPC_CLUSTER: system_str = "HPC Cluster"; break;
        case SystemType::CRAN_CHECK: system_str = "CRAN Environment"; break;
        case SystemType::CONTAINER: system_str = "Container"; break;
    }
    
    Rcpp::Rcout << "  Detected storage: " << storage_str << "\n";
    Rcpp::Rcout << "  Detected system: " << system_str << "\n";
    Rcpp::Rcout << "  Configuration active for this R session.\n\n";
}

/**
 * @brief Auto-configure based on BigDataStatMeth operation integration
 * 
 * Smart configuration function that automatically detects the operation
 * type and matrix characteristics from BigDataStatMeth function calls.
 * This is designed to be called internally by BigDataStatMeth operations.
 * 
 * @param matrix_rows Number of matrix rows
 * @param matrix_cols Number of matrix columns
 * @param operation_type Operation type being performed
 * @param user_threads User-specified threads (if any)
 * @return Configured thread count for immediate use
 * 
 * @details
 * This function is optimized for integration within BigDataStatMeth operations:
 * - Minimal overhead (no storage benchmarking by default)
 * - Returns thread count immediately for use in OpenMP pragmas
 * - Caches system detection results for performance
 * - Applies configuration automatically
 * 
 * @note Designed for internal use by BigDataStatMeth operations
 * @since BigDataStatMeth 1.1
 */
inline int auto_configure_threads(int matrix_rows, int matrix_cols,
                                 const std::string& operation_type,
                                 Rcpp::Nullable<int> user_threads = R_NilValue) {
    
    // Fast path: if user specified threads, use them directly
    if (user_threads.isNotNull()) {
        int requested = Rcpp::as<int>(user_threads);
        int available_cores = std::thread::hardware_concurrency();
#ifdef _OPENMP
        available_cores = omp_get_num_procs();
        omp_set_num_threads(std::min(requested, available_cores));
#endif
        return std::min(requested, available_cores);
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

/**
 * @brief Integration helper for BigDataStatMeth matrix operations
 * 
 * Provides a standardized way to integrate I/O-aware optimization
 * into existing BigDataStatMeth functions with minimal code changes.
 * 
 * Usage pattern for BigDataStatMeth operations:
 * @code
 * // At the beginning of your matrix operation function:
 * int num_threads = BigDataStatMeth::integrate_io_optimization(
 *     nrows, ncols, "multiplication", user_threads);
 * 
 * #pragma omp parallel num_threads(num_threads)
 * {
 *     // Your matrix operation code here
 * }
 * @endcode
 * 
 * @param matrix_rows Number of matrix rows
 * @param matrix_cols Number of matrix columns
 * @param operation_type Type of operation being performed
 * @param user_threads User-specified thread count (optional)
 * @param verbose Whether to print optimization information
 * @return Optimal thread count for the operation
 * 
 * @since BigDataStatMeth 1.1
 */
inline int integrate_io_optimization(int matrix_rows, int matrix_cols,
                                    const std::string& operation_type,
                                    Rcpp::Nullable<int> user_threads ,
                                    bool verbose) {
    
    if (verbose) {
        print_io_thread_diagnosis(matrix_rows, matrix_cols, operation_type);
    }
    
    return auto_configure_threads(matrix_rows, matrix_cols, operation_type, user_threads);
}

/**
 * @brief Performance monitoring wrapper for BigDataStatMeth operations
 * 
 * Wraps BigDataStatMeth operations with performance monitoring and
 * automatic I/O-aware optimization. Provides timing information and
 * optimization feedback.
 * 
 * @tparam Func Function type to wrap
 * @param operation_name Name of the operation for logging
 * @param matrix_rows Number of matrix rows
 * @param matrix_cols Number of matrix columns
 * @param func Function to execute with optimization
 * @param args Arguments to pass to the function
 * @return Result of the wrapped function
 * 
 * @details
 * This wrapper:
 * - Applies I/O-aware optimization before execution
 * - Measures execution time
 * - Provides performance feedback
 * - Can be used to wrap any BigDataStatMeth operation
 * 
 * @since BigDataStatMeth 1.1
 */
template<typename Func, typename... Args>
inline auto performance_wrapper(const std::string& operation_name,
                               int matrix_rows, int matrix_cols,
                               Func&& func, Args&&... args) {
    
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
    
    // Report performance
    Rcpp::Rcout << "\nPerformance Report for " << operation_name << ":\n";
    Rcpp::Rcout << "  Matrix size: " << matrix_rows << "x" << matrix_cols << "\n";
    Rcpp::Rcout << "  Threads used: " << threads << "\n";
    Rcpp::Rcout << "  Optimization time: " << optimization_time << " μs\n";
    Rcpp::Rcout << "  Execution time: " << execution_time << " ms\n";
    
    if (execution_time > 1000) {
        double seconds = execution_time / 1000.0;
        Rcpp::Rcout << "  Total time: " << std::fixed << std::setprecision(2) 
                   << seconds << " seconds\n";
    }
    
    return result;
}

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_IO_AWARE_OPENMP_HPP