/**
 * @file system_analysis.hpp
 * @brief Advanced system analysis and detection for BigDataStatMeth
 * 
 * Header-only system detection functions that analyze hardware,
 * storage, memory, and environment characteristics for optimal
 * thread configuration.
 */

#ifndef BIGDATASTATMETH_SYSTEM_ANALYSIS_HPP
#define BIGDATASTATMETH_SYSTEM_ANALYSIS_HPP

namespace BigDataStatMeth {

// /**
//  * @brief Advanced storage detection with benchmark capabilities
//  * 
//  * Enhanced storage detection that can optionally perform I/O benchmarks
//  * to more accurately determine storage characteristics.
//  * 
//  * @param enable_benchmark Whether to perform I/O benchmark test
//  * @param benchmark_size Size of benchmark test in KB
//  * @return StorageType enum indicating detected storage type
//  */
// inline StorageType detect_storage_type_advanced(bool enable_benchmark, 
//                                                 size_t benchmark_size) {
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
//             if (ofs.is_open()) {
//                 ofs.write(test_data.data(), test_data.size());
//                 ofs.flush();
//                 ofs.close();
//                 
//                 auto end = std::chrono::high_resolution_clock::now();
//                 auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//                 
//                 // Calculate throughput (MB/s)
//                 double mb_per_sec = (benchmark_size / 1024.0) / (duration.count() / 1e6);
//                 
//                 // Clean up
//                 std::remove(temp_file.c_str());
//                 
//                 // Classify based on throughput
//                 if (mb_per_sec > 500) {
//                     return StorageType::NVME_SSD;
//                 } else if (mb_per_sec > 100) {
//                     return StorageType::SATA_SSD;
//                 } else {
//                     return StorageType::HDD;
//                 }
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
//     if (mounts.is_open()) {
//         std::string line;
//         
//         // Check for high-performance storage indicators
//         while (std::getline(mounts, line)) {
//             if (line.find("nvme") != std::string::npos) {
//                 return StorageType::NVME_SSD;
//             } else if (line.find("ssd") != std::string::npos) {
//                 return StorageType::SATA_SSD;
//             } else if (line.find("tmpfs") != std::string::npos && 
//                 line.find("/tmp") != std::string::npos) {
//                 return StorageType::NVME_SSD;  // RAM disk, treat as fastest
//             }
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

// /**
//  * @brief Detect system type and environment characteristics
//  * 
//  * Comprehensive system environment detection that identifies the type
//  * of system and runtime environment to apply appropriate optimizations.
//  * 
//  * @return SystemType enum indicating detected system type
//  */
// inline SystemType detect_system_type() {
//     
//     // Check for CRAN environment first (highest priority)
//     const char* thread_limit = std::getenv("OMP_THREAD_LIMIT");
//     const char* r_check = std::getenv("_R_CHECK_PACKAGE_NAME_");
//     if ((thread_limit && std::atoi(thread_limit) <= 2) || r_check) {
//         return SystemType::CRAN_CHECK;
//     }
//     
//     // Check for HPC job scheduler environments
//     if (std::getenv("SLURM_JOB_ID") || std::getenv("PBS_JOBID") ||
//         std::getenv("LSB_JOBID") || std::getenv("SGE_JOB_ID") ||
//         std::getenv("SLURM_CPUS_PER_TASK") || std::getenv("PBS_NCPUS")) {
//         return SystemType::HPC_CLUSTER;
//     }
//     
//     // Check for container environments
// #ifdef __linux__
//     std::ifstream cgroup("/proc/1/cgroup");
//     if (cgroup.is_open()) {
//         std::string line;
//         while (std::getline(cgroup, line)) {
//             if (line.find("docker") != std::string::npos ||
//                 line.find("lxc") != std::string::npos ||
//                 line.find("containerd") != std::string::npos) {
//                 return SystemType::CONTAINER;
//             }
//         }
//     }
//     
//     // Check for Singularity/Apptainer
//     if (std::getenv("SINGULARITY_CONTAINER") || std::getenv("APPTAINER_CONTAINER")) {
//         return SystemType::CONTAINER;
//     }
// #endif
//     
//     // Determine based on core count
//     int available_cores = std::thread::hardware_concurrency();
// #ifdef _OPENMP
//     available_cores = omp_get_num_procs();
// #endif
//     
//     if (available_cores > 16) {
//         return SystemType::SERVER;
//     } else {
//         return SystemType::DESKTOP;
//     }
// }

// /**
//  * @brief Memory-aware thread limiting for large systems
//  * 
//  * Calculates safe thread limits based on available memory and matrix size
//  * to prevent memory exhaustion and excessive swapping.
//  * 
//  * @param matrix_elements Total number of matrix elements
//  * @param element_size Size of each matrix element in bytes
//  * @param safety_factor Memory safety factor
//  * @return Maximum threads considering memory constraints
//  */
// inline int calculate_memory_safe_threads(size_t matrix_elements, 
//                                          size_t element_size,
//                                          double safety_factor) {
//     
//     // Get available system memory
//     size_t available_memory = 0;
//     
// #ifdef __linux__
//     std::ifstream meminfo("/proc/meminfo");
//     if (meminfo.is_open()) {
//         std::string line;
//         while (std::getline(meminfo, line)) {
//             if (line.find("MemAvailable:") == 0) {
//                 size_t kb;
//                 sscanf(line.c_str(), "MemAvailable: %zu kB", &kb);
//                 available_memory = kb * 1024;  // Convert to bytes
//                 break;
//             }
//         }
//     }
// #endif
//     
//     // Fallback: estimate based on physical memory
//     if (available_memory == 0) {
//         // Conservative estimate: assume 75% of 8GB per 16 cores
//         int cores = std::thread::hardware_concurrency();
//         available_memory = static_cast<size_t>(cores / 16.0 * 8ULL * 1024 * 1024 * 1024 * 0.75);
//     }
//     
//     // Calculate memory requirements
//     size_t matrix_memory = matrix_elements * element_size;
//     size_t per_operation_memory = matrix_memory * 3;  // A + B + C matrices
//     size_t per_thread_overhead = 64 * 1024 * 1024;   // 64MB per thread
//     
//     size_t safe_memory = static_cast<size_t>(available_memory * safety_factor);
//     
//     if (per_operation_memory > safe_memory) {
//         return 1;  // Even single-threaded may be tight
//     }
//     
//     size_t remaining_memory = safe_memory - per_operation_memory;
//     int max_threads = static_cast<int>(remaining_memory / per_thread_overhead);
//     
//     return std::max(1, max_threads);
// }

/**
 * @brief Get comprehensive system information
 * 
 * Returns detailed information about system capabilities and configuration.
 * 
 * @param include_benchmark Include storage benchmark in analysis
 * @return SystemInfo structure with comprehensive system data
 */
inline SystemInfo get_system_info_cpp(bool include_benchmark) {
    
    SystemInfo info;
    
    // Basic system information
    int available_cores = std::thread::hardware_concurrency();
    int current_threads = 1;
    
#ifdef _OPENMP
    available_cores = omp_get_num_procs();
    current_threads = omp_get_max_threads();
#endif
    
    // Detect system characteristics
    SystemType sys_type = detect_system_type();
    StorageType storage_type = detect_storage_type_advanced(include_benchmark, 1024);
    
    // Fill system info structure
    info.version = BIGDATASTATMETH_VERSION_STRING;
    info.cpu_cores = available_cores;
    info.current_openmp_threads = current_threads;
    info.storage_io_capacity = static_cast<int>(storage_type);
    info.memory_safe_threads_5kx5k = calculate_memory_safe_threads(25000000, 8, 0.8);
    info.openmp_enabled = BIGDATASTATMETH_HAS_OPENMP;
    info.is_hpc = (sys_type == SystemType::HPC_CLUSTER);
    info.is_cran = (sys_type == SystemType::CRAN_CHECK);
    info.is_container = (sys_type == SystemType::CONTAINER);
    
    // Convert enums to strings
    switch (sys_type) {
    case SystemType::DESKTOP: info.system_type = "Desktop/Laptop"; break;
    case SystemType::SERVER: info.system_type = "Server"; break;
    case SystemType::HPC_CLUSTER: info.system_type = "HPC Cluster"; break;
    case SystemType::CRAN_CHECK: info.system_type = "CRAN Environment"; break;
    case SystemType::CONTAINER: info.system_type = "Container"; break;
    }
    
    switch (storage_type) {
    case StorageType::HDD: info.storage_type = "Hard Disk Drive"; break;
    case StorageType::SATA_SSD: info.storage_type = "SATA SSD"; break;
    case StorageType::NVME_SSD: info.storage_type = "NVMe SSD"; break;
    case StorageType::NETWORK: info.storage_type = "Network Storage"; break;
    case StorageType::UNKNOWN: info.storage_type = "Unknown"; break;
    }
    
    // Feature availability
    info.features.push_back("openmp");
    info.features.push_back("io_optimization");
    info.features.push_back("storage_detection");
    info.features.push_back("system_analysis");
    info.features.push_back("performance_diagnostics");
    
    // Memory information
#ifdef __linux__
    std::ifstream meminfo("/proc/meminfo");
    if (meminfo.is_open()) {
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.find("MemAvailable:") == 0) {
                size_t kb;
                sscanf(line.c_str(), "MemAvailable: %zu kB", &kb);
                info.available_memory_gb = kb / 1024.0 / 1024.0;
                break;
            }
        }
    }
#endif
    
    return info;
}

/**
 * @brief Get current configuration information
 * 
 * Returns detailed information about current OpenMP and system configuration.
 * 
 * @return ConfigInfo structure with current configuration data
 */
inline ConfigInfo get_config_info_cpp() {
    
    ConfigInfo config;
    
    int current_threads = 1;
    int available_cores = 1;
    
#ifdef _OPENMP
    current_threads = omp_get_max_threads();
    available_cores = omp_get_num_procs();
#endif
    
    // Enhanced system detection
    SystemType sys_type = detect_system_type();
    StorageType storage_type = detect_storage_type_advanced(false, 1024);
    
    // Fill configuration structure
    config.version = BIGDATASTATMETH_VERSION_STRING;
    config.current_threads = current_threads;
    config.available_cores = available_cores;
    config.is_hpc = (sys_type == SystemType::HPC_CLUSTER);
    config.is_cran = (sys_type == SystemType::CRAN_CHECK);
    config.is_container = (sys_type == SystemType::CONTAINER);
    config.optimization_active = (current_threads < available_cores) || 
        (current_threads <= 32 && available_cores > 64);
    config.thread_efficiency = static_cast<double>(current_threads) / available_cores;
    
    // Convert enums to strings
    switch (sys_type) {
    case SystemType::DESKTOP: config.system_type = "Desktop"; break;
    case SystemType::SERVER: config.system_type = "Server"; break;
    case SystemType::HPC_CLUSTER: config.system_type = "HPC Cluster"; break;
    case SystemType::CRAN_CHECK: config.system_type = "CRAN Environment"; break;
    case SystemType::CONTAINER: config.system_type = "Container"; break;
    }
    
    switch (storage_type) {
    case StorageType::HDD: config.storage_type = "HDD"; break;
    case StorageType::SATA_SSD: config.storage_type = "SATA SSD"; break;
    case StorageType::NVME_SSD: config.storage_type = "NVMe SSD"; break;
    case StorageType::NETWORK: config.storage_type = "Network Storage"; break;
    case StorageType::UNKNOWN: config.storage_type = "Unknown"; break;
    }
    
    // Collect environment variables
    const char* env_names[] = {
        "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "SLURM_CPUS_PER_TASK", 
        "PBS_NCPUS", "SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID", "SGE_JOB_ID",
        "SINGULARITY_CONTAINER", "APPTAINER_CONTAINER"
    };
    
    for (const char* env_name : env_names) {
        const char* env_value = std::getenv(env_name);
        if (env_value) {
            std::string entry = std::string(env_name) + "=" + std::string(env_value);
            config.environment_vars.push_back(entry);
        }
    }
    
    // Memory information
#ifdef __linux__
    std::ifstream meminfo("/proc/meminfo");
    if (meminfo.is_open()) {
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.find("MemAvailable:") == 0) {
                size_t kb;
                sscanf(line.c_str(), "MemAvailable: %zu kB", &kb);
                config.available_memory_kb = static_cast<double>(kb);
                break;
            }
        }
    }
#endif
    
    return config;
}

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_SYSTEM_ANALYSIS_HPP