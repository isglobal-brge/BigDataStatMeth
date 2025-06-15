/**
 * @file io_diagnostics.hpp
 * @brief Performance diagnostics and testing functions for BigDataStatMeth
 * 
 * Header-only diagnostic functions that provide comprehensive testing,
 * benchmarking, and analysis capabilities for both C++ and R users.
 */

#ifndef BIGDATASTATMETH_IO_DIAGNOSTICS_HPP
#define BIGDATASTATMETH_IO_DIAGNOSTICS_HPP

namespace BigDataStatMeth {

/**
 * @brief C++ performance testing function
 * 
 * Performs comprehensive performance testing with different thread configurations.
 * Header-only version for C++ developers.
 * 
 * @param test_size Size of test matrix
 * @param operation_type Type of operation to test
 * @param iterations Number of test iterations
 * @return PerformanceResult structure with detailed results
 */
inline PerformanceResult test_performance_cpp(int test_size, 
                                              const std::string& operation_type, 
                                              int iterations) {
    
    PerformanceResult result;
    result.operation_type = operation_type;
    result.test_size = test_size;
    result.iterations = iterations;
    result.times_ms.reserve(iterations);
    
    int current_threads = 1;
#ifdef _OPENMP
    current_threads = omp_get_max_threads();
#endif
    result.threads_used = current_threads;
    
    for (int i = 0; i < iterations; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        
        if (operation_type == "multiplication") {
            // Simple matrix multiplication benchmark
            std::vector<std::vector<double>> A(test_size, std::vector<double>(test_size, 1.0));
            std::vector<std::vector<double>> B(test_size, std::vector<double>(test_size, 1.0));
            std::vector<std::vector<double>> C(test_size, std::vector<double>(test_size, 0.0));
            
#ifdef _OPENMP
#pragma omp parallel for num_threads(current_threads)
#endif
            for (int row = 0; row < test_size; ++row) {
                for (int col = 0; col < test_size; ++col) {
                    double sum = 0.0;
                    for (int k = 0; k < test_size; ++k) {
                        sum += A[row][k] * B[k][col];
                    }
                    C[row][col] = sum;
                }
            }
            
            // Force computation to complete
            volatile double checksum = C[0][0];
            (void)checksum;
            
        } else if (operation_type == "crossprod") {
            // Cross-product benchmark
            std::vector<std::vector<double>> A(test_size, std::vector<double>(test_size, 1.0));
            std::vector<std::vector<double>> C(test_size, std::vector<double>(test_size, 0.0));
            
#ifdef _OPENMP
#pragma omp parallel for num_threads(current_threads)
#endif
            for (int i = 0; i < test_size; ++i) {
                for (int j = 0; j < test_size; ++j) {
                    double sum = 0.0;
                    for (int k = 0; k < test_size; ++k) {
                        sum += A[k][i] * A[k][j];
                    }
                    C[i][j] = sum;
                }
            }
            
            volatile double checksum = C[0][0];
            (void)checksum;
            
        } else if (operation_type == "svd") {
            // Simplified SVD-like computation
            std::vector<std::vector<double>> A(test_size, std::vector<double>(test_size, 1.0));
            std::vector<double> v(test_size, 1.0);
            
            // Simplified power iteration
            for (int iter = 0; iter < 10; ++iter) {
                std::vector<double> new_v(test_size, 0.0);
                
#ifdef _OPENMP
#pragma omp parallel for num_threads(current_threads)
#endif
                for (int i = 0; i < test_size; ++i) {
                    for (int j = 0; j < test_size; ++j) {
                        new_v[i] += A[i][j] * v[j];
                    }
                }
                
                // Normalize
                double norm = 0.0;
                for (int i = 0; i < test_size; ++i) {
                    norm += new_v[i] * new_v[i];
                }
                norm = std::sqrt(norm);
                
                for (int i = 0; i < test_size; ++i) {
                    v[i] = new_v[i] / norm;
                }
            }
            
            volatile double checksum = v[0];
            (void)checksum;
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        result.times_ms.push_back(static_cast<double>(duration.count()));
    }
    
    // Calculate statistics
    result.mean_time_ms = 0.0;
    for (double t : result.times_ms) {
        result.mean_time_ms += t;
    }
    result.mean_time_ms /= result.times_ms.size();
    
    result.min_time_ms = *std::min_element(result.times_ms.begin(), result.times_ms.end());
    result.max_time_ms = *std::max_element(result.times_ms.begin(), result.times_ms.end());
    
    // Calculate standard deviation
    double variance = 0.0;
    for (double t : result.times_ms) {
        variance += (t - result.mean_time_ms) * (t - result.mean_time_ms);
    }
    variance /= result.times_ms.size();
    result.std_dev_ms = std::sqrt(variance);
    result.relative_std_dev = result.std_dev_ms / result.mean_time_ms;
    
    // Performance metrics
    if (operation_type == "multiplication") {
        // Approximate FLOPS for matrix multiplication
        double flops = 2.0 * test_size * test_size * test_size;
        result.ops_per_second = flops / (result.mean_time_ms / 1000.0);
    }
    
    return result;
}

/**
 * @brief C++ auto-tuning function
 * 
 * Performs empirical testing to find optimal thread configuration.
 * Header-only version for C++ developers.
 * 
 * @param sample_rows Number of rows for testing
 * @param sample_cols Number of columns for testing
 * @param operation_type Primary operation type to optimize for
 * @param max_threads Maximum threads to test
 * @param test_iterations Number of iterations per thread count
 * @return Rcpp::List with optimal configuration and results
 */
inline std::pair<int, std::vector<PerformanceResult>> auto_tune_cpp(
    int sample_rows, int sample_cols,
    const std::string& operation_type,
    int max_threads, int test_iterations) {
    
    int available_cores = std::thread::hardware_concurrency();
#ifdef _OPENMP
    available_cores = omp_get_num_procs();
#endif
    
    if (max_threads == 0) {
        max_threads = available_cores;
    } else {
        max_threads = std::min(max_threads, available_cores);
    }
    
    std::vector<PerformanceResult> results;
    std::vector<double> performance_scores;
    
    // Store original configuration
    int original_threads = 1;
#ifdef _OPENMP
    original_threads = omp_get_max_threads();
#endif
    
    // Test different thread counts
    for (int threads = 1; threads <= max_threads; threads++) {
#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif
        
        PerformanceResult test_result = test_performance_cpp(
            std::min(sample_rows, 2000), operation_type, test_iterations);
        
        results.push_back(test_result);
        performance_scores.push_back(1000.0 / test_result.mean_time_ms);
    }
    
    // Find optimal thread count
    auto max_it = std::max_element(performance_scores.begin(), performance_scores.end());
    int optimal_threads = static_cast<int>(std::distance(performance_scores.begin(), max_it) + 1);
    
    // Apply optimal configuration
#ifdef _OPENMP
    omp_set_num_threads(optimal_threads);
    omp_set_dynamic(false);
#endif
    
    return std::make_pair(optimal_threads, results);
}

/**
 * @brief C++ system validation function
 * 
 * Comprehensive validation of the I/O-aware optimization system.
 * Header-only version for C++ developers.
 * 
 * @param quick_test Run only essential tests
 * @param include_benchmarks Include performance benchmarks
 * @return Map with validation results
 */
inline std::map<std::string, bool> validate_optimization_cpp(bool quick_test, 
                                                            bool include_benchmarks) {
    
    std::map<std::string, bool> validation_results;
    
    // Test 1: System Detection
    try {
        SystemInfo system_info = get_system_info_cpp();
        validation_results["system_detection"] = true;
        validation_results["openmp_available"] = system_info.openmp_enabled;
    } catch (...) {
        validation_results["system_detection"] = false;
        validation_results["openmp_available"] = false;
    }
    
    // Test 2: Thread Optimization Logic
    try {
        int threads1 = get_io_aware_threads(R_NilValue, 1000, 1000, "multiplication");
        int threads2 = get_io_aware_threads(R_NilValue, 5000, 5000, "svd");
        validation_results["optimization_logic"] = (threads1 > 0 && threads2 > 0);
    } catch (...) {
        validation_results["optimization_logic"] = false;
    }
    
    // Test 3: Configuration Application
    try {
        ConfigInfo original_config = get_config_info_cpp();
        configure_bigdata_io_threads(2000, 2000, "multiplication");
        ConfigInfo new_config = get_config_info_cpp();
        validation_results["configuration_test"] = true;
        validation_results["configuration_changed"] = 
            (new_config.current_threads != original_config.current_threads);
    } catch (...) {
        validation_results["configuration_test"] = false;
        validation_results["configuration_changed"] = false;
    }
    
    if (!quick_test) {
        // Test 4: I/O Intensity Analysis
        try {
            double intensity1 = calculate_io_intensity_ratio(1000000, "multiplication");
            double intensity2 = calculate_io_intensity_ratio(1000000, "svd");
            validation_results["intensity_analysis"] = (intensity1 > 0 && intensity2 > intensity1);
        } catch (...) {
            validation_results["intensity_analysis"] = false;
        }
    }
    
    if (include_benchmarks) {
        // Test 5: Performance Benchmarks
        try {
            PerformanceResult perf_mult = test_performance_cpp(1000, "multiplication", 2);
            PerformanceResult perf_cross = test_performance_cpp(1000, "crossprod", 2);
            validation_results["benchmarks"] = (perf_mult.mean_time_ms > 0 && perf_cross.mean_time_ms > 0);
        } catch (...) {
            validation_results["benchmarks"] = false;
        }
    }
    
    return validation_results;
}

/**
 * @brief Comprehensive I/O diagnosis with C++ interface
 * 
 * Provides detailed analysis and prints diagnostic information.
 * Can be called from both C++ and R.
 * 
 * @param matrix_rows Number of matrix rows
 * @param matrix_cols Number of matrix columns
 * @param operation_type Operation type for analysis
 * @param enable_storage_benchmark Enable storage benchmarking
 * @param print_output Whether to print diagnostic output
 */
inline void diagnose_io_cpp(int matrix_rows, int matrix_cols, 
                           const std::string& operation_type,
                           bool enable_storage_benchmark,
                           bool print_output) {
    
    size_t elements = static_cast<size_t>(matrix_rows) * matrix_cols;
    StorageType storage_type = detect_storage_type_advanced(enable_storage_benchmark);
    SystemType system_type = detect_system_type();
    double io_intensity = calculate_io_intensity_ratio_enhanced(elements, operation_type, storage_type);
    int storage_capacity = static_cast<int>(storage_type);
    int recommended = get_io_aware_threads(R_NilValue, matrix_rows, matrix_cols, operation_type);
    int memory_safe = calculate_memory_safe_threads(elements);
    
    int available_cores = 1;
#ifdef _OPENMP
    available_cores = omp_get_num_procs();
#endif
    
    if (!print_output) return;
    
    // Print detailed diagnosis
    std::cout << "\n=== BigDataStatMeth I/O-Aware Thread Analysis ===\n";
    std::cout << "Matrix: " << matrix_rows << "x" << matrix_cols 
              << " (" << elements << " elements)\n";
    std::cout << "Operation: " << operation_type << "\n";
    std::cout << "System: " << available_cores << " cores available";
    
    if (system_type == SystemType::HPC_CLUSTER) std::cout << " (HPC environment)";
    if (system_type == SystemType::CRAN_CHECK) std::cout << " (CRAN environment)";
    if (system_type == SystemType::CONTAINER) std::cout << " (Container environment)";
    std::cout << "\n";
    
    // Storage type information
    std::string storage_str;
    switch (storage_type) {
        case StorageType::HDD: storage_str = "HDD"; break;
        case StorageType::SATA_SSD: storage_str = "SATA SSD"; break;
        case StorageType::NVME_SSD: storage_str = "NVMe SSD"; break;
        case StorageType::NETWORK: storage_str = "Network Storage"; break;
        case StorageType::UNKNOWN: storage_str = "Unknown"; break;
    }
    
    std::cout << "\nAnalysis:\n";
    std::cout << "  Storage type: " << storage_str << "\n";
    std::cout << "  I/O intensity ratio: " << std::fixed << std::setprecision(2) 
              << io_intensity << "\n";
    std::cout << "  Storage I/O capacity: " << storage_capacity << " threads\n";
    std::cout << "  Memory-safe limit: " << memory_safe << " threads\n";
    std::cout << "  Recommended threads: " << recommended << "\n";
    
    // Explain the strategy
    std::cout << "\nOptimization strategy: ";
    if (io_intensity > 4.0) {
        std::cout << "I/O BOUND - Minimize thread contention\n";
        std::cout << "  → Using minimal threads to avoid I/O critical section bottlenecks\n";
    } else if (io_intensity > 2.0) {
        std::cout << "BALANCED - I/O and compute optimization\n";
        std::cout << "  → Balancing I/O capacity with computational parallelism\n";
    } else if (io_intensity > 0.8) {
        std::cout << "COMPUTE BOUND - CPU-focused parallelization\n";
        std::cout << "  → Using more threads while respecting I/O constraints\n";
    } else {
        std::cout << "CPU INTENSIVE - Aggressive parallelization\n";
        std::cout << "  → Maximum parallelization with libgomp safety limits\n";
    }
    
    // Performance expectations
    if (available_cores > 64 && recommended < available_cores/3) {
        std::cout << "\nNote: Large system detected (" << available_cores 
                  << " cores) - threads limited to prevent:\n";
        std::cout << "  • libgomp memory overflow\n";
        std::cout << "  • I/O critical section contention\n";
        std::cout << "  • Memory bandwidth saturation\n";
        
        double improvement = static_cast<double>(available_cores) / recommended * 1.5;
        std::cout << "  Expected performance improvement: " 
                  << std::fixed << std::setprecision(1) << improvement << "x faster\n";
    }
    
    // Memory constraints warning
    if (memory_safe < available_cores / 2) {
        std::cout << "\nMemory constraint detected:\n";
        std::cout << "  • Matrix size may strain available memory\n";
        std::cout << "  • Consider using smaller block sizes or chunked processing\n";
    }
    
    std::cout << "================================================\n\n";
}

/**
 * @brief Quick configuration check for C++ users
 * 
 * Returns current optimization status and configuration.
 * 
 * @return ConfigInfo structure with current configuration
 */
inline ConfigInfo check_config_cpp() {
    return get_config_info_cpp();
}

/**
 * @brief Performance comparison for C++ users
 * 
 * Compares performance with and without I/O optimization.
 * 
 * @param matrix_rows Test matrix rows
 * @param matrix_cols Test matrix columns
 * @param operation_type Operation type to test
 * @param iterations Number of test iterations
 * @return Pair of results: (unoptimized, optimized)
 */
inline std::pair<PerformanceResult, PerformanceResult> compare_performance_cpp(
    int matrix_rows, int matrix_cols, const std::string& operation_type,
    int iterations) {
    
    // Store original configuration
    ConfigInfo original_config = get_config_info_cpp();
    
    // Test 1: Without optimization (use all cores)
    int max_cores = std::thread::hardware_concurrency();
#ifdef _OPENMP
    max_cores = omp_get_num_procs();
    omp_set_num_threads(max_cores);
#endif
    
    PerformanceResult unoptimized = test_performance_cpp(
        std::min(matrix_rows, 2000), operation_type, iterations);
    
    // Test 2: With I/O optimization
    configure_bigdata_io_threads(matrix_rows, matrix_cols, operation_type);
    
    PerformanceResult optimized = test_performance_cpp(
        std::min(matrix_rows, 2000), operation_type, iterations);
    
    return std::make_pair(unoptimized, optimized);
}

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_IO_DIAGNOSTICS_HPP