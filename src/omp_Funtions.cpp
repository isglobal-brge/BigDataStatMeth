#include <BigDataStatMeth.hpp>

//' Get BigDataStatMeth system information and OpenMP configuration
//'
//' This function provides detailed information about the current system
//' and OpenMP configuration, including automatic detection of HPC clusters,
//' Linux servers, and other system types.
//'
//' @return List containing system information:
//' \itemize{
//'   \item system_type - Type of system detected
//'   \item os_name - Operating system name
//'   \item num_cores - Number of CPU cores
//'   \item recommended_threads - Recommended number of threads
//'   \item current_threads - Currently configured threads
//'   \item is_hpc - Whether HPC environment was detected
//'   \item is_linux_server - Whether Linux server was detected
//'   \item is_conservative_mode - Whether conservative mode is active
//'   \item detection_method - Method used for system detection
//'   \item openmp_schedule - Recommended OpenMP schedule
//' }
//' @examples
//' # Get system information
//' info <- bdGetSystemInfo()
//' print(info)
//' 
//' # Check if we're on an HPC system
//' if (info$is_hpc) {
//'   message("HPC environment detected - using smart optimization")
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List bdGetSystemInfo()
{
   try {
       const SystemInfo& sys_info = get_system_info();
       
       return Rcpp::List::create(
           Rcpp::Named("system_type") = sys_info.system_type,
           Rcpp::Named("os_name") = sys_info.os_name,
           Rcpp::Named("num_cores") = sys_info.num_cores,
           Rcpp::Named("recommended_threads") = sys_info.recommended_threads,
           Rcpp::Named("current_threads") = getDTthreads(INT_MAX, false),
           Rcpp::Named("is_hpc") = sys_info.is_hpc,
           Rcpp::Named("is_linux_server") = sys_info.is_linux_server,
           Rcpp::Named("is_conservative_mode") = sys_info.is_conservative_mode,
           Rcpp::Named("detection_method") = sys_info.detection_method,
           Rcpp::Named("openmp_schedule") = BigDataStatMeth::get_recommended_omp_schedule()
       );
       
   } catch(std::exception& ex) {
       Rcpp::Rcerr << "Error in bdGetSystemInfo: " << ex.what() << std::endl;
       return R_NilValue;
   }
}

//' Configure OpenMP threads with intelligent optimization
//'
//' Allows users to specify the exact number of threads while applying
//' system-specific optimizations for best performance.
//'
//' @param threads Number of threads to use (required)
//' @param apply_optimizations Whether to apply system-specific optimizations (default: TRUE)
//' @examples
//' # On a 128-core server, use only 32 threads with optimizations
//' bdConfigureOpenMP(32)
//' 
//' # Use 48 threads without system optimizations
//' bdConfigureOpenMP(48, apply_optimizations = FALSE)
//' 
//' # Check what was configured
//' bdGetSystemInfo()
//' @export
// [[Rcpp::export]]
void bdConfigureOpenMP(int threads, bool apply_optimizations = true)
{
   try {
       const SystemInfo& sys_info = get_system_info();
       
       if (threads <= 0 || threads > sys_info.num_cores * 2) {
           Rcpp::Rcerr << "Invalid number of threads: " << threads << std::endl;
           Rcpp::Rcerr << "Must be between 1 and " << (sys_info.num_cores * 2) << std::endl;
           return;
       }
       
#ifdef _OPENMP
       omp_set_num_threads(threads);
       
       if (apply_optimizations) {
           apply_hpc_optimized_settings(threads);
           
           Rcpp::Rcout << "OpenMP configured with optimizations:\n";
           Rcpp::Rcout << "  System: " << sys_info.system_type << "\n";
           Rcpp::Rcout << "  Available cores: " << sys_info.num_cores << "\n";
           Rcpp::Rcout << "  Configured threads: " << threads << "\n";
           
           if (sys_info.is_hpc || sys_info.is_linux_server) {
               Rcpp::Rcout << "  Optimizations applied:\n";
               Rcpp::Rcout << "    ✓ Static scheduling\n";
               Rcpp::Rcout << "    ✓ Disabled dynamic adjustment\n";
               Rcpp::Rcout << "    ✓ Thread affinity optimization\n";
           }
           
       } else {
           Rcpp::Rcout << "OpenMP configured without optimizations:\n";
           Rcpp::Rcout << "  Threads: " << threads << "\n";
           Rcpp::Rcout << "  Using default OpenMP settings\n";
       }
       
       // Performance recommendation
       if (sys_info.is_linux_server && threads > sys_info.num_cores) {
           Rcpp::Rcout << "\n⚠️  Warning: Using more threads than physical cores may cause performance degradation\n";
       }
       
       if (threads > 64) {
           Rcpp::Rcout << "\n💡 Tip: For very high thread counts, consider NUMA-aware job allocation\n";
       }
       
#else
       Rcpp::Rcout << "OpenMP not available - thread configuration ignored\n";
#endif
       
   } catch(std::exception& ex) {
       Rcpp::Rcerr << "Error in bdConfigureOpenMP: " << ex.what() << std::endl;
   }
}

//' Apply intelligent HPC optimizations
//'
//' Optimizes OpenMP settings for HPC environments while respecting
//' the resources allocated by job schedulers (SLURM, PBS, etc.).
//' Uses ALL allocated threads with optimized configuration.
//'
//' @examples
//' # Apply HPC optimizations (respects allocated resources)
//' bdOptimizeForHPC()
//' 
//' # Check what was configured
//' bdGetSystemInfo()
//' @export
// [[Rcpp::export]]
void bdOptimizeForHPC()
{
   try {
       const SystemInfo& sys_info = get_system_info();
       int optimal_threads = get_optimal_threads(R_NilValue);
       
       apply_hpc_optimized_settings(optimal_threads);
       
       Rcpp::Rcout << "Smart HPC Optimization Applied:\n";
       Rcpp::Rcout << "  System: " << sys_info.system_type << "\n";
       Rcpp::Rcout << "  Threads: " << optimal_threads << " (respects job allocation)\n";
       
       if (sys_info.is_hpc) {
           Rcpp::Rcout << "  HPC Optimizations:\n";
           Rcpp::Rcout << "    ✓ Static scheduling (better for HPC)\n";
           Rcpp::Rcout << "    ✓ Disabled dynamic thread adjustment\n";
           Rcpp::Rcout << "    ✓ Consistent thread configuration\n";
           
           // Show job scheduler info if available
           const char* slurm_cpus = std::getenv("SLURM_CPUS_PER_TASK");
           const char* pbs_cpus = std::getenv("PBS_NCPUS");
           const char* job_id = std::getenv("SLURM_JOB_ID");
           if (!job_id) job_id = std::getenv("PBS_JOBID");
           
           if (job_id) {
               Rcpp::Rcout << "    ✓ Job ID: " << job_id << "\n";
           }
           if (slurm_cpus) {
               Rcpp::Rcout << "    ✓ SLURM CPUs allocated: " << slurm_cpus << "\n";
           }
           if (pbs_cpus) {
               Rcpp::Rcout << "    ✓ PBS CPUs allocated: " << pbs_cpus << "\n";
           }
       } else if (sys_info.is_linux_server) {
           Rcpp::Rcout << "  Linux Server Optimizations:\n";
           Rcpp::Rcout << "    ✓ Static scheduling\n";
           Rcpp::Rcout << "    ✓ Optimized for server workloads\n";
           Rcpp::Rcout << "    ✓ NUMA-aware configuration\n";
       }
       
       Rcpp::Rcout << "\nThis should provide optimal performance while using ALL allocated resources.\n";
       
   } catch(std::exception& ex) {
       Rcpp::Rcerr << "Error in bdOptimizeForHPC: " << ex.what() << std::endl;
   }
}

//' Print comprehensive system diagnostics
//'
//' Prints detailed diagnostic information about the current system,
//' OpenMP configuration, and performance settings. Useful for
//' troubleshooting and optimization.
//'
//' @examples
//' # Print full diagnostics
//' bdPrintDiagnostics()
//' @export
// [[Rcpp::export]]
void bdPrintDiagnostics()
{
   try {
       BigDataStatMeth::print_system_diagnostics();
       
       // Additional diagnostic information
       Rcpp::Rcout << "=== Environment Variables ===\n";
       
       const char* env_vars[] = {
           "OMP_NUM_THREADS", "OMP_SCHEDULE", "OMP_PROC_BIND", "OMP_PLACES",
           "SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "PBS_JOBID", "PBS_NCPUS",
           "R_DATATABLE_NUM_THREADS", "R_DATATABLE_NUM_PROCS_PERCENT",
           "BD_FORCE_CONSERVATIVE", "BD_SYSTEM_TYPE"
       };
       
       for (const char* var : env_vars) {
           const char* value = std::getenv(var);
           if (value) {
               Rcpp::Rcout << var << " = " << value << "\n";
           } else {
               Rcpp::Rcout << var << " = (not set)\n";
           }
       }
       
       Rcpp::Rcout << "===============================\n";
       
   } catch(std::exception& ex) {
       Rcpp::Rcerr << "Error in bdPrintDiagnostics: " << ex.what() << std::endl;
   }
}

//' Test OpenMP performance
//'
//' Performs a simple OpenMP performance test to measure the effectiveness
//' of current thread settings. Returns the execution time in milliseconds.
//'
//' @return Numeric value representing execution time in milliseconds
//' @examples
//' # Test current performance
//' time_ms <- bdTestOpenMPPerformance()
//' cat("OpenMP test completed in", time_ms, "milliseconds\n")
//' 
//' # Compare different configurations
//' bdConfigureOpenMP(32)
//' time_32 <- bdTestOpenMPPerformance()
//' 
//' bdConfigureOpenMP(64)
//' time_64 <- bdTestOpenMPPerformance()
//' 
//' cat("32 threads:", time_32, "ms\n")
//' cat("64 threads:", time_64, "ms\n")
//' @export
// [[Rcpp::export]]
double bdTestOpenMPPerformance()
{
   try {
       double time_ms = BigDataStatMeth::test_openmp_performance();
       
       const SystemInfo& sys_info = get_system_info();
       int current_threads = getDTthreads(INT_MAX, false);
       
       Rcpp::Rcout << "OpenMP Performance Test Results:\n";
       Rcpp::Rcout << "  System: " << sys_info.system_type << "\n";
       Rcpp::Rcout << "  Threads: " << current_threads << "\n";
       Rcpp::Rcout << "  Time: " << time_ms << " ms\n";
       
       if (sys_info.is_hpc || sys_info.is_linux_server) {
           Rcpp::Rcout << "  Mode: Optimized for " << sys_info.system_type << "\n";
       } else {
           Rcpp::Rcout << "  Mode: Standard\n";
       }
       
       return time_ms;
       
   } catch(std::exception& ex) {
       Rcpp::Rcerr << "Error in bdTestOpenMPPerformance: " << ex.what() << std::endl;
       return -1.0;
   }
}

//' Benchmark different thread configurations on current system
//'
//' Tests performance with different thread counts to find the optimal
//' configuration for the current workload and system.
//'
//' @param thread_counts Vector of thread counts to test (optional)
//' @param iterations Number of iterations per test (default: 3)
//' @examples
//' # Test common configurations
//' results <- bdBenchmarkThreads()
//' 
//' # Test specific thread counts
//' results <- bdBenchmarkThreads(c(16, 32, 48, 64))
//' print(results)
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame bdBenchmarkThreads(Rcpp::Nullable<Rcpp::IntegerVector> thread_counts = R_NilValue, 
                                 int iterations = 3)
{
   try {
       const SystemInfo& sys_info = get_system_info();
       
       std::vector<int> test_threads;
       
       if (thread_counts.isNotNull()) {
           Rcpp::IntegerVector user_threads = Rcpp::as<Rcpp::IntegerVector>(thread_counts);
           for (int i = 0; i < user_threads.size(); i++) {
               if (user_threads[i] > 0 && user_threads[i] <= sys_info.num_cores * 2) {
                   test_threads.push_back(user_threads[i]);
               }
           }
       } else {
           // Default test configurations based on system
           if (sys_info.num_cores >= 64) {
               test_threads = {1, 8, 16, 32, 48, 64};
           } else if (sys_info.num_cores >= 32) {
               test_threads = {1, 4, 8, 16, 24, 32};
           } else {
               test_threads = {1, 2, 4, 8, sys_info.num_cores};
           }
           
           // Remove duplicates and values > available cores
           std::sort(test_threads.begin(), test_threads.end());
           test_threads.erase(std::unique(test_threads.begin(), test_threads.end()), test_threads.end());
           test_threads.erase(std::remove_if(test_threads.begin(), test_threads.end(),
                             [&](int t) { return t > sys_info.num_cores; }), test_threads.end());
       }
       
       if (test_threads.empty()) {
           Rcpp::Rcerr << "No valid thread configurations to test" << std::endl;
           return R_NilValue;
       }
       
       Rcpp::Rcout << "Benchmarking thread configurations on " << sys_info.system_type << "...\n";
       Rcpp::Rcout << "Available cores: " << sys_info.num_cores << "\n\n";
       
       std::vector<int> results_threads;
       std::vector<double> results_times;
       std::vector<double> results_speedup;
       
       double baseline_time = 0.0;
       
       for (size_t i = 0; i < test_threads.size(); i++) {
           int t = test_threads[i];
           
           // Configure OpenMP for this test
           bdConfigureOpenMP(t, true); // Apply optimizations
           
           // Run benchmark
           double total_time = 0.0;
           for (int iter = 0; iter < iterations; iter++) {
               total_time += BigDataStatMeth::test_openmp_performance();
           }
           double avg_time = total_time / iterations;
           
           if (i == 0) {
               baseline_time = avg_time;
           }
           
           double speedup = baseline_time / avg_time;
           
           results_threads.push_back(t);
           results_times.push_back(avg_time);
           results_speedup.push_back(speedup);
           
           Rcpp::Rcout << "Threads: " << t << " | Time: " << avg_time << " ms | Speedup: " << speedup << "x\n";
       }
       
       // Find optimal configuration
       auto min_time_idx = std::min_element(results_times.begin(), results_times.end()) - results_times.begin();
       int optimal_threads = results_threads[min_time_idx];
       
       Rcpp::Rcout << "\n✅ Optimal configuration: " << optimal_threads << " threads\n";
       Rcpp::Rcout << "Speedup vs 1 thread: " << results_speedup[min_time_idx] << "x\n";
       
       // Apply optimal configuration
       bdConfigureOpenMP(optimal_threads, true);
       
       return Rcpp::DataFrame::create(
           Rcpp::Named("threads") = results_threads,
           Rcpp::Named("avg_time_ms") = results_times,
           Rcpp::Named("speedup") = results_speedup,
           Rcpp::Named("efficiency") = Rcpp::NumericVector(results_speedup.begin(), results_speedup.end()) / 
                                      Rcpp::NumericVector(results_threads.begin(), results_threads.end())
       );
       
   } catch(std::exception& ex) {
       Rcpp::Rcerr << "Error in bdBenchmarkThreads: " << ex.what() << std::endl;
       return R_NilValue;
   }
}

//' Get recommended thread configuration for current system
//'
//' Analyzes the current system and provides intelligent recommendations
//' for optimal thread configuration based on available resources.
//'
//' @param max_threads Maximum threads user wants to consider (optional)
//' @return List with recommendations
//' @examples
//' # Get recommendations for current system
//' rec <- bdGetRecommendations()
//' print(rec)
//' 
//' # Get recommendations limiting to 48 threads
//' rec <- bdGetRecommendations(max_threads = 48)
//' @export
// [[Rcpp::export]]
Rcpp::List bdGetRecommendations(Rcpp::Nullable<int> max_threads = R_NilValue)
{
   try {
       const SystemInfo& sys_info = get_system_info();
       
       int user_max = sys_info.num_cores;
       if (max_threads.isNotNull()) {
           user_max = std::min(Rcpp::as<int>(max_threads), sys_info.num_cores);
       }
       
       std::vector<std::string> recommendations;
       
       // System-specific recommendations
       if (sys_info.is_hpc) {
           recommendations.push_back("HPC detected: Use threads allocated by job scheduler");
           recommendations.push_back("Set OMP_NUM_THREADS to match SLURM_CPUS_PER_TASK or PBS_NCPUS");
           recommendations.push_back("Enable thread affinity: OMP_PROC_BIND=close");
           
       } else if (sys_info.is_linux_server) {
           if (sys_info.num_cores > 64) {
               recommendations.push_back("Large server detected: Consider NUMA topology");
               recommendations.push_back("Optimal range: " + std::to_string(sys_info.num_cores/4) + 
                                       " to " + std::to_string(sys_info.num_cores/2) + " threads");
           } else {
               recommendations.push_back("Server detected: Can use most available cores");
               recommendations.push_back("Optimal range: " + std::to_string(sys_info.num_cores/2) + 
                                       " to " + std::to_string(sys_info.num_cores) + " threads");
           }
           recommendations.push_back("Enable static scheduling for consistent performance");
           
       } else {
           recommendations.push_back("Desktop system: Can use most cores efficiently");
           recommendations.push_back("Recommended: " + std::to_string(sys_info.recommended_threads) + " threads");
       }
       
       // Performance considerations
       if (user_max > 32) {
           recommendations.push_back("For >32 threads: Monitor for diminishing returns");
           recommendations.push_back("Consider memory bandwidth limitations");
       }
       
       // Optimal configurations
       std::vector<int> optimal_configs;
       if (sys_info.is_linux_server && sys_info.num_cores > 64) {
           // Large server: provide several options
           optimal_configs.push_back(std::min(user_max, sys_info.num_cores / 4));
           optimal_configs.push_back(std::min(user_max, sys_info.num_cores / 2));
           optimal_configs.push_back(std::min(user_max, sys_info.num_cores * 3 / 4));
       } else {
           // Standard recommendations
           optimal_configs.push_back(std::min(user_max, sys_info.recommended_threads));
           if (sys_info.recommended_threads < user_max) {
               optimal_configs.push_back(std::min(user_max, sys_info.num_cores));
           }
       }
       
       return Rcpp::List::create(
           Rcpp::Named("system_type") = sys_info.system_type,
           Rcpp::Named("total_cores") = sys_info.num_cores,
           Rcpp::Named("user_max_threads") = user_max,
           Rcpp::Named("optimal_configs") = optimal_configs,
           Rcpp::Named("recommendations") = recommendations,
           Rcpp::Named("current_threads") = getDTthreads(INT_MAX, false)
       );
       
   } catch(std::exception& ex) {
       Rcpp::Rcerr << "Error in bdGetRecommendations: " << ex.what() << std::endl;
       return R_NilValue;
   }
}

//' Check CRAN compliance and thread limitations
//'
//' This function checks if the current environment has CRAN-like restrictions
//' and shows how threading is being limited for compliance.
//'
//' @return List with CRAN compliance information
//' @examples
//' # Check CRAN compliance status
//' cran_info <- bdCheckCRANCompliance()
//' print(cran_info)
//' 
//' # See if threads are being limited
//' if (cran_info$is_cran_environment) {
//'   message("Running in CRAN-like environment - threads limited to 2")
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List bdCheckCRANCompliance()
{
    try {
        const SystemInfo& sys_info = get_system_info();
        bool is_cran = is_cran_environment();
        
        // Get what threads would be without CRAN limits
        int unlimited_threads = sys_info.num_cores;
        if (sys_info.is_hpc) {
            const char* slurm_cpus = std::getenv("SLURM_CPUS_PER_TASK");
            if (slurm_cpus) {
                unlimited_threads = std::atoi(slurm_cpus);
            }
        }
        
        int actual_threads = get_optimal_threads(R_NilValue);
        bool is_limited = (actual_threads < unlimited_threads);
        
        std::vector<std::string> limitations;
        std::vector<std::string> recommendations;
        
        if (is_cran) {
            limitations.push_back("CRAN check environment detected");
            limitations.push_back("Threads automatically limited to 2 for compliance");
            recommendations.push_back("This is normal for package checking");
            recommendations.push_back("Full performance available in user environments");
        } else if (is_limited) {
            const char* omp_limit = std::getenv("OMP_THREAD_LIMIT");
            if (omp_limit) {
                limitations.push_back("OMP_THREAD_LIMIT=" + std::string(omp_limit) + " detected");
                recommendations.push_back("Remove OMP_THREAD_LIMIT for full performance");
            }
        } else {
            recommendations.push_back("No CRAN limitations detected");
            recommendations.push_back("Full threading performance available");
        }
        
        return Rcpp::List::create(
            Rcpp::Named("is_cran_environment") = is_cran,
            Rcpp::Named("is_thread_limited") = is_limited,
            Rcpp::Named("available_cores") = sys_info.num_cores,
            Rcpp::Named("unlimited_threads") = unlimited_threads,
            Rcpp::Named("actual_threads") = actual_threads,
            Rcpp::Named("limitations") = limitations,
            Rcpp::Named("recommendations") = recommendations,
            Rcpp::Named("system_type") = sys_info.system_type
        );
        
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdCheckCRANCompliance: " << ex.what() << std::endl;
        return R_NilValue;
    }
}

//' Override CRAN limitations (for user environments only)
//'
//' This function allows users to override CRAN limitations in their own
//' environments while maintaining compliance during package checking.
//' 
//' @param force_threads Number of threads to force (optional)
//' @param ignore_cran_limits Whether to ignore CRAN-like limitations (default: FALSE)
//' @examples
//' # Check current limitations
//' bdCheckCRANCompliance()
//' 
//' # In user environment, override for full performance
//' # (This should NOT be used during package checking)
//' bdOverrideCRANLimits(ignore_cran_limits = TRUE)
//' 
//' # Force specific number of threads
//' bdOverrideCRANLimits(force_threads = 64, ignore_cran_limits = TRUE)
//' @export
// [[Rcpp::export]]
void bdOverrideCRANLimits(Rcpp::Nullable<int> force_threads = R_NilValue, 
                          bool ignore_cran_limits = false)
{
    try {
        bool is_cran = is_cran_environment();
        
        if (is_cran && !ignore_cran_limits) {
            Rcpp::Rcout << "CRAN check environment detected.\n";
            Rcpp::Rcout << "Use ignore_cran_limits=TRUE only in user environments.\n";
            Rcpp::Rcout << "Current threads limited to 2 for CRAN compliance.\n";
            return;
        }
        
        if (ignore_cran_limits) {
            Rcpp::Rcout << "⚠️  WARNING: Ignoring CRAN limitations.\n";
            Rcpp::Rcout << "This should only be used in user environments, not during package checking.\n\n";
        }
        
        const SystemInfo& sys_info = get_system_info();
        int target_threads;
        
        if (force_threads.isNotNull()) {
            target_threads = Rcpp::as<int>(force_threads);
            target_threads = std::min(target_threads, sys_info.num_cores * 2); // Sanity limit
        } else {
            // Use full system capability
            if (sys_info.is_hpc) {
                const char* slurm_cpus = std::getenv("SLURM_CPUS_PER_TASK");
                if (slurm_cpus) {
                    target_threads = std::atoi(slurm_cpus);
                } else {
                    target_threads = sys_info.num_cores;
                }
            } else {
                target_threads = sys_info.num_cores;
            }
        }
        
        // Don't apply CRAN limits if override is requested
        if (!ignore_cran_limits) {
            target_threads = apply_cran_limits(target_threads);
        }
        
#ifdef _OPENMP
        omp_set_num_threads(target_threads);
        apply_hpc_optimized_settings(target_threads);
#endif
        
        Rcpp::Rcout << "Thread configuration updated:\n";
        Rcpp::Rcout << "  System: " << sys_info.system_type << "\n";
        Rcpp::Rcout << "  Threads: " << target_threads << "\n";
        
        if (ignore_cran_limits && is_cran) {
            Rcpp::Rcout << "  CRAN limits: IGNORED (user override)\n";
        } else if (!ignore_cran_limits && is_cran) {
            Rcpp::Rcout << "  CRAN limits: APPLIED (2 threads max)\n";
        } else {
            Rcpp::Rcout << "  CRAN limits: NOT APPLICABLE\n";
        }
        
    } catch(std::exception& ex) {
        Rcpp::Rcerr << "Error in bdOverrideCRANLimits: " << ex.what() << std::endl;
    }
}