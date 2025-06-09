/**
 * @file openme-utils.hpp
 * @brief OpenMP threading utilities and configuration management
 * with system detection
 *
 * This file provides comprehensive thread management utilities for OpenMP
 * parallelization in the BigDataStatMeth package. It handles thread
 * initialization, configuration, forking behavior, and runtime thread
 * count determination based on system resources and environment variables.
 * There are functions for system detection to optimize performance
 * on HPC clusters and Linux servers while maintaining full compatibility
 * with the existing BigDataStatMeth codebase. Automatically respects CRAN 
 * policies while allowing full performance in user environments.
 *
 * Key features:
 * - Dynamic thread count management
 * - Environment variable configuration
 * - Fork-safe thread handling
 * - Thread throttling for performance optimization
 * - System resource detection and limits
 * - Automatic HPC/Linux server detection
 * - Smart thread configuration respecting job allocations
 * - Environment-specific optimizations
 * - Performance monitoring and diagnostics
 *
 * Key improvements for HPC/Linux servers:
 * - Respects job scheduler allocations (SLURM, PBS)
 * - Optimized OpenMP configuration without limiting threads
 * - Static scheduling for better performance
 * - Thread affinity optimizations
 *
 * Environment variables:
 * - R_DATATABLE_NUM_THREADS: Direct thread count control
 * - R_DATATABLE_NUM_PROCS_PERCENT: CPU utilization percentage
 * - R_DATATABLE_THROTTLE: Thread activation threshold
 * - BD_FORCE_CONSERVATIVE: Force conservative mode
 * - BD_SYSTEM_TYPE: Override system detection
 * - SLURM_CPUS_PER_TASK: SLURM CPU allocation
 * - PBS_NCPUS: PBS CPU allocation
 *
 */

#include <Rcpp.h>
#include "H5Cpp.h"

#ifdef _OPENMP
#include <pthread.h>
#endif
#include <errno.h>     // errno
#include <ctype.h>     // isspace

#ifndef _INITVARS
#define _INITVARS

/**
 * @brief System information structure
 */
struct SystemInfo {
    std::string system_type;
    std::string os_name;
    int num_cores;
    int recommended_threads;
    bool is_hpc;
    bool is_linux_server;
    bool is_conservative_mode;
    std::string detection_method;
};

/**
 * @brief Global thread management variables (existing)
 * @{
 */
static int  DTthreads = -1;    //!< Current thread count, initialized to -1
static int  DTthrottle = -1;   //!< Thread activation threshold
static bool RestoreAfterFork = true; //!< Thread restoration behavior after fork
/** @} */

/**
 * @brief Global variables for system detection
 * @{
 */
static SystemInfo g_system_info;
static bool g_system_detected = false;
static bool g_conservative_forced = false;
/** @} */

/**
 * @brief Check if we're running in a CRAN-like environment
 * @return bool True if CRAN restrictions should be applied
 */
static bool is_cran_environment() {
    // CRAN sets OMP_THREAD_LIMIT=2
    const char* thread_limit = std::getenv("OMP_THREAD_LIMIT");
    if (thread_limit && std::atoi(thread_limit) <= 2) {
        return true;
    }
    
    // Check for CRAN-specific environment variables
    if (std::getenv("_R_CHECK_EXAMPLE_TIMING_") != nullptr ||
        std::getenv("_R_CHECK_TESTS_NLINES_") != nullptr ||
        std::getenv("R_CMD_CHECK") != nullptr) {
        return true;
    }
    
    // Check for typical CRAN server patterns
    const char* hostname = std::getenv("HOSTNAME");
    if (hostname) {
        std::string host_str(hostname);
        std::transform(host_str.begin(), host_str.end(), host_str.begin(), ::tolower);
        if (host_str.find("cran") != std::string::npos ||
            host_str.find("check") != std::string::npos ||
            host_str.find("winbuilder") != std::string::npos) {
            return true;
        }
    }
    
    return false;
}

/**
 * @brief Apply CRAN-compliant thread limits
 * @param requested_threads Number of threads user/system wants
 * @return int CRAN-compliant number of threads
 */
static int apply_cran_limits(int requested_threads) {
    if (is_cran_environment()) {
        // CRAN compliance: maximum 2 threads
        return std::min(requested_threads, 2);
    }
    
    // Check OMP_THREAD_LIMIT (CRAN sets this to 2)
    const char* omp_limit = std::getenv("OMP_THREAD_LIMIT");
    if (omp_limit) {
        int limit = std::atoi(omp_limit);
        if (limit > 0) {
            return std::min(requested_threads, limit);
        }
    }
    
    // In user environments, no artificial limits
    return requested_threads;
}


/**
 * @brief Detects the current system type and characteristics (SMART VERSION)
 * @return SystemInfo Structure with system details
 */
static SystemInfo detect_system_type() {
    SystemInfo info;
    
    // Initialize defaults
    info.num_cores = 1;
    info.is_hpc = false;
    info.is_linux_server = false;
    info.is_conservative_mode = false;
    info.detection_method = "default";
    
#ifdef _OPENMP
    info.num_cores = omp_get_num_procs();
#endif
    
    // Detect OS
#ifdef _WIN32
    info.os_name = "Windows";
#elif defined(__APPLE__)
    info.os_name = "macOS";
#elif defined(__linux__)
    info.os_name = "Linux";
#else
    info.os_name = "Unknown";
#endif
    
    // Check for CRAN environment first
    bool is_cran = is_cran_environment();
    if (is_cran) {
        info.system_type = "CRAN Check Environment";
        info.recommended_threads = 2; // CRAN limit
        info.detection_method = "CRAN compliance";
        info.is_conservative_mode = true; // Force conservative for CRAN
        return info;
    }
    
    // Check for forced conservative mode
    const char* force_conservative = std::getenv("BD_FORCE_CONSERVATIVE");
    if (force_conservative && (std::strcmp(force_conservative, "1") == 0 || 
        std::strcmp(force_conservative, "true") == 0)) {
        info.is_conservative_mode = true;
        info.system_type = "Forced Conservative";
        info.recommended_threads = apply_cran_limits(1);
        info.detection_method = "forced conservative";
        return info;
    }
    
    // HPC Detection
    info.is_hpc = (std::getenv("SLURM_JOB_ID") != nullptr) ||
        (std::getenv("PBS_JOBID") != nullptr) ||
        (std::getenv("LSB_JOBID") != nullptr) ||
        (std::getenv("SGE_JOB_ID") != nullptr) ||
        (std::getenv("SLURM_CPUS_PER_TASK") != nullptr) ||
        (std::getenv("PBS_NCPUS") != nullptr);
    
    if (info.is_hpc) {
        info.system_type = "HPC Cluster";
        info.detection_method = "HPC scheduler detection";
        
        // Get scheduler allocation
        const char* slurm_cpus = std::getenv("SLURM_CPUS_PER_TASK");
        const char* pbs_cpus = std::getenv("PBS_NCPUS");
        const char* omp_threads = std::getenv("OMP_NUM_THREADS");
        
        int requested_threads = info.num_cores; // fallback
        
        if (slurm_cpus) {
            requested_threads = std::atoi(slurm_cpus);
            info.detection_method += " (SLURM)";
        } else if (pbs_cpus) {
            requested_threads = std::atoi(pbs_cpus);
            info.detection_method += " (PBS)";
        } else if (omp_threads) {
            requested_threads = std::atoi(omp_threads);
            info.detection_method += " (OMP_NUM_THREADS)";
        }
        
        // Apply CRAN limits if necessary
        info.recommended_threads = apply_cran_limits(requested_threads);
        
        // Only conservative if CRAN-limited
        info.is_conservative_mode = (info.recommended_threads < requested_threads);
        
        return info;
    }
    
    // Linux Server Detection
    if (info.os_name == "Linux" && info.num_cores > 16) {
        info.is_linux_server = true;
        info.system_type = "Linux Server";
        
        int requested_threads = info.num_cores;
        info.recommended_threads = apply_cran_limits(requested_threads);
        info.detection_method = "Linux server heuristic";
        
        // Conservative only if CRAN-limited
        info.is_conservative_mode = (info.recommended_threads < requested_threads);
        
        return info;
    }
    
    // Standard system types
    int requested_threads;
    if (info.os_name == "macOS") {
        info.system_type = "Mac Desktop";
        requested_threads = std::max(1, info.num_cores - 1);
        info.detection_method = "macOS detection";
    } else if (info.os_name == "Windows") {
        info.system_type = "Windows Desktop";
        requested_threads = std::max(1, info.num_cores / 2);
        info.detection_method = "Windows detection";
    } else if (info.os_name == "Linux") {
        info.system_type = "Linux Desktop";
        requested_threads = std::max(1, info.num_cores / 2);
        info.detection_method = "Linux desktop detection";
    } else {
        info.system_type = "Unknown System";
        requested_threads = std::max(1, info.num_cores / 2);
        info.detection_method = "fallback";
    }
    
    // Apply CRAN limits
    info.recommended_threads = apply_cran_limits(requested_threads);
    info.is_conservative_mode = (info.recommended_threads < requested_threads);
    
    return info;
}


/**
 * @brief Gets system information (cached after first call)
 * @return const SystemInfo& Reference to system information
 */
static const SystemInfo& get_system_info() {
    if (!g_system_detected) {
        g_system_info = detect_system_type();
        g_system_detected = true;
    }
    return g_system_info;
}

/**
 * @brief Retrieves integer value from environment variable (existing function)
 */
static int getIntEnv(const char *name, int def)
{
    const char *val = getenv(name);
    if (val==NULL) return def;
    size_t nchar = strlen(val);
    if (nchar==0) return def;
    char *end;
    errno = 0;
    long int ans = strtol(val, &end, 10);
    while (isspace(*end)) end++;
    if (errno || (size_t)(end-val)!=nchar || ans<1 || ans>INT_MAX) {
        Rcpp::warning("Ignoring invalid %s==\"%s\". Not an integer >= 1. Please remove any characters that are not a digit [0-9].", name, val);
        return def;
    }
    return (int)ans;
}
    
/**
 * @brief Minimum of two integers (existing)
 */
static inline int imin(int a, int b) { return a < b ? a : b; }

/**
 * @brief Maximum of two integers (existing)
 */
static inline int imax(int a, int b) { return a > b ? a : b; }

#endif 

#ifndef _INITDTTHREADS
#define _INITDTTHREADS

/**
 * @brief Thread initialization with system-aware configuration
 * 
 */
inline void initDTthreads() {
    // Get system information
    const SystemInfo& sys_info = get_system_info();
    
    int ans = getIntEnv("R_DATATABLE_NUM_THREADS", INT_MIN);
    
    if (ans >= 1) {
        ans = imin(ans, omp_get_num_procs());
    } else {
        // Logic for different system types
        if (sys_info.is_conservative_mode) {
            // For forced conservative mode only
            int perc = getIntEnv("R_DATATABLE_NUM_PROCS_PERCENT", 25);
            perc = std::min(perc, 50);
            ans = imax(omp_get_num_procs() * perc / 100, 1);
            ans = std::min(ans, 2); // Maximum 2 threads for conservative
        } else {
            // Standard logic for all other systems (including HPC)
            int perc = getIntEnv("R_DATATABLE_NUM_PROCS_PERCENT", 50);
            if (perc <= 1 || perc > 100) {
                Rcpp::warning("Ignoring invalid R_DATATABLE_NUM_PROCS_PERCENT==%d. If used it must be an integer between 2 and 100. Default is 50. See ?setDTtheads.", perc);
                perc = 50;
            }
            ans = imax(omp_get_num_procs() * perc / 100, 1);
        }
    }
    
    // Apply standard limits
    ans = imin(ans, omp_get_thread_limit());
    ans = imin(ans, omp_get_max_threads());
    ans = imin(ans, getIntEnv("OMP_THREAD_LIMIT", INT_MAX));
    ans = imin(ans, getIntEnv("OMP_NUM_THREADS", INT_MAX));
    ans = imax(ans, 1);
    
    DTthreads = ans;
    DTthrottle = imax(1, getIntEnv("R_DATATABLE_THROTTLE", 
                                   sys_info.is_conservative_mode ? 512 : 1024));
}
#endif 

#ifndef _GETDTHREADS
#define _GETDTHREADS

/**
 * @brief Thread count determination 
 */
inline int getDTthreads(const int64_t n, const bool throttle) {
    initDTthreads();
    
    if (n < 1) return 1;
    
    int64_t ans;
    if (throttle) {
        ans = 1 + (n - 1) / DTthrottle;
    } else {
        ans = n;
    }
    
    int final_threads = (ans >= DTthreads) ? DTthreads : (int)ans;
    
    return final_threads;
}
#endif 

#ifndef _MYGETENV
#define _MYGETENV
/**
 * @brief Safe environment variable retrieval with default
 */
static const char *mygetenv(const char *name, const char *unset) {
    const char *ans = getenv(name);
    return (ans==NULL || ans[0]=='\0') ? unset : ans;
}
#endif 

#ifndef _GETDTHREADS_R
#define _GETDTHREADS_R

/**
 * @brief R interface with system information
 */
inline SEXP getDTthreads_R(SEXP verbose) {
    if(!IS_TRUE_OR_FALSE(verbose))
        Rf_error("%s must be TRUE or FALSE", "verbose");
        
    if (LOGICAL(verbose)[0]) {
        const SystemInfo& sys_info = get_system_info();
        
#ifndef _OPENMP
        Rprintf("This installation of BigDataStatMeth has not been compiled with OpenMP support.\n");
#else
        Rprintf("  OpenMP version (_OPENMP)       %d\n", _OPENMP);
#endif
        // System information
        Rprintf("  System Type                    %s\n", sys_info.system_type.c_str());
        Rprintf("  Operating System               %s\n", sys_info.os_name.c_str());
        Rprintf("  Detection Method               %s\n", sys_info.detection_method.c_str());
        Rprintf("  Conservative Mode              %s\n", sys_info.is_conservative_mode ? "YES" : "NO");
        if (sys_info.is_hpc) {
            Rprintf("  HPC Environment                DETECTED\n");
        }
        if (sys_info.is_linux_server) {
            Rprintf("  Linux Server                   DETECTED\n");
        }
        
        // Standard information
        Rprintf("  omp_get_num_procs()            %d\n", omp_get_num_procs());
        Rprintf("  Recommended Threads            %d\n", sys_info.recommended_threads);
        Rprintf("  R_DATATABLE_NUM_PROCS_PERCENT  %s\n", mygetenv("R_DATATABLE_NUM_PROCS_PERCENT", "unset (default 50)"));
        Rprintf("  R_DATATABLE_NUM_THREADS        %s\n", mygetenv("R_DATATABLE_NUM_THREADS", "unset"));
        Rprintf("  R_DATATABLE_THROTTLE           %s\n", mygetenv("R_DATATABLE_THROTTLE", "unset (default 1024)"));
        Rprintf("  omp_get_thread_limit()         %d\n", omp_get_thread_limit());
        Rprintf("  omp_get_max_threads()          %d\n", omp_get_max_threads());
        Rprintf("  OMP_THREAD_LIMIT               %s\n", mygetenv("OMP_THREAD_LIMIT", "unset"));
        Rprintf("  OMP_NUM_THREADS                %s\n", mygetenv("OMP_NUM_THREADS", "unset"));
        Rprintf("  RestoreAfterFork               %s\n", RestoreAfterFork ? "true" : "false");
        Rprintf("  BigDataStatMeth is using %d threads with throttle==%d.\n", getDTthreads(INT_MAX, false), DTthrottle);
        
        // Performance tips
        if (sys_info.is_hpc || sys_info.is_linux_server) {
            Rprintf("\n  PERFORMANCE NOTE: Smart optimization applied for %s\n", sys_info.system_type.c_str());
            Rprintf("  This should provide optimal performance while using all allocated resources.\n");
        }
    }
    return Rf_ScalarInteger(getDTthreads(INT_MAX, false));
}
#endif

// Existing fork handling code (unchanged)
#ifndef _PREFORKDTTHREADS
#define _PREFORKDTTHREADS
    static int pre_fork_DTthreads = 0;
#endif

#ifndef _WHENFORK
#define _WHENFORK
    inline void when_fork() {
        pre_fork_DTthreads = DTthreads;
        DTthreads = 1;
    }
#endif

#ifndef _AFTERFORK
#define _AFTERFORK
    inline void after_fork() {
        if (RestoreAfterFork) DTthreads = pre_fork_DTthreads;
    }
#endif

#ifndef _AVOID_OPENMP_HANG
#define _AVOID_OPENMP_HANG
    inline void avoid_openmp_hang_within_fork() {
#ifdef _OPENMP
        pthread_atfork(&when_fork, &after_fork, NULL);
#endif
    }
#endif
    
#ifndef _GET_FINAL_THREADS
#define _GET_FINAL_THREADS
    /**
     * @brief Thread determination
     */
    inline unsigned int get_number_threads(Rcpp::Nullable<int> threads, Rcpp::Nullable<bool> bparal) {
        
        unsigned int ithreads = std::thread::hardware_concurrency();
        
        if( bparal.isNotNull() ) {
            if(Rcpp::as<bool>(bparal) == false) {
                ithreads = 1;
                return(ithreads);
            }
        }
        
        if(threads.isNotNull()) {
            int requested_threads = Rcpp::as<int>(threads);
            if (requested_threads <= (int)ithreads){
                ithreads = requested_threads;
            }
        } else {
            ithreads = getDTthreads(0, false);
        }    
        
        return(ithreads);
    }
#endif

// NEW: Smart HPC optimization functions

#ifndef _SMART_HPC_OPTIMIZATION
#define _SMART_HPC_OPTIMIZATION

/**
 * @brief Force conservative mode (useful for debugging)
 */
inline void bd_force_conservative_mode(bool enable = true) {
    g_conservative_forced = enable;
    g_system_detected = false; // Force re-detection
    // Note: We rely on R's Sys.setenv() for environment variable setting
    // since setenv/unsetenv are not portable across all systems
}

/**
 * @brief Apply HPC-optimized OpenMP settings without limiting threads
 */
inline void apply_hpc_optimized_settings(int num_threads) {
#ifdef _OPENMP
    const SystemInfo& sys_info = get_system_info();
    
    // Set the requested number of threads
    omp_set_num_threads(num_threads);
    
    if (sys_info.is_hpc) {
        // HPC-specific optimizations (but use all allocated threads)
        omp_set_dynamic(false);        // Disable dynamic adjustment
        omp_set_nested(false);         // Disable nested parallelism
        omp_set_schedule(omp_sched_static, 0); // Use static scheduling
        
    } else if (sys_info.is_linux_server) {
        // Linux server optimizations
        omp_set_dynamic(false);
        omp_set_nested(false);
        omp_set_schedule(omp_sched_static, 0);
    } else {
        // Desktop optimizations
        omp_set_dynamic(true);
        omp_set_schedule(omp_sched_dynamic, 1000);
    }
#endif
}

/**
 * @brief Get optimal number of threads respecting system allocation
 *        CRAN-compatible thread optimization
 */
inline int get_optimal_threads(Rcpp::Nullable<int> user_threads) {
    
    const SystemInfo& sys_info = get_system_info();
    
    if (user_threads.isNotNull()) {
        // User specified threads - but apply CRAN limits
        int requested = Rcpp::as<int>(user_threads);
        return apply_cran_limits(requested);
    }
    
    // For HPC, try to use scheduler allocation (with CRAN limits)
    if (sys_info.is_hpc) {
        const char* slurm_cpus = std::getenv("SLURM_CPUS_PER_TASK");
        const char* pbs_cpus = std::getenv("PBS_NCPUS");
        const char* omp_threads = std::getenv("OMP_NUM_THREADS");
        
        if (slurm_cpus) {
            return apply_cran_limits(std::atoi(slurm_cpus));
        } else if (pbs_cpus) {
            return apply_cran_limits(std::atoi(pbs_cpus));
        } else if (omp_threads) {
            return apply_cran_limits(std::atoi(omp_threads));
        }
    }
    
    // Fallback to getDTthreads (CRAN-compliant)
    return getDTthreads(INT_MAX, false);
}

#endif