/**
 * @file SystemInfo.hpp
 * @brief System information utilities for memory and CPU detection
 * @details Provides cross-platform system information including RAM and CPU cores.
 * 
 * This header provides utilities for querying system resources, useful for:
 * - Memory management decisions (can we allocate this matrix?)
 * - Parallel processing configuration (how many threads?)
 * - User warnings about resource constraints
 * 
 * Supported platforms:
 * - Windows (Vista+)
 * - Linux (with /proc filesystem)
 * - macOS (10.4+)
 * 
 * All functions are inline and header-only for easy integration.
 * 
 * @author BigDataStatMeth team
 * @date 2026-02-20
 */

#ifndef BIGDATASTATMETH_SYSTEMINFO_HPP
#define BIGDATASTATMETH_SYSTEMINFO_HPP

#include <cstddef>
#include <string>
#include <fstream>
#include <sstream>

// Platform-specific headers
#ifdef _WIN32
    #include <windows.h>
#elif defined(__linux__)
    #include <unistd.h>
    #include <fstream>
#elif defined(__APPLE__)
    #include <unistd.h>
    #include <sys/types.h>
    #include <sys/sysctl.h>
    #include <mach/mach.h>
    #include <mach/mach_host.h>
#endif

namespace BigDataStatMeth {

/**
 * @class SystemInfo
 * @brief System resource information utilities
 * 
 * @details
 * Provides static methods to query system resources in a cross-platform manner.
 * All values are returned in megabytes (MB) for consistency.
 * 
 * **Thread safety:** All methods are thread-safe as they only read system state.
 * 
 * **Error handling:** If detection fails, methods return 0. Check return value
 * before using.
 * 
 * **Example usage:**
 * @code{.cpp}
 * // Check if we can allocate 5GB
 * size_t needed_mb = 5000;
 * if (SystemInfo::canAllocate(needed_mb)) {
 *     // Safe to proceed
 *     Eigen::MatrixXd mat(nrows, ncols);
 * } else {
 *     Rf_error("Insufficient memory: need %zu MB", needed_mb);
 * }
 * 
 * // Get optimal thread count
 * int threads = SystemInfo::getCPUCores();
 * omp_set_num_threads(threads);
 * @endcode
 */
class SystemInfo {
public:
    
    /**
     * @brief Get total system RAM in megabytes
     * 
     * @return Total RAM in MB, or 0 if detection fails
     * 
     * @details
     * Returns the total physical RAM installed in the system.
     * This value is constant for a given system.
     * 
     * **Platform implementation:**
     * - Windows: Uses GlobalMemoryStatusEx
     * - Linux: Reads MemTotal from /proc/meminfo
     * - macOS: Uses sysctl hw.memsize
     * 
     * @note On virtualized systems, returns the RAM allocated to the VM,
     * not the host's total RAM.
     */
    static inline size_t getTotalRAM_MB() {
        
#ifdef _WIN32
        // Windows: Use GlobalMemoryStatusEx
        MEMORYSTATUSEX statex;
        statex.dwLength = sizeof(statex);
        if (GlobalMemoryStatusEx(&statex)) {
            return static_cast<size_t>(statex.ullTotalPhys / (1024 * 1024));
        }
        return 0;
        
#elif defined(__linux__)
        // Linux: Read /proc/meminfo
        std::ifstream meminfo("/proc/meminfo");
        if (!meminfo.is_open()) {
            return 0;
        }
        
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.compare(0, 9, "MemTotal:") == 0) {
                std::istringstream iss(line.substr(9));
                size_t kb;
                if (iss >> kb) {
                    return kb / 1024;  // Convert KB to MB
                }
            }
        }
        return 0;
        
#elif defined(__APPLE__)
        // macOS: Use sysctl
        int64_t mem_bytes;
        size_t len = sizeof(mem_bytes);
        if (sysctlbyname("hw.memsize", &mem_bytes, &len, NULL, 0) == 0) {
            return static_cast<size_t>(mem_bytes / (1024 * 1024));
        }
        return 0;
        
#else
        // Unknown platform
        return 0;
#endif
    }
    
    
    /**
     * @brief Get available (free) system RAM in megabytes
     * 
     * @return Available RAM in MB, or 0 if detection fails
     * 
     * @details
     * Returns the amount of RAM that can be allocated without swapping.
     * This value changes dynamically as processes allocate/free memory.
     * 
     * **Platform implementation:**
     * - Windows: Uses GlobalMemoryStatusEx (ullAvailPhys)
     * - Linux: Reads MemAvailable from /proc/meminfo
     * - macOS: Uses vm_statistics64 (free + inactive pages)
     * 
     * @note
     * - This is an estimate; actual allocatable memory may be less
     * - On Linux, uses MemAvailable (kernel 3.14+) which is more accurate
     *   than just MemFree
     * - Value can change rapidly; use immediately after calling
     * 
     * @warning Do not cache this value. Query it each time you need it.
     */
    static inline size_t getAvailableRAM_MB() {
        
#ifdef _WIN32
        // Windows: Use GlobalMemoryStatusEx
        MEMORYSTATUSEX statex;
        statex.dwLength = sizeof(statex);
        if (GlobalMemoryStatusEx(&statex)) {
            return static_cast<size_t>(statex.ullAvailPhys / (1024 * 1024));
        }
        return 0;
        
#elif defined(__linux__)
        // Linux: Read /proc/meminfo
        // Prefer MemAvailable (kernel 3.14+) over MemFree
        std::ifstream meminfo("/proc/meminfo");
        if (!meminfo.is_open()) {
            return 0;
        }
        
        size_t mem_available_kb = 0;
        size_t mem_free_kb = 0;
        bool found_available = false;
        
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.compare(0, 13, "MemAvailable:") == 0) {
                std::istringstream iss(line.substr(13));
                if (iss >> mem_available_kb) {
                    found_available = true;
                    break;
                }
            } else if (line.compare(0, 8, "MemFree:") == 0) {
                std::istringstream iss(line.substr(8));
                iss >> mem_free_kb;
            }
        }
        
        // Prefer MemAvailable, fallback to MemFree
        size_t kb = found_available ? mem_available_kb : mem_free_kb;
        return kb / 1024;  // Convert KB to MB
        
#elif defined(__APPLE__)
        // macOS: Use vm_statistics64
        vm_statistics64_data_t vm_stat;
        mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
        
        if (host_statistics64(mach_host_self(), HOST_VM_INFO64,
                             (host_info64_t)&vm_stat, &count) == KERN_SUCCESS) {
            // Free + inactive pages (both can be allocated)
            // natural_t page_size;
            // host_page_size(mach_host_self(), &page_size);
            // uint64_t free_bytes = (vm_stat.free_count + vm_stat.inactive_count) * page_size;
            vm_size_t page_size;
            host_page_size(mach_host_self(), &page_size);
            uint64_t free_bytes = (vm_stat.free_count + vm_stat.inactive_count) * 
                static_cast<uint64_t>(page_size);
            
            return static_cast<size_t>(free_bytes / (1024 * 1024));
        }
        return 0;
        
#else
        // Unknown platform
        return 0;
#endif
    }
    
    
    /**
     * @brief Check if a given amount of memory can be safely allocated
     * 
     * @param size_mb Size in megabytes to check
     * @param safety_margin_pct Percentage of available RAM to keep free (default 20\%)
     * @return true if allocation is likely safe, false otherwise
     * 
     * @details
     * Checks if the requested memory can be allocated while maintaining a safety
     * margin of free RAM. This helps prevent system instability from exhausting
     * all available memory.
     * 
     * **Formula:**
     * @code
     * can_allocate = (size_mb < available_ram * (1 - safety_margin))
     * @endcode
     * 
     * **Safety margin rationale:**
     * - Default 20%: Leaves room for OS and other processes
     * - Can be adjusted based on use case
     * - Set to 0 for aggressive allocation (not recommended)
     * 
     * @param safety_margin_pct Safety margin as percentage (0-100)
     * 
     * @note
     * - This is a heuristic check, not a guarantee
     * - Allocation can still fail due to fragmentation
     * - Always handle allocation failures with try-catch
     * 
     * @warning Setting safety_margin_pct = 0 may cause system instability
     * 
     * @example
     * @code{.cpp}
     * // Conservative check (20\% margin)
     * if (SystemInfo::canAllocate(5000)) {
     *     // Proceed with 5GB allocation
     * }
     * 
     * // Aggressive check (5\% margin)
     * if (SystemInfo::canAllocate(5000, 5)) {
     *     // More likely to succeed but riskier
     * }
     * @endcode
     */
    static inline bool canAllocate(size_t size_mb, double safety_margin_pct = 20.0) {
        
        size_t available_mb = getAvailableRAM_MB();
        
        // If detection failed, be conservative
        if (available_mb == 0) {
            return false;
        }
        
        // Clamp safety margin to [0, 100]
        if (safety_margin_pct < 0.0) safety_margin_pct = 0.0;
        if (safety_margin_pct > 100.0) safety_margin_pct = 100.0;
        
        // Calculate usable RAM (available - safety margin)
        double usable_mb = available_mb * (1.0 - safety_margin_pct / 100.0);
        
        return size_mb <= static_cast<size_t>(usable_mb);
    }
    
    
    /**
     * @brief Get number of CPU cores (logical processors)
     * 
     * @return Number of CPU cores, or 0 if detection fails
     * 
     * @details
     * Returns the number of logical processors (includes hyperthreading).
     * Useful for configuring parallel processing (e.g., OpenMP threads).
     * 
     * **Platform implementation:**
     * - Windows: Uses GetSystemInfo
     * - Linux: Reads from sysconf(_SC_NPROCESSORS_ONLN)
     * - macOS: Uses sysctl hw.ncpu
     * 
     * @note
     * - Returns logical cores (with hyperthreading), not physical cores
     * - On systems with CPU pinning/affinity, may return fewer cores
     * - Use this value as a maximum; actual optimal thread count may be less
     * 
     * @warning Don't blindly use all cores; leave room for OS and other processes
     * 
     * @example
     * @code{.cpp}
     * #ifdef _OPENMP
     * int threads = SystemInfo::getCPUCores();
     * // Use 80% of cores to avoid oversubscription
     * threads = std::max(1, static_cast<int>(threads * 0.8));
     * omp_set_num_threads(threads);
     * #endif
     * @endcode
     */
    static inline size_t getCPUCores() {
        
#ifdef _WIN32
        // Windows: Use GetSystemInfo
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        return static_cast<size_t>(sysinfo.dwNumberOfProcessors);
        
#elif defined(__linux__) || defined(__APPLE__)
        // Linux/macOS: Use sysconf
        long cores = sysconf(_SC_NPROCESSORS_ONLN);
        if (cores > 0) {
            return static_cast<size_t>(cores);
        }
        return 0;
        
#else
        // Unknown platform
        return 0;
#endif
    }
    
    
    /**
     * @brief Get OS name as string
     * 
     * @return OS name ("Windows", "Linux", "macOS", or "Unknown")
     * 
     * @details
     * Simple OS identification for logging/debugging purposes.
     * 
     * @note This is compile-time detection, not runtime
     */
    static inline std::string getOSName() {
#ifdef _WIN32
        return "Windows";
#elif defined(__linux__)
        return "Linux";
#elif defined(__APPLE__)
        return "macOS";
#else
        return "Unknown";
#endif
    }
    
}; // class SystemInfo

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_SYSTEMINFO_HPP
