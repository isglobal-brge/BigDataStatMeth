/**
 * @file hdf5_systemInfo.cpp
 * @brief R wrappers for SystemInfo class
 * @details Exposes system information utilities to R users
 */

#include <BigDataStatMeth.hpp>
// #include <SystemInfo.hpp>

//' Get total system RAM
//'
//' @description
//' Returns the total physical RAM installed in the system.
//'
//' @return Numeric value with total RAM in gigabytes (GB)
//'
//' @details
//' This function queries the operating system to determine total RAM.
//' Works on Windows, Linux, and macOS.
//'
//' The value returned is the physical RAM available to the system:
//' - On physical machines: actual installed RAM
//' - On virtual machines: RAM allocated to the VM
//' - On containers: RAM limit set for the container
//'
//' @examples
//' \donttest{
//' # Check total RAM
//' total <- get_total_ram()
//' cat("System has", total, "GB of RAM\n")
//' 
//' # Returns 16.0 on a 16GB system
//' }
//'
//' @seealso \code{\link{get_available_ram}}, \code{\link{get_cpu_cores}}
//'
//' @export
// [[Rcpp::export]]
double get_total_ram() {
    size_t mb = BigDataStatMeth::SystemInfo::getTotalRAM_MB();
    return static_cast<double>(mb) / 1024.0;  // Convert MB to GB
}


//' Get available (free) system RAM
//'
//' @description
//' Returns the amount of RAM currently available for allocation.
//'
//' @return Numeric value with available RAM in gigabytes (GB)
//'
//' @details
//' This function returns the RAM that can be allocated without swapping.
//' The value changes dynamically as processes allocate and free memory.
//'
//' **Important notes:**
//' - Value can change rapidly; don't cache it
//' - On Linux, uses MemAvailable (more accurate than MemFree)
//' - Includes memory that can be reclaimed from caches
//' - Actual allocatable memory may be slightly less
//'
//' **Use case:**
//' Check available RAM before loading large datasets into memory.
//'
//' @examples
//' \donttest{
//' available <- get_available_ram()
//' cat("Available RAM:", round(available, 2), "GB\n")
//' 
//' # Use it to decide how much data to load
//' fn <- tempfile(fileext = ".h5")
//' X  <- hdf5_create_matrix(fn, "data/M",
//'                           data = matrix(rnorm(1000), 100, 10))
//' 
//' size_gb <- prod(dim(X)) * 8 / 1e9
//' if (get_available_ram() > size_gb * 1.2) {
//'   mat <- as.matrix(X)
//' } else {
//'   mat <- X[1:50, ]
//' }
//' 
//' hdf5_close_all()
//' unlink(fn)
//' }
//'
//' @seealso \code{\link{get_total_ram}}, \code{\link{can_allocate}}
//'
//' @export
// [[Rcpp::export]]
double get_available_ram() {
    size_t mb = BigDataStatMeth::SystemInfo::getAvailableRAM_MB();
    return static_cast<double>(mb) / 1024.0;  // Convert MB to GB
}


//' Check if memory allocation is safe
//'
//' @description
//' Checks whether a given amount of memory can be safely allocated
//' while maintaining a safety margin.
//'
//' @param size_gb Size in gigabytes (GB) to check
//' @param safety_margin_pct Percentage of available RAM to keep free
//'   (default 20 percent)
//'
//' @return Logical. TRUE if allocation is likely safe, FALSE otherwise
//'
//' @details
//' This function checks if the requested memory can be allocated while
//' keeping a safety margin of free RAM. This helps prevent:
//' - System instability from memory exhaustion
//' - Swapping (which degrades performance)
//' - Out-of-memory errors from other processes
//'
//' **Formula:**
//' \code{can_allocate = (size_gb < available_ram * (1 - safety_margin / 100))}
//'
//' **Safety margin guidelines:**
//' \itemize{
//'   \item 20 percent (default): Conservative, recommended for most cases
//'   \item 10 percent: Moderate, for controlled environments
//'   \item 5 percent: Aggressive, only if you know what you're doing
//'   \item 0 percent: Maximum risk, not recommended
//' }
//'
//' @note
//' This is a heuristic check, not a guarantee. Allocation can still fail
//' due to memory fragmentation or competing processes.
//'
//' @examples
//' \donttest{
//' # Check if 1 GB can be safely allocated
//' if (can_allocate(1)) {
//'   message("1 GB allocation is safe")
//' } else {
//'   message("Not enough RAM for 1 GB allocation")
//' }
//' 
//' # Use it to decide how much data to load
//' fn <- tempfile(fileext = ".h5")
//' X  <- hdf5_create_matrix(fn, "data/M",
//'                           data = matrix(rnorm(1000), 100, 10))
//' 
//' size_gb <- prod(dim(X)) * 8 / 1e9   # estimate in GB
//' if (can_allocate(size_gb)) {
//'   mat <- as.matrix(X)
//' } else {
//'   mat <- X[1:50, ]   # load subset
//' }
//' 
//' hdf5_close_all()
//' unlink(fn)
//' }
//'
//' @seealso \code{\link{get_available_ram}}
//'
//' @export
// [[Rcpp::export]]
bool can_allocate(double size_gb, double safety_margin_pct = 20.0) {
    size_t size_mb = static_cast<size_t>(size_gb * 1024.0);
    return BigDataStatMeth::SystemInfo::canAllocate(size_mb, safety_margin_pct);
}


//' Get number of CPU cores
//'
//' @description
//' Returns the number of logical CPU cores (processors) available.
//'
//' @return Integer with number of CPU cores
//'
//' @details
//' This function returns the number of logical processors, which includes
//' cores from hyperthreading/SMT. Useful for configuring parallel processing.
//'
//' **Typical values:**
//' - 4-core CPU without hyperthreading: 4
//' - 4-core CPU with hyperthreading: 8
//' - 8-core CPU with hyperthreading: 16
//'
//' **Usage for parallelization:**
//' Don't blindly use all cores. A common practice is to use 80-90 percent of
//' available cores to leave room for the OS and other processes.
//'
//' @note
//' - Returns logical cores (with hyperthreading), not physical cores
//' - On systems with CPU pinning, may return fewer cores
//' - Value reflects cores available to the process
//'
//' @examples
//' \donttest{
//' # Get CPU cores
//' cores <- get_cpu_cores()
//' cat("System has", cores, "CPU cores\n")
//'
//' # Configure parallel processing (use 80 percent of cores)
//' threads <- max(1, floor(cores * 0.8))
//' options(BigDataStatMeth.threads = threads)
//' }
//'
//' @seealso \code{\link{get_total_ram}}
//'
//' @export
// [[Rcpp::export]]
int get_cpu_cores() {
    size_t cores = BigDataStatMeth::SystemInfo::getCPUCores();
    return static_cast<int>(cores);
}


//' Get system information summary
//'
//' @description
//' Returns a comprehensive summary of system resources.
//'
//' @return Named list with system information:
//' \describe{
//'   \item{os}{Operating system name}
//'   \item{total_ram_gb}{Total RAM in GB}
//'   \item{available_ram_gb}{Available RAM in GB}
//'   \item{ram_used_pct}{Percentage of RAM currently used}
//'   \item{cpu_cores}{Number of CPU cores}
//' }
//'
//' @details
//' Convenience function that calls all system info methods and
//' returns a summary. Useful for debugging and logging.
//'
//' @examples
//' \donttest{
//' # Get full system info
//' info <- system_info()
//' print(info)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List system_info() {
    
    double total_gb = get_total_ram();
    double available_gb = get_available_ram();
    int cores = get_cpu_cores();
    
    // Calculate RAM usage percentage
    double used_pct = 0.0;
    if (total_gb > 0) {
        used_pct = ((total_gb - available_gb) / total_gb) * 100.0;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("os") = BigDataStatMeth::SystemInfo::getOSName(),
        Rcpp::Named("total_ram_gb") = total_gb,
        Rcpp::Named("available_ram_gb") = available_gb,
        Rcpp::Named("ram_used_pct") = used_pct,
        Rcpp::Named("cpu_cores") = cores
    );
}
