# R/memory_utils.R
# Utilidades para detectar y gestionar memoria del sistema
# NOTA: Usa implementación C++ (SystemInfo class) para máxima eficiencia

#' Get dynamic memory thresholds based on system RAM
#'
#' @description
#' Calculates appropriate memory thresholds for conversions based on
#' total system RAM. Returns conservative values suitable for most systems.
#'
#' @param total_ram Numeric. Total system RAM in GB. If NULL (default),
#'   auto-detects using C++ implementation.
#'
#' @return Named list with thresholds in MB:
#' \describe{
#'   \item{silent}{Size below which conversions happen silently}
#'   \item{warning}{Size requiring user confirmation}
#'   \item{force}{Size requiring force=TRUE}
#'   \item{blocked}{Size that cannot be converted}
#' }
#'
#' @details
#' **Calculation logic:**
#' 
#' Uses percentage of total RAM with safety margins:
#' \itemize{
#'   \item Silent: 15\% of RAM (up to 3GB max)
#'   \item Warning: 30\% of RAM (up to 8GB max)
#'   \item Force: 50\% of RAM (up to 16GB max)
#'   \item Blocked: 80\% of RAM (hard limit)
#' }
#'
#' **Fallback values** (if RAM detection fails):
#' \itemize{
#'   \item Silent: 2 GB
#'   \item Warning: 4 GB
#'   \item Force: 8 GB
#'   \item Blocked: 16 GB
#' }
#'
#' @examples
#' \donttest{
#' # Auto-detect
#' thresholds <- get_memory_thresholds()
#' # $silent: 3000 MB, $warning: 8000 MB, etc.
#'
#' # Manual specification (e.g., for 32GB system)
#' thresholds <- get_memory_thresholds(total_ram = 32)
#' }
#'
#' @seealso \code{\link{get_total_ram}}, \code{\link{memory_info}}
#'
#' @export
get_memory_thresholds <- function(total_ram = NULL) {
  
  # Try to detect total RAM if not provided
  # Uses C++ implementation for efficiency
  if (is.null(total_ram) || is.na(total_ram)) {
    total_ram <- get_total_ram()
  }
  
  # If detection failed, use conservative defaults
  if (is.na(total_ram) || total_ram <= 0) {
    return(list(
      silent  = 2048,   # 2 GB
      warning = 4096,   # 4 GB
      force   = 8192,   # 8 GB
      blocked = 16384   # 16 GB
    ))
  }
  
  # Calculate thresholds as percentage of total RAM
  # Use MB for consistency with existing code
  total_mb <- total_ram * 1024
  
  thresholds <- list(
    # Silent: 15% of RAM, capped at 3GB
    silent = min(total_mb * 0.15, 3072),
    
    # Warning: 30% of RAM, capped at 8GB
    warning = min(total_mb * 0.30, 8192),
    
    # Force required: 50% of RAM, capped at 16GB
    force = min(total_mb * 0.50, 16384),
    
    # Blocked: 80% of RAM (no cap, but practical limit)
    blocked = total_mb * 0.80
  )
  
  # Round to nearest 100 MB for cleaner numbers
  thresholds <- lapply(thresholds, function(x) round(x / 100) * 100)
  
  return(thresholds)
}


#' Print system memory information
#'
#' @description
#' Displays system memory information and current conversion thresholds.
#' Useful for debugging and understanding memory limits.
#'
#' @return Invisible list with memory information
#'
#' @examples
#' \donttest{
#' memory_info()
#' 
#' # Conversion Thresholds:
#' #   Silent:  2.4 GB
#' #   Warning: 4.8 GB
#' #   Force:   8.0 GB
#' #   Blocked: 12.8 GB
#' }
#'
#' @seealso \code{\link{system_info}}, \code{\link{get_total_ram}}
#'
#' @export
memory_info <- function() {
  
  # Get system info from C++ implementation
  info <- system_info()
  
  thresholds <- get_memory_thresholds(info$total_ram_gb)
  
  cat("System Information:\n")
  cat(sprintf("  OS:          %s\n", info$os))
  cat(sprintf("  CPU Cores:   %d\n", info$cpu_cores))
  
  cat("\nMemory:\n")
  if (!is.na(info$total_ram_gb) && info$total_ram_gb > 0) {
    cat(sprintf("  Total RAM:   %.1f GB\n", info$total_ram_gb))
  } else {
    cat("  Total RAM:   Unable to detect\n")
  }
  
  if (!is.na(info$available_ram_gb) && info$available_ram_gb > 0) {
    cat(sprintf("  Available:   %.1f GB", info$available_ram_gb))
    if (!is.na(info$ram_used_pct)) {
      cat(sprintf(" (%.0f%% used)", info$ram_used_pct))
    }
    cat("\n")
  } else {
    cat("  Available:   Unable to detect\n")
  }
  
  cat("\nConversion Thresholds:\n")
  cat(sprintf("  Silent:  %.1f GB\n", thresholds$silent / 1024))
  cat(sprintf("  Warning: %.1f GB\n", thresholds$warning / 1024))
  cat(sprintf("  Force:   %.1f GB\n", thresholds$force / 1024))
  cat(sprintf("  Blocked: %.1f GB\n", thresholds$blocked / 1024))
  
  invisible(list(
    os = info$os,
    cpu_cores = info$cpu_cores,
    total_gb = info$total_ram_gb,
    available_gb = info$available_ram_gb,
    ram_used_pct = info$ram_used_pct,
    thresholds_mb = thresholds
  ))
}


#' Get recommended number of threads for parallel operations
#'
#' @description
#' Returns a recommended number of threads for parallel operations,
#' based on available CPU cores and system load.
#'
#' @param use_fraction Numeric. Fraction of available cores to use (default 0.8).
#'   Using all cores (1.0) may overload the system.
#'
#' @return Integer with recommended number of threads (minimum 1)
#'
#' @details
#' This function uses the C++ implementation to detect CPU cores and
#' returns a conservative estimate leaving room for the OS and other processes.
#'
#' **Default behavior** (use_fraction = 0.8):
#' - 4 cores → 3 threads
#' - 8 cores → 6 threads
#' - 16 cores → 13 threads
#'
#' @examples
#' \donttest{
#' # Get recommended threads
#' threads <- get_recommended_threads()
#' 
#' # Use for OpenMP operations
#' options(BigDataStatMeth.threads = threads)
#' 
#' # More aggressive (use 90% of cores)
#' threads <- get_recommended_threads(use_fraction = 0.9)
#' }
#'
#' @seealso \code{\link{get_cpu_cores}}
#'
#' @export
get_recommended_threads <- function(use_fraction = 0.8) {
  
  if (use_fraction < 0.1 || use_fraction > 1.0) {
    stop("use_fraction must be between 0.1 and 1.0")
  }
  
  cores <- get_cpu_cores()
  
  if (cores <= 0) {
    # Detection failed, return conservative default
    return(1L)
  }
  
  # Calculate recommended threads
  threads <- max(1L, as.integer(floor(cores * use_fraction)))
  
  return(threads)
}
