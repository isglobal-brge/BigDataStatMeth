# Global options for HDF5Matrix operations
# Provides centralized control over parallelization, block processing and compression

#' @title HDF5Matrix Global Options
#' @description
#' Internal environment to store global computation options for HDF5Matrix operations.
#' @keywords internal
.hdf5matrix_options <- new.env(parent = emptyenv())

# Set defaults
.hdf5matrix_options$paral       <- NULL  # NULL = auto-detect
.hdf5matrix_options$block_size  <- NULL  # NULL = auto-calculate
.hdf5matrix_options$threads     <- NULL  # NULL = use all available
.hdf5matrix_options$compression <- NULL  # NULL = use default (6)


#' Set or get HDF5Matrix computation options
#'
#' @description
#' Configure global settings for parallelization, block processing and compression in
#' HDF5Matrix operations. These settings affect all HDF5Matrix computations
#' unless explicitly overridden in individual method calls.
#'
#' @param paral Logical or NULL. Enable OpenMP parallelization?
#'   \itemize{
#'     \item \code{TRUE}: Force parallel execution
#'     \item \code{FALSE}: Force serial execution
#'     \item \code{NULL}: Let BigDataStatMeth auto-detect (default)
#'   }
#' @param block_size Integer or NULL. Number of elements per block for
#'   block-wise processing.
#'   \itemize{
#'     \item Integer > 0: Use this block size
#'     \item \code{NULL}: Auto-calculate based on matrix dimensions (default)
#'   }
#' @param threads Integer or NULL. Number of OpenMP threads to use.
#'   \itemize{
#'     \item Integer > 0: Use this many threads
#'     \item \code{NULL}: Use all available threads (default)
#'   }
#' @param compression Integer (0-9) or NULL. gzip compression level for
#'   created datasets.
#'   \itemize{
#'     \item \code{0}: No compression (fastest, largest files)
#'     \item \code{1-3}: Light compression (fast, moderate savings)
#'     \item \code{6}: Balanced compression (default, 60-80\% space savings)
#'     \item \code{7-9}: Maximum compression (slowest, best ratio)
#'     \item \code{NULL}: Use built-in default of 6
#'   }
#'
#' @return
#' When called with arguments: invisibly returns a list of all current options.
#' When called without arguments: returns a list of all current options.
#'
#' @details
#' BigDataStatMeth achieves high performance through two key mechanisms:
#'
#' **Block-wise processing:**
#' Large matrices are processed in chunks that fit in memory. The \code{block_size}
#' parameter controls chunk size. Smaller blocks use less memory but require more
#' I/O operations. Larger blocks are faster but require more RAM.
#'
#' **OpenMP parallelization:**
#' Operations are distributed across CPU cores. The \code{paral} and \code{threads}
#' parameters control this. Parallelization provides near-linear speedup for
#' compute-intensive operations.
#'
#' **Compression:**
#' Datasets are created with gzip compression (level 6 by default). This reduces
#' disk usage by 60-80\% at the cost of additional CPU time for compress/decompress.
#' For benchmarks or workflows where speed is critical, set \code{compression = 0}.
#' For long-term storage or large datasets, keep the default.
#'
#' **Priority:**
#' Options set here serve as defaults. Individual method calls can override:
#' \code{A$multiply(B, paral = TRUE, threads = 4, block_size = 2000)}
#'
#' **Recommendations:**
#' \itemize{
#'   \item For interactive analysis: Leave defaults (\code{NULL}) - auto-detect works well
#'   \item For scripts/HPC: Set explicitly based on your hardware and data size
#'   \item For huge datasets (>10GB): Reduce \code{block_size} to fit in RAM
#'   \item For many-core systems: Set \code{threads} explicitly (auto may be too aggressive)
#'   \item For benchmarks: Set \code{compression = 0} to eliminate gzip overhead
#' }
#'
#' @examples
#' # View current options
#' hdf5matrix_options()
#'
#' # Enable parallelization with 8 threads
#' hdf5matrix_options(paral = TRUE, threads = 8)
#'
#' # Set block size to 1000 elements
#' hdf5matrix_options(block_size = 1000)
#'
#' # Disable compression for benchmarking
#' hdf5matrix_options(compression = 0)
#'
#' # Reset to defaults
#' hdf5matrix_options(paral = NULL, threads = NULL, block_size = NULL, compression = NULL)
#'
#' # Now all operations use these settings
#' C <- A %*% B              # Uses options set above
#' D <- crossprod(A)         # Uses options set above
#'
#' # Override for specific operation
#' E <- A$multiply(B, paral = FALSE)  # Forces serial, ignores global option
#'
#' @export
hdf5matrix_options <- function(paral = NULL, block_size = NULL, threads = NULL,
                                compression = NULL) {
  # If called with no arguments, just return current options
  if (missing(paral) && missing(block_size) && missing(threads) &&
      missing(compression)) {
    return(list(
      paral       = .hdf5matrix_options$paral,
      block_size  = .hdf5matrix_options$block_size,
      threads     = .hdf5matrix_options$threads,
      compression = .hdf5matrix_options$compression
    ))
  }

  # Validate inputs
  if (!is.null(paral) && !is.logical(paral)) {
    stop("paral must be TRUE, FALSE, or NULL")
  }

  if (!is.null(block_size)) {
    if (!is.numeric(block_size) || block_size <= 0) {
      stop("block_size must be a positive integer or NULL")
    }
    block_size <- as.integer(block_size)
  }

  if (!is.null(threads)) {
    if (!is.numeric(threads) || threads <= 0) {
      stop("threads must be a positive integer or NULL")
    }
    threads <- as.integer(threads)
  }

  if (!is.null(compression)) {
    if (!is.numeric(compression) || compression < 0 || compression > 9) {
      stop("compression must be an integer 0-9 or NULL")
    }
    compression <- as.integer(compression)
  }

  # Set options
  if (!missing(paral))       .hdf5matrix_options$paral       <- paral
  if (!missing(block_size))  .hdf5matrix_options$block_size  <- block_size
  if (!missing(threads))     .hdf5matrix_options$threads     <- threads
  if (!missing(compression)) .hdf5matrix_options$compression <- compression

  # Return current state invisibly
  invisible(list(
    paral       = .hdf5matrix_options$paral,
    block_size  = .hdf5matrix_options$block_size,
    threads     = .hdf5matrix_options$threads,
    compression = .hdf5matrix_options$compression
  ))
}


#' Get global option value with fallback
#'
#' @description
#' Internal helper to retrieve global option value. If the option is NULL
#' (not set), returns the provided default. If a non-NULL value is explicitly
#' passed (from a method call), that takes priority over the global option.
#'
#' @param name Option name ("paral", "block_size", "threads", or "compression")
#' @param default Fallback value if option is NULL
#' @param override Value passed to method call (takes priority if not NULL/missing)
#'
#' @return The effective value to use
#'
#' @keywords internal
.get_option <- function(name, default = NULL, override = NULL) {
  # Priority: override > global option > default
  if (!is.null(override)) {
    return(override)
  }

  global_val <- .hdf5matrix_options[[name]]
  if (!is.null(global_val)) {
    return(global_val)
  }

  return(default)
}


#' Show current HDF5Matrix performance settings
#'
#' @description
#' Display current global options in a user-friendly format.
#'
#' @return Invisibly returns the options list
#'
#' @examples
#' show_hdf5matrix_options()
#'
#' @export
show_hdf5matrix_options <- function() {
  opts <- hdf5matrix_options()

  cat("\nHDF5Matrix Global Options\n")
  cat("=========================\n\n")

  cat("Parallelization:\n")
  cat("  paral:      ", if (is.null(opts$paral)) "NULL (auto-detect)"
                        else if (opts$paral) "TRUE (enabled)"
                        else "FALSE (disabled)", "\n")
  cat("  threads:    ", if (is.null(opts$threads)) "NULL (use all available)"
                        else paste(opts$threads, "threads"), "\n\n")

  cat("Block Processing:\n")
  cat("  block_size: ", if (is.null(opts$block_size)) "NULL (auto-calculate)"
                        else paste(opts$block_size, "elements"), "\n\n")

  cat("Compression:\n")
  cat("  compression:", if (is.null(opts$compression)) "NULL (default = 6, gzip balanced)"
                        else if (opts$compression == 0) "0 (disabled)"
                        else paste0(opts$compression, " (gzip level)"), "\n\n")

  cat("Note: These are defaults. Individual methods can override.\n")
  cat("Example: A$multiply(B, paral = TRUE, threads = 4, compression = 0)\n\n")

  invisible(opts)
}
