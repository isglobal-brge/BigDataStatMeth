# R/S3_conversion.R
# Conversión de HDF5Matrix a objetos en memoria

#' Convert HDF5Matrix to in-memory matrix
#'
#' @description
#' Reads entire HDF5 dataset into memory as a standard R matrix.
#' **WARNING**: This loads all data into RAM. For large datasets (>1GB),
#' this may cause memory exhaustion.
#'
#' @param x An \code{HDF5Matrix} object
#' @param force Logical. If \code{TRUE}, skip size warnings. Default \code{FALSE}.
#' @param max_size_mb Numeric. Maximum size in MB to convert without warning.
#'   Default is \code{NULL} (auto-detect based on system RAM). Use \code{Inf} to disable.
#' @param ... Additional arguments (currently unused)
#'
#' @return Standard R matrix with data from HDF5 file
#'
#' @details
#' **Size thresholds and behavior:**
#' 
#' \describe{
#'   \item{Small datasets (< max_size_mb):}{Convert silently}
#'   \item{Medium datasets (max_size_mb to 2GB):}{Show warning, require confirmation}
#'   \item{Large datasets (> 2GB):}{Show error, require force=TRUE}
#'   \item{Huge datasets (> 8GB):}{Block conversion even with force=TRUE}
#' }
#'
#' **Memory estimation:**
#' The function estimates memory usage as:
#' \code{nrow * ncol * 8 bytes (for numeric)}
#' 
#' Actual memory usage may be higher due to:
#' - R's internal overhead
#' - Temporary copies during conversion
#' - Other objects in memory
#'
#' **Recommendations:**
#' \itemize{
#'   \item For large datasets, use subsetting instead: \code{X[1:1000, ]}
#'   \item For analysis, use HDF5Matrix methods directly (they work on-disk)
#'   \item Only convert to memory when absolutely necessary
#' }
#'
#' @examples
#' \donttest{
#' # Small dataset - converts silently
#' X <- hdf5_create_matrix("test.h5", "data/X", nrow = 100, ncol = 50)
#' mat <- as.matrix(X)  # Works fine
#'
#' # Large dataset - shows warning
#' Y <- hdf5_create_matrix("test.h5", "data/Y", nrow = 10000, ncol = 5000)
#' mat <- as.matrix(Y)  # Warning: 381 MB - proceed? [y/N]
#'
#' # Very large - requires force=TRUE
#' Z <- hdf5_create_matrix("test.h5", "data/Z", nrow = 50000, ncol = 10000)
#' mat <- as.matrix(Z, force = TRUE)  # Explicitly allow
#'
#' # Better: Use subsetting for large datasets
#' subset <- Z[1:1000, 1:100]  # Only 0.8 MB
#' }
#'
#' @seealso 
#' \code{\link{[.HDF5Matrix}} for subsetting,
#' \code{\link{as.data.frame.HDF5Matrix}} for data frame conversion
#'
#' @export
as.matrix.HDF5Matrix <- function(x, force = FALSE, max_size_mb = NULL, ...) {
  
  if (!x$is_valid()) {
    stop("HDF5Matrix is closed or invalid")
  }
  
  # Get dimensions
  dims <- dim(x)
  nrows <- dims[1]
  ncols <- dims[2]
  
  # Estimate memory size (assuming numeric/double: 8 bytes per element)
  size_bytes <- nrows * ncols * 8
  size_mb <- size_bytes / (1024^2)
  size_gb <- size_mb / 1024
  
  # Get dynamic thresholds based on system RAM
  # User can override with max_size_mb parameter
  if (is.null(max_size_mb)) {
    thresholds <- get_memory_thresholds()
    SILENT_THRESHOLD  <- thresholds$silent
    MEDIUM_THRESHOLD  <- thresholds$warning
    LARGE_THRESHOLD   <- thresholds$force
    HUGE_THRESHOLD    <- thresholds$blocked
  } else {
    # User-specified threshold
    SILENT_THRESHOLD  <- max_size_mb
    MEDIUM_THRESHOLD  <- max_size_mb * 2
    LARGE_THRESHOLD   <- max_size_mb * 4
    HUGE_THRESHOLD    <- max_size_mb * 8
  }
  
  # Check using C++ implementation for accuracy
  available_gb <- get_available_ram()
  
  # Check size and apply warnings
  if (size_mb > HUGE_THRESHOLD) {
    # Absolutely too large
    stop(
      "Dataset is too large for memory conversion: ", 
      sprintf("%.1f GB", size_gb), "\n",
      "  Dimensions: ", nrows, " x ", ncols, "\n",
      sprintf("  System limit: %.1f GB (80%% of %.1f GB total RAM)\n", 
              HUGE_THRESHOLD / 1024, get_total_ram()),
      sprintf("  Available RAM: %.1f GB\n", available_gb),
      "  This would likely crash R.\n",
      "  Solutions:\n",
      "    1. Use subsetting: X[1:1000, ] to get a smaller portion\n",
      "    2. Work with HDF5Matrix directly (no conversion needed)\n",
      "    3. Use block-wise processing\n",
      call. = FALSE
    )
  }
  
  if (size_mb > LARGE_THRESHOLD && !force) {
    # Very large - requires explicit force
    # Use C++ canAllocate for accurate check
    can_do_it <- can_allocate(size_gb, safety_margin_pct = 20)
    
    stop(
      "Dataset is very large: ", sprintf("%.1f GB", size_gb), "\n",
      "  Dimensions: ", nrows, " x ", ncols, "\n",
      sprintf("  Threshold: %.1f GB\n", LARGE_THRESHOLD / 1024),
      sprintf("  Available RAM: %.1f GB\n", available_gb),
      if (!can_do_it) "  Insufficient RAM even with 20% safety margin\n" else "",
      "  Loading this into memory may cause problems.\n",
      "  If you're certain you have enough RAM, use:\n",
      "    as.matrix(X, force = TRUE)\n",
      "  Or better, use subsetting:\n",
      "    X[1:1000, ]  # Get first 1000 rows instead\n",
      call. = FALSE
    )
  }
  
  if (size_mb > MEDIUM_THRESHOLD && !force) {
    # Medium-large - interactive confirmation
    message(
      sprintf("Dataset size: %.1f MB (%.2f GB)", size_mb, size_gb),
      "\nDimensions: ", nrows, " x ", ncols,
      sprintf("\nAvailable RAM: %.1f GB", available_gb),
      "\nThis will load all data into memory."
    )
    
    if (interactive()) {
      response <- readline(prompt = "Proceed? [y/N]: ")
      if (!tolower(response) %in% c("y", "yes")) {
        message("Conversion cancelled.")
        return(invisible(NULL))
      }
    } else {
      # Non-interactive: require force=TRUE
      stop(
        "In non-interactive mode, use force=TRUE to convert large datasets:\n",
        "  as.matrix(X, force = TRUE)",
        call. = FALSE
      )
    }
  }
  
  # Final safety check using C++ implementation
  if (!can_allocate(size_gb, safety_margin_pct = 10)) {
    warning(
      "Low memory warning: Only ", sprintf("%.1f GB", available_gb), 
      " RAM available for ", sprintf("%.1f GB", size_gb), " dataset.\n",
      "Allocation may fail or cause swapping (slow performance)."
    )
  }
  
  # If we get here, proceed with conversion
  if (size_mb > SILENT_THRESHOLD) {
    message("Loading ", sprintf("%.1f MB", size_mb), " into memory...")
  }
  
  # Read entire dataset
  # Use subsetting to read all: X[, ] triggers rcpp_hdf5dataset_subset()
  result <- x[, , drop = FALSE]

  # Propagate dimnames (already handled in [.HDF5Matrix, but kept explicit
  # here for clarity and in case as.matrix is called via another path)
  dn <- dimnames(x)
  if (!is.null(dn)) dimnames(result) <- dn

  if (size_mb > SILENT_THRESHOLD) {
    message("Done. Matrix loaded successfully.")
  }

  return(result)
}


#' Convert HDF5Matrix to data.frame
#'
#' @description
#' Reads entire HDF5 dataset into memory as a data.frame.
#' **WARNING**: This loads all data into RAM.
#'
#' @param x An \code{HDF5Matrix} object
#' @param force Logical. If \code{TRUE}, skip size warnings.
#' @param max_size_mb Numeric. Maximum size in MB to convert without warning.
#' @param row.names Logical or character vector. Row names to use.
#' @param optional Logical. Passed to \code{as.data.frame}.
#' @param ... Additional arguments passed to \code{as.data.frame}
#'
#' @return data.frame with data from HDF5 file
#'
#' @details
#' First converts to matrix using \code{as.matrix.HDF5Matrix} (with same
#' size checks), then to data.frame. All memory warnings apply.
#'
#' @examples
#' \donttest{
#' X <- hdf5_create_matrix("test.h5", "data/X", nrow = 100, ncol = 5)
#' df <- as.data.frame(X)
#' }
#'
#' @seealso \code{\link{as.matrix.HDF5Matrix}}
#'
#' @export
as.data.frame.HDF5Matrix <- function(x, 
                                     force = FALSE, 
                                     max_size_mb = NULL,
                                     row.names = NULL, 
                                     optional = FALSE, 
                                     ...) {
  
  # Convert to matrix first (includes size checks)
  mat <- as.matrix(x, force = force, max_size_mb = max_size_mb)
  
  if (is.null(mat)) {
    return(invisible(NULL))
  }
  
  # Convert matrix to data.frame
  as.data.frame(mat, row.names = row.names, optional = optional, ...)
}


#' Get memory size of HDF5Matrix without loading
#'
#' @description
#' Estimates how much memory the dataset would occupy if loaded into RAM.
#'
#' @param x An \code{HDF5Matrix} object
#' @param unit Character. Unit for size: "bytes", "KB", "MB", "GB". Default "MB".
#'
#' @return Numeric value with estimated memory size
#'
#' @examples
#' \donttest{
#' X <- hdf5_create_matrix("test.h5", "data/X", nrow = 10000, ncol = 5000)
#' object_size(X)              # 381.5 MB
#' object_size(X, unit = "GB") # 0.37 GB
#' }
#'
#' @export
object_size <- function(x, unit = c("MB", "bytes", "KB", "GB")) {
  UseMethod("object_size")
}

#' @export
object_size.HDF5Matrix <- function(x, unit = c("MB", "bytes", "KB", "GB")) {
  unit <- match.arg(unit)
  
  if (!x$is_valid()) {
    stop("HDF5Matrix is closed or invalid")
  }
  
  dims <- dim(x)
  size_bytes <- dims[1] * dims[2] * 8  # Assume numeric/double
  
  size <- switch(unit,
    "bytes" = size_bytes,
    "KB"    = size_bytes / 1024,
    "MB"    = size_bytes / (1024^2),
    "GB"    = size_bytes / (1024^3)
  )
  
  return(size)
}
