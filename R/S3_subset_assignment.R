# S3 method - Subsetting assignment
# [<- for HDF5Matrix objects


#' Subsetting assignment for HDF5Matrix objects
#'
#' @param x An \code{HDF5Matrix} object
#' @param i Row indices (numeric, logical, or missing)
#' @param j Column indices (numeric, logical, or missing)
#' @param ... Ignored
#' @param value Values to assign (scalar, vector, or matrix)
#'
#' @return The modified \code{HDF5Matrix} object (invisibly)
#'
#' @details
#' Writes data to the HDF5 dataset backing the \code{HDF5Matrix} object.
#' Supports:
#' \itemize{
#'   \item Scalar assignment: \code{X[i, j] <- 5}
#'   \item Vector assignment: \code{X[i, ] <- c(1, 2, 3)}
#'   \item Matrix assignment: \code{X[1:3, 1:3] <- matrix(...)}
#'   \item Full replacement: \code{X[] <- matrix(...)}
#' }
#'
#' The value is automatically recycled or reshaped to match the target dimensions.
#' Changes are written immediately to disk.
#'
#' @examples
#' \donttest{
#' 
#' tmp <- tempfile(fileext = ".h5")
#' 
#' # Create a matrix
#' X  <- hdf5_create_matrix(tmp, "data/X", data = matrix(rnorm(100), 10, 10))
#' 
#' X <- hdf5_matrix(tmp, "data/X")
#'
#' # Assign scalar
#' X[1, 1] <- 42
#'
#' # Assign row
#' X[2, ] <- 1:10
#'
#' # Assign block
#' X[1:3, 1:3] <- matrix(0, 3, 3)
#' 
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @export
`[<-.HDF5Matrix` <- function(x, i, j, ..., value) {
  if (!x$is_valid()) {
    stop("HDF5Matrix object is closed or invalid")
  }

  dims <- x$dim()
  nrows <- dims[1L]
  ncols <- dims[2L]

  # Determine indices
  i_missing <- missing(i)
  j_missing <- missing(j)

  # Case 1: X[] <- value (replace entire matrix)
  if (i_missing && j_missing) {
    # Value must be a matrix with matching dimensions or a scalar
    if (length(value) == 1L) {
      # Scalar: create matrix filled with scalar
      value <- matrix(value, nrow = nrows, ncol = ncols)
    } else if (is.matrix(value)) {
      if (nrow(value) != nrows || ncol(value) != ncols) {
        stop(sprintf("Replacement matrix must be %d x %d", nrows, ncols))
      }
    } else {
      # Vector: recycle to matrix
      if (length(value) != nrows * ncols) {
        stop(sprintf("Replacement vector must have length %d", nrows * ncols))
      }
      value <- matrix(value, nrow = nrows, ncol = ncols, byrow = FALSE)
    }

    rcpp_hdf5dataset_write_all(x$.__enclos_env__$private$ptr, value)

    # Write dimnames if value carries them
    dn <- dimnames(value)
    if (!is.null(dn)) {
      info <- x$info()
      rn   <- if (!is.null(dn[[1]])) dn[[1]] else character(0)
      cn   <- if (!is.null(dn[[2]])) dn[[2]] else character(0)
      if (length(rn) > 0 || length(cn) > 0)
        # bdWrite_hdf5_dimnames(info$filename, info$group, info$dataset, rn, cn)
          rcpp_hdf5dataset_write_dimnames(x$.__enclos_env__$private$ptr, rn, cn)
    }
    return(invisible(x))
  }

  # Parse indices
  if (i_missing) {
    i_seq <- seq_len(nrows)
  } else if (is.logical(i)) {
    if (length(i) != nrows) {
      stop("Logical index for rows must have length ", nrows)
    }
    i_seq <- which(i)
  } else if (is.numeric(i)) {
    i_seq <- as.integer(i)
    if (any(i_seq < 1L | i_seq > nrows)) {
      stop("Row indices out of bounds")
    }
  } else {
    stop("Row indices must be numeric, logical, or missing")
  }

  if (j_missing) {
    j_seq <- seq_len(ncols)
  } else if (is.logical(j)) {
    if (length(j) != ncols) {
      stop("Logical index for columns must have length ", ncols)
    }
    j_seq <- which(j)
  } else if (is.numeric(j)) {
    j_seq <- as.integer(j)
    if (any(j_seq < 1L | j_seq > ncols)) {
      stop("Column indices out of bounds")
    }
  } else {
    stop("Column indices must be numeric, logical, or missing")
  }

  n_i <- length(i_seq)
  n_j <- length(j_seq)

  # Prepare value matrix
  if (length(value) == 1L) {
    # Scalar: replicate to match dimensions
    value_mat <- matrix(value, nrow = n_i, ncol = n_j)
  } else if (is.matrix(value)) {
    if (nrow(value) != n_i || ncol(value) != n_j) {
      stop(sprintf("Replacement matrix must be %d x %d", n_i, n_j))
    }
    value_mat <- value
  } else {
    # Vector: recycle to matrix
    expected_length <- n_i * n_j
    if (length(value) != expected_length) {
      # R recycles vectors, but be explicit about dimension match
      if (expected_length %% length(value) != 0L) {
        warning("Length of replacement value is not a multiple of target size")
      }
      value <- rep_len(value, expected_length)
    }
    value_mat <- matrix(value, nrow = n_i, ncol = n_j, byrow = FALSE)
  }

  # Check if indices are contiguous (can use single write_block call)
  i_contiguous <- (length(i_seq) == 0L) || 
                  (length(i_seq) == 1L) ||
                  (all(diff(i_seq) == 1L))
  j_contiguous <- (length(j_seq) == 0L) || 
                  (length(j_seq) == 1L) ||
                  (all(diff(j_seq) == 1L))

  if (i_contiguous && j_contiguous) {
    # Single contiguous block - efficient
    row_offset <- if (length(i_seq) > 0L) min(i_seq) else 1L
    col_offset <- if (length(j_seq) > 0L) min(j_seq) else 1L

    rcpp_hdf5dataset_write_block(
      x$.__enclos_env__$private$ptr,
      value_mat,
      row_offset, col_offset,
      n_i, n_j
    )
  } else {
    # Non-contiguous indices - write element by element or row by row
    # This is less efficient but handles arbitrary index patterns

    for (idx_i in seq_along(i_seq)) {
      for (idx_j in seq_along(j_seq)) {
        rcpp_hdf5dataset_write_block(
          x$.__enclos_env__$private$ptr,
          matrix(value_mat[idx_i, idx_j], 1, 1),
          i_seq[idx_i], j_seq[idx_j],
          1L, 1L
        )
      }
    }
  }

  invisible(x)
}
