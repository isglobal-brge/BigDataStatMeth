# S3 methods - Matrix algebra
# %*%, crossprod, tcrossprod for HDF5Matrix objects
# Updated to use global options from hdf5matrix_options()



#' @rawNamespace export("%*%")
#' @rawNamespace S3method("%*%",default)
#' @rawNamespace S3method("%*%",HDF5Matrix)
`%*%` <- function(x, y) 
    UseMethod("%*%")

#' @noRd
`%*%.default` <- function(x, y) 
    base::`%*%`(x, y)

#' Matrix multiplication of HDF5Matrix objects
#'
#' @param x An \code{HDF5Matrix} object (left-hand side)
#' @param y An \code{HDF5Matrix} object (right-hand side)
#' @return A new \code{HDF5Matrix} containing the result
#'
#' @details
#' Uses the BigDataStatMeth block-wise multiplication algorithm with optional
#' OpenMP parallelization.
#'
#' **Performance settings:**
#'
#' This method uses global options set via \code{\link{hdf5matrix_options}}.
#' For explicit control over parallelization and block size, use the R6 method:
#' \code{x$multiply(y, paral = TRUE, block_size = 2000, threads = 8)}
#'
#' **Examples of global configuration:**
#' \preformatted{
#' # Enable parallelization for all operations
#' hdf5matrix_options(paral = TRUE, threads = 8)
#' C <- A \%*\% B  # Now uses 8 threads
#'
#' # Reduce memory usage for huge matrices
#' hdf5matrix_options(block_size = 1000)
#' D <- A \%*\% B  # Now uses smaller blocks
#'
#' # Reset to defaults
#' hdf5matrix_options(paral = NULL, block_size = NULL, threads = NULL)
#' }
#'
#' The result is stored in the same HDF5 file as \code{x}, under a temporary
#' dataset name. Use \code{$info()} to see the location.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' rhdf5::h5createFile(tmp)
#' rhdf5::h5createGroup(tmp, "data")
#' rhdf5::h5write(matrix(rnorm(100), 10, 10), tmp, "data/A")
#' rhdf5::h5write(matrix(rnorm(100), 10, 10), tmp, "data/B")
#'
#' A <- hdf5_matrix(tmp, "data/A")
#' B <- hdf5_matrix(tmp, "data/B")
#'
#' # Standard R syntax
#' C <- A %*% B
#' dim(C)  # 10 x 10
#'
#' # With performance tuning
#' hdf5matrix_options(paral = TRUE, threads = 4)
#' D <- A %*% B  # Uses 4 threads
#'
#' A$close(); B$close(); C$close(); D$close()
#' unlink(tmp)
#' }
#'
#' @seealso
#' \code{\link{hdf5matrix_options}} for global performance settings,
#' \code{\link{HDF5Matrix}} for the R6 \code{$multiply()} method with explicit parameters
#'
# #' @export
#' @noRd
`%*%.HDF5Matrix` <- function(x, y) {
  if (!inherits(y, "HDF5Matrix")) {
    stop("Both operands must be HDF5Matrix objects")
  }
  
  # Use global options - .get_option returns NULL if not set (auto-detect in C++)
  x$multiply(y, 
             transpose_self  = FALSE, 
             transpose_other = FALSE,
             paral      = .get_option("paral"),
             block_size = .get_option("block_size"),
             threads    = .get_option("threads"),
             compression = .get_option("compression", default = NULL)
             )
}


#' Cross product of HDF5Matrix objects
#'
#' @param x An \code{HDF5Matrix} object
#' @param y An \code{HDF5Matrix} object, or \code{NULL} (default) to compute
#'   \code{t(x) \%*\% x}
#' @param ... Ignored
#' @return A new \code{HDF5Matrix} containing the result
#'
#' @details
#' Computes \eqn{t(x) \times y} (or \eqn{t(x) \times x} when \code{y = NULL}).
#' Uses the dedicated BigDataStatMeth block-wise cross-product algorithm, which
#' is more efficient than explicitly computing \code{t(x) \%*\% y}.
#'
#' **Performance settings:**
#'
#' This method uses global options set via \code{\link{hdf5matrix_options}}.
#'
#' **Symmetric optimization:**
#'
#' When \code{y = NULL} or \code{y} refers to the same dataset as \code{x},
#' the symmetric optimisation (\code{bisSymetric = TRUE}) is applied automatically,
#' providing significant speedup.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' rhdf5::h5createFile(tmp)
#' rhdf5::h5createGroup(tmp, "data")
#' rhdf5::h5write(matrix(rnorm(60), 6, 10), tmp, "data/X")
#' rhdf5::h5write(matrix(rnorm(60), 6, 10), tmp, "data/Y")
#'
#' X <- hdf5_matrix(tmp, "data/X")
#' Y <- hdf5_matrix(tmp, "data/Y")
#'
#' # t(X) %*% X  (10 x 10, symmetric)
#' C1 <- crossprod(X)
#' dim(C1)
#'
#' # t(X) %*% Y  (10 x 10)
#' C2 <- crossprod(X, Y)
#' dim(C2)
#'
#' # With parallelization
#' hdf5matrix_options(paral = TRUE, threads = 8)
#' C3 <- crossprod(X)  # Uses 8 threads
#'
#' X$close(); Y$close(); C1$close(); C2$close(); C3$close()
#' unlink(tmp)
#' }
#'
#' @seealso
#' \code{\link{hdf5matrix_options}} for global performance settings
#'
#' @rawNamespace S3method(crossprod, HDF5Matrix)
crossprod.HDF5Matrix <- function(x, y = NULL, ...) {
  if (!x$is_valid()) stop("x is closed or invalid")

  if (is.null(y)) {
    y <- x
  } else {
    if (!inherits(y, "HDF5Matrix")) stop("y must be an HDF5Matrix object")
    if (!y$is_valid()) stop("y is closed or invalid")
  }

  x_ptr <- x$.__enclos_env__$private$ptr
  y_ptr <- y$.__enclos_env__$private$ptr

  # Use global options
  result_info <- rcpp_hdf5dataset_crossprod(
    x_ptr, 
    y_ptr, 
    paral       = .get_option("paral"),
    block_size  = .get_option("block_size"),
    threads     = .get_option("threads"),
    compression = .get_option("compression", default = NULL)
  )
  
  hdf5_matrix(result_info$filename, result_info$path)
}


#' Transposed cross product of HDF5Matrix objects
#'
#' @param x An \code{HDF5Matrix} object
#' @param y An \code{HDF5Matrix} object, or \code{NULL} (default) to compute
#'   \code{x \%*\% t(x)}
#' @param ... Ignored
#' @return A new \code{HDF5Matrix} containing the result
#'
#' @details
#' Computes \eqn{x \times t(y)} (or \eqn{x \times t(x)} when \code{y = NULL}).
#' Uses the dedicated BigDataStatMeth block-wise transposed cross-product
#' algorithm, which is more efficient than explicitly computing \code{x \%*\% t(y)}.
#'
#' **Performance settings:**
#'
#' This method uses global options set via \code{\link{hdf5matrix_options}}.
#'
#' **Symmetric optimization:**
#'
#' When \code{y = NULL} or \code{y} refers to the same dataset as \code{x},
#' the symmetric optimisation is applied automatically, providing significant speedup.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' rhdf5::h5createFile(tmp)
#' rhdf5::h5createGroup(tmp, "data")
#' rhdf5::h5write(matrix(rnorm(60), 6, 10), tmp, "data/X")
#' rhdf5::h5write(matrix(rnorm(40), 6, 4),  tmp, "data/Y")
#'
#' X <- hdf5_matrix(tmp, "data/X")
#' Y <- hdf5_matrix(tmp, "data/Y")
#'
#' # X %*% t(X)  (6 x 6, symmetric)
#' C1 <- tcrossprod(X)
#' dim(C1)
#'
#' # X %*% t(Y)  (6 x 6)
#' C2 <- tcrossprod(X, Y)
#' dim(C2)
#'
#' # With custom block size
#' hdf5matrix_options(block_size = 500)
#' C3 <- tcrossprod(X)  # Uses block_size = 500
#'
#' X$close(); Y$close(); C1$close(); C2$close(); C3$close()
#' unlink(tmp)
#' }
#'
#' @seealso
#' \code{\link{hdf5matrix_options}} for global performance settings
#'
#' @rawNamespace S3method(tcrossprod, HDF5Matrix)
tcrossprod.HDF5Matrix <- function(x, y = NULL, ...) {
  if (!x$is_valid()) stop("x is closed or invalid")

  if (is.null(y)) {
    y <- x
  } else {
    if (!inherits(y, "HDF5Matrix")) stop("y must be an HDF5Matrix object")
    if (!y$is_valid()) stop("y is closed or invalid")
  }

  x_ptr <- x$.__enclos_env__$private$ptr
  y_ptr <- y$.__enclos_env__$private$ptr

  # Use global options
  result_info <- rcpp_hdf5dataset_tcrossprod(
    x_ptr, 
    y_ptr,
    paral       = .get_option("paral"),
    block_size  = .get_option("block_size"),
    threads     = .get_option("threads"),
    compression = .get_option("compression", default = NULL)
  )
  
  hdf5_matrix(result_info$filename, result_info$path)
}
