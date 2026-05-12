# S3 methods - Matrix algebra
# %*%, crossprod, tcrossprod for HDF5Matrix objects
# Updated to use global options from hdf5matrix_options()



# #' @rawNamespace export("%*%")
# #' @rawNamespace S3method("%*%",default)
# #' @rawNamespace S3method("%*%",HDF5Matrix)
# `%*%` <- function(x, y) 
#     UseMethod("%*%")

#' Matrix multiplication for HDF5Matrix
#'
#' @description
#' S3 generic for \code{\%*\%}. Dispatches to \code{\%*\%.HDF5Matrix}
#' for \code{HDF5Matrix} objects, and to \code{base::\%*\%} for all others.
#'
#' @param x Left-hand side matrix.
#' @param y Right-hand side matrix.
#' @return Result matrix.
#' @name %*%
#' @rdname matmul
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
#' \code{HDF5Matrix} for the R6 \code{$multiply()} method with explicit parameters
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



# -- crossprod.HDF5Matrix ------------------------------------------------------

#' Cross product of HDF5Matrix objects
#'
#' @description
#' S3 generic for \code{crossprod()}. Dispatches to
#' \code{\link{crossprod.HDF5Matrix}} for \code{HDF5Matrix} objects,
#' and to \code{base::crossprod()} for all others.
#'
#' @param x A matrix or \code{HDF5Matrix}.
#' @param y A matrix, \code{HDF5Matrix}, or \code{NULL}.
#' @param ... Additional arguments passed to the method.
#' @return Result of the cross product.
#' @name crossprod
#' @rdname crossprod.HDF5Matrix
#' @rawNamespace export(crossprod)
#' @rawNamespace S3method(crossprod, HDF5Matrix)
#' @rawNamespace S3method(crossprod, default)
crossprod <- function(x, y = NULL, ...) UseMethod("crossprod")

#' @noRd
crossprod.default <- function(x, y = NULL, ...) base::crossprod(x, y)

#' Cross product of HDF5Matrix objects
#'
#' @param x         An \code{HDF5Matrix} object.
#' @param y         An \code{HDF5Matrix} object, or \code{NULL} (default) to
#'   compute \code{t(x) \%*\% x}.
#' @param outgroup  Character or \code{NULL}. HDF5 group where the result is
#'   stored. Default \code{"OUTPUT"}.
#' @param outdataset Character or \code{NULL}. Dataset name for the result.
#'   Default \code{"CrossProd_x"} (single matrix) or
#'   \code{"CrossProd_x_x_y"} (two matrices).
#' @param ...       Ignored.
#' @return A new \code{HDF5Matrix} pointing to the result dataset.
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
#' fn <- tempfile(fileext = ".h5")
#' X <- hdf5_create_matrix(fn, "INPUT/X", data = matrix(rnorm(60), 6, 10))
#' Y <- hdf5_create_matrix(fn, "INPUT/Y", data = matrix(rnorm(60), 6, 10))
#'
#' # t(X) %*% X  → stored in OUTPUT/CrossProd_X
#' C1 <- crossprod(X)
#' dim(C1)
#'
#' # t(X) %*% Y  → stored in OUTPUT/CrossProd_X_x_Y
#' C2 <- crossprod(X, Y)
#'
#' # Custom output location
#' C3 <- crossprod(X, outgroup = "RESULTS", outdataset = "my_crossprod")
#'
#' hdf5_close_all()
#' unlink(fn)
#' }
#'
#' @seealso
#' \code{\link{hdf5matrix_options}} for global performance settings
#'
#' @rawNamespace S3method(crossprod, HDF5Matrix)
crossprod.HDF5Matrix <- function(x, y = NULL, outgroup = NULL, outdataset = NULL, ...) {
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
      x_ptr, y_ptr,
      paral       = .get_option("paral"),
      block_size  = .get_option("block_size"),
      threads     = .get_option("threads"),
      compression = .get_option("compression", default = NULL),
      outgroup    = outgroup,
      outdataset  = outdataset
  )
  
  hdf5_matrix(result_info$filename, result_info$path)
}



# -- tcrossprod.HDF5Matrix -----------------------------------------------------

#' Transposed cross product of HDF5Matrix objects
#'
#' @description
#' S3 generic for \code{tcrossprod()}. Dispatches to
#' \code{\link{tcrossprod.HDF5Matrix}} for \code{HDF5Matrix} objects,
#' and to \code{base::tcrossprod()} for all others.
#'
#' @param x A matrix or \code{HDF5Matrix}.
#' @param y A matrix, \code{HDF5Matrix}, or \code{NULL}.
#' @param ... Additional arguments passed to the method.
#' @return Result of the cross product.
#' @name tcrossprod
#' @rdname tcrossprod.HDF5Matrix
#' @rawNamespace export(tcrossprod)
#' @rawNamespace S3method(tcrossprod, HDF5Matrix)
#' @rawNamespace S3method(tcrossprod, default)
tcrossprod <- function(x, y = NULL, ...) UseMethod("tcrossprod")

#' @noRd
tcrossprod.default <- function(x, y = NULL, ...) base::tcrossprod(x, y)

#' Transposed cross product of HDF5Matrix objects
#'
#' @param x         An \code{HDF5Matrix} object.
#' @param y         An \code{HDF5Matrix} object, or \code{NULL} (default) to
#'   compute \code{x \%*\% t(x)}.
#' @param outgroup  Character or \code{NULL}. HDF5 group where the result is
#'   stored. Default \code{"OUTPUT"}.
#' @param outdataset Character or \code{NULL}. Dataset name for the result.
#'   Default \code{"tCrossProd_x"} (single matrix) or
#'   \code{"tCrossProd_x_x_y"} (two matrices).
#' @param ...       Ignored.
#' @return A new \code{HDF5Matrix} pointing to the result dataset.
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
#' fn <- tempfile(fileext = ".h5")
#' X <- hdf5_create_matrix(fn, "INPUT/X", data = matrix(rnorm(60), 6, 10))
#' Y <- hdf5_create_matrix(fn, "INPUT/Y", data = matrix(rnorm(60), 6, 10))
#'
#' # t(X) %*% X  → stored in OUTPUT/CrossProd_X
#' C1 <- tcrossprod(X)
#' dim(C1)
#'
#' # t(X) %*% Y  → stored in OUTPUT/CrossProd_X_x_Y
#' C2 <- tcrossprod(X, Y)
#'
#' # Custom output location
#' C3 <- tcrossprod(X, outgroup = "RESULTS", outdataset = "my_tcrossprod")
#'
#' hdf5_close_all()
#' unlink(fn)
#' }
#'
#' @seealso
#' \code{\link{hdf5matrix_options}} for global performance settings
#'
#' @rawNamespace S3method(tcrossprod, HDF5Matrix)
tcrossprod.HDF5Matrix <- function(x, y = NULL, outgroup = NULL, outdataset = NULL, ...) {
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
    compression = .get_option("compression", default = NULL),
    outgroup    = outgroup,
    outdataset  = outdataset
  )
  
  hdf5_matrix(result_info$filename, result_info$path)
}
