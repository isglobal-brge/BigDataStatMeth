# S3 methods - Core
# dim, print, str for HDF5Matrix objects

#' @name HDF5Matrix-S3
#' @title S3 methods for HDF5Matrix
#' @description Standard R generic methods for \code{HDF5Matrix} objects,
#'   allowing them to be used identically to in-memory matrices.
NULL


#' Dimensions of an HDF5Matrix
#'
#' @param x An \code{HDF5Matrix} object
#' @return Integer vector \code{c(nrows, ncols)}
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' 
#' X  <- hdf5_create_matrix(tmp, "data/matrix", data = matrix(rnorm(100), 10, 10))
#' X <- hdf5_matrix(tmp, "data/matrix")
#'
#' dim(X)    # c(10, 10)
#' nrow(X)   # 10
#' ncol(X)   # 10
#'
#' X$close()
#' unlink(tmp)
#' }
#'
#' @export
dim.HDF5Matrix <- function(x) {
  x$dim()
}

#' Length of an HDF5Matrix
#'
#' @description
#' Returns the total number of elements in an \code{HDF5Matrix} object,
#' defined as \code{prod(dim(x))} — consistent with the behaviour of
#' \code{base::length()} for ordinary R matrices.
#'
#' @param x An \code{HDF5Matrix} object.
#' @return A single integer: \code{nrow(x) * ncol(x)}.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' X <- hdf5_create_matrix(tmp, "data/m", data = matrix(1:100, 10, 10))
#' length(X)   # 100
#' close(X); unlink(tmp)
#' }
#'
#' @export
length.HDF5Matrix <- function(x) {
    prod(dim(x))
}

#' Print an HDF5Matrix object
#'
#' @param x An \code{HDF5Matrix} object
#' @param ... Ignored
#' @return Invisible \code{x}
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' 
#' X  <- hdf5_create_matrix(tmp, "data/matrix", data = matrix(rnorm(100), 10, 10))
#' X <- hdf5_matrix(tmp, "data/matrix")
#'
#' print(X)
#' X   # same as print(X)
#'
#' X$close()
#' unlink(tmp)
#' }
#'
#' @export
print.HDF5Matrix <- function(x, ...) {
  x$print()
}


#' Structure of an HDF5Matrix object
#'
#' @param object An \code{HDF5Matrix} object
#' @param ... Ignored
#' @return Invisible \code{object}
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' 
#' X  <- hdf5_create_matrix(tmp, "data/matrix", data = matrix(rnorm(100), 10, 10))
#' X <- hdf5_matrix(tmp, "data/matrix")
#'
#' str(X)
#'
#' X$close()
#' unlink(tmp)
#' }
#'
#' @export
str.HDF5Matrix <- function(object, ...) {
  cat("HDF5Matrix: ")
  if (object$is_valid()) {
    dims <- dim(object)
    info <- object$info()
    cat(dims[1], "x", dims[2], info$type, "\n")
    cat("  File:", info$filename, "\n")
    cat("  Path:", paste0(info$group, "/", info$dataset), "\n")
  } else {
    cat("CLOSED\n")
  }
  invisible(object)
}
# S3 close() method for HDF5Matrix
# Added to S3_core.R


#' Close HDF5Matrix
#'
#' @description
#' Close an HDF5Matrix object and release file resources immediately.
#' This is the standard R interface for resource cleanup.
#'
#' @param con An \code{HDF5Matrix} object
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisible \code{NULL}
#'
#' @details
#' Closes the HDF5 dataset without waiting for garbage collection.
#' After calling \code{close()}, \code{is_valid()} returns \code{FALSE} and
#' any further operations on the object will fail. The file is immediately
#' accessible by other tools such as HDFView.
#'
#' **Both syntaxes work:**
#' \itemize{
#'   \item \code{close(X)} - Standard R generic (recommended)
#'   \item \code{X$close()} - R6 method (still supported)
#' }
#'
#' For emergency closure of all open HDF5 objects in the session, see
#' \code{\link{hdf5_close_all}}.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#'
#' X  <- hdf5_create_matrix(tmp, "data/matrix", data = matrix(rnorm(100), 10, 10))
#' 
#' # Open matrix
#' X <- hdf5_matrix(tmp, "data/matrix")
#' data <- X[1:5, 1:5]
#'
#' # Close using S3 method (recommended)
#' close(X)
#'
#' # Or using R6 method (still works)
#' # X$close()
#'
#' X$is_valid()  # FALSE
#'
#' unlink(tmp)
#' }
#'
#' @seealso 
#' \code{\link{hdf5_matrix}} for opening datasets,
#' \code{\link{hdf5_close_all}} for closing all open objects
#'
#' @export
close.HDF5Matrix <- function(con, ...) {
  if (!inherits(con, "HDF5Matrix")) {
    stop("Object is not an HDF5Matrix")
  }
  
  # Call R6 close() method
  con$close()
  
  invisible(NULL)
}


#' Check if HDF5Matrix is open
#'
#' @description
#' Check whether an HDF5Matrix object is still valid and open.
#'
#' @param x An \code{HDF5Matrix} object
#'
#' @return Logical. \code{TRUE} if object is valid and open, \code{FALSE} otherwise.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' 
#' X  <- hdf5_create_matrix(tmp, "data/matrix", data = matrix(rnorm(100), 10, 10))
#'
#' X <- hdf5_matrix(tmp, "data/matrix")
#' is_open(X)  # TRUE
#'
#' close(X)
#' is_open(X)  # FALSE
#'
#' unlink(tmp)
#' }
#'
#' @export
is_open <- function(x) {
  UseMethod("is_open")
}

#' @export
is_open.HDF5Matrix <- function(x) {
  x$is_valid()
}
