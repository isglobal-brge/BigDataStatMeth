# HDF5Matrix - Public constructor


#' Open an HDF5 dataset as an HDF5Matrix object
#'
#' @description
#' Opens an existing dataset in an HDF5 file and returns an \code{HDF5Matrix}
#' object that can be used with standard R syntax: \code{dim()}, \code{[i,j]},
#' \code{\%*\%}, \code{crossprod()}, and \code{tcrossprod()}.
#'
#' @param filename Path to an existing HDF5 file
#' @param path Full path to the dataset within the file, in the form
#'   \code{"group/dataset"} or \code{"group/subgroup/dataset"}
#'
#' @return An \code{HDF5Matrix} object
#'
#' @details
#' The HDF5 file remains open while the object exists, which improves
#' performance for repeated operations. Call \code{$close()} to release
#' the file lock explicitly, or use \code{rm()} and \code{gc()} for
#' automatic cleanup. In an emergency, \code{\link{hdf5_close_all}} closes
#' all open \code{HDF5Matrix} objects.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#'
#' # Create dataset using BigDataStatMeth API
#' X <- hdf5_create_matrix(tmp, "data/expression",
#'                          data = matrix(rnorm(200), 20, 10))
#'
#' dim(X)         # 20 x 10
#' X[1:5, 1:3]   # subset
#' crossprod(X)   # t(X) %*% X
#'
#' close(X)
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{hdf5_close_all}}, \code{\link{[.HDF5Matrix}},
#'   \code{\link{crossprod.HDF5Matrix}}, \code{\link{tcrossprod.HDF5Matrix}}
#'
#' @export
hdf5_matrix <- function(filename, path) {
  obj <- HDF5Matrix$new(filename, path)
  class(obj) <- c("HDF5Matrix", class(obj))
  obj
}
