# S3_pseudoinverse.R
#
# S3 generic pseudoinverse() and HDF5Matrix method.
# pseudoinverse() is not in base R, so we define the generic here.

#' Moore-Penrose pseudoinverse
#'
#' @description
#' Generic function for computing the Moore-Penrose pseudoinverse.
#' The \code{HDF5Matrix} method computes the pseudoinverse entirely on
#' disk using block-wise SVD; the full matrix is never loaded into RAM.
#'
#' Delegates to \code{bdpseudoinv_hdf5()}. Result stored in the same
#' HDF5 file under \code{OUTPUT/<dataset>_pinv} by default.
#'
#' @param x   An object. For \code{HDF5Matrix}, the dataset stored on disk.
#' @param \dots Additional arguments passed to \code{x$pseudoinverse()}:
#'   \code{outgroup}, \code{outdataset}, \code{overwrite}, \code{threads},
#'   \code{compression}.
#' @return For \code{HDF5Matrix}: a new \code{HDF5Matrix} containing the
#'   pseudoinverse.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' m   <- matrix(c(1,2,3,4,5,6), 3, 2)
#' X   <- hdf5_create_matrix(tmp, "data/A", data = m)
#' P   <- pseudoinverse(X)
#' dim(P)   # 2 x 3
#' close(X); close(P)
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{solve.HDF5Matrix}}, \code{\link{svd.HDF5Matrix}}
#' @export
pseudoinverse <- function(x, ...) UseMethod("pseudoinverse")

#' @rdname pseudoinverse
#' @export
pseudoinverse.HDF5Matrix <- function(x, ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")
    x$pseudoinverse(...)
}

#' @rdname pseudoinverse
#' @export
pseudoinverse.default <- function(x, ...) {
    stop("pseudoinverse() is not implemented for objects of class '",
         class(x)[1], "'")
}
