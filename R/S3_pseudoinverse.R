# S3_pseudoinverse.R
#
# S3 generic pseudoinverse() with HDF5Matrix and matrix methods.
# pseudoinverse() is not in base R, so we define the generic here.

#' Moore-Penrose pseudoinverse
#'
#' @description
#' Generic function for computing the Moore-Penrose pseudoinverse.
#'
#' The \code{HDF5Matrix} method delegates to \code{bdpseudoinv_hdf5()}.
#' This method reads the complete dataset into memory, computes the
#' pseudoinverse via a direct LAPACK SVD, and writes the result back to
#' HDF5. It is HDF5-backed but not block-wise: the full matrix and its
#' pseudoinverse must both fit in available RAM during the computation.
#' Result stored in the same HDF5 file under \code{OUTPUT/<dataset>_pinv}
#' by default.
#'
#' The \code{matrix} method delegates to \code{bdpseudoinv()} and computes
#' the pseudoinverse of an in-memory matrix, returning a plain R matrix.
#'
#' Both methods use the same singular-value tolerance (\eqn{10^{-9}}):
#' singular values at or below the tolerance are treated as zero,
#' following the standard Moore-Penrose construction
#' \eqn{A^+ = V \Sigma^+ U^\top}.
#'
#' @param x   An object: an \code{HDF5Matrix}, or a numeric \code{matrix}.
#' @param threads Optional integer. Number of threads for parallel computation
#'   when \code{x} is a plain \code{matrix}. If \code{NULL} (default), uses the
#'   maximum available threads. Ignored for \code{HDF5Matrix}, where the thread
#'   count is passed via \code{\\dots} instead (see below).
#' @param \\dots Additional arguments.
#'   For \code{HDF5Matrix}: passed to \code{x$pseudoinverse()}, i.e.
#'   \code{outgroup}, \code{outdataset}, \code{overwrite}, \code{threads},
#'   \code{compression}.
#' @return For \code{HDF5Matrix}: a new \code{HDF5Matrix} containing the
#'   pseudoinverse. For \code{matrix}: a plain numeric matrix.
#'
#' @examples
#' \donttest{
#' ## In-memory matrix
#' m <- matrix(c(1,2,3,4,5,6), 3, 2)
#' pseudoinverse(m)
#'
#' ## HDF5Matrix
#' tmp <- tempfile(fileext = ".h5")
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
pseudoinverse.matrix <- function(x, threads = NULL, ...) {
    bdpseudoinv(x, threads = threads)
}

#' @rdname pseudoinverse
#' @export
pseudoinverse.default <- function(x, ...) {
    stop("pseudoinverse() is not implemented for objects of class '",
         class(x)[1], "'")
}
