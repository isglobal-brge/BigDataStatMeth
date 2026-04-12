# S3_eigen.R
#
# S3 generic override and HDF5Matrix method for eigen().
#
# Design note:
#   base::eigen(x, symmetric, only.values, EISPACK) has NO "..." in its
#   formals — extra arguments like k, which, compute_vectors, threads
#   cause "unused argument" errors when dispatched through base::eigen.
#   The generic must therefore be overridden, exactly as done for qr()
#   in S3_factorizations.R (Phase 9).
#
#   The .default method mirrors the base::eigen signature and forwards
#   to base::eigen for all non-HDF5Matrix objects.


#' Spectral decomposition
#'
#' @description
#' Overrides \code{base::eigen()} to dispatch on \code{\link{HDF5Matrix}}
#' objects. For plain R matrices the call is forwarded to
#' \code{base::eigen()}.
#'
#' @param x          An \code{\link{HDF5Matrix}} or any R object accepted
#'   by \code{base::eigen()}.
#' @param symmetric  Logical. Whether to assume \code{x} is symmetric.
#'   Passed to \code{base::eigen()} for plain R objects.
#' @param \dots      For \code{HDF5Matrix}: additional arguments forwarded
#'   to \code{x$eigen()} — \code{k}, \code{which}, \code{compute_vectors},
#'   \code{tolerance}, \code{max_iter}, \code{overwrite}, \code{threads}.
#'   Ignored for \code{eigen.default}.
#'
#' @return For \code{HDF5Matrix}: a named list with elements
#'   \code{values} (numeric vector), \code{vectors}
#'   (\code{HDF5Matrix} or \code{NULL}), \code{values_imag}, and
#'   \code{is_symmetric}.
#'   For other objects: the result of \code{base::eigen()}.
#'
#' @seealso \code{\link{svd.HDF5Matrix}}, \code{\link{prcomp.HDF5Matrix}}
#' @export
eigen <- function(x, symmetric, ...) UseMethod("eigen")

#' @rdname eigen
#' @export
eigen.default <- function(x, symmetric = !isSymmetric(x),
                           only.values = FALSE, EISPACK = FALSE, ...) {
    base::eigen(x, symmetric = symmetric,
                only.values = only.values, EISPACK = EISPACK)
}


#' @rdname eigen
#' @description
#' Computes eigenvalues (and optionally eigenvectors) of an
#' \code{\link{HDF5Matrix}} stored on disk. All computation is block-wise;
#' the full matrix is never loaded into RAM.
#'
#' Delegates to \code{bdEigen_hdf5()} (Spectra / Eigen backend). Output
#' datasets are written to the same HDF5 file under
#' \code{EIGEN/<dataset>/values} and \code{EIGEN/<dataset>/vectors}.
#'
#' @param k               Integer or \code{NULL}. Number of eigenvalues
#'   to compute (default: all).
#' @param which           Character. Which eigenvalues to return:
#'   \code{"LM"} (largest magnitude, default), \code{"SM"}, \code{"LR"},
#'   \code{"SR"}, \code{"LI"}, \code{"SI"}.
#' @param compute_vectors Logical. Also compute eigenvectors
#'   (default \code{TRUE}).
#' @param tolerance       Numeric or \code{NULL}. Convergence tolerance.
#' @param max_iter        Integer or \code{NULL}. Maximum Lanczos iterations.
#' @param overwrite       Logical. Overwrite existing output datasets.
#' @param threads         Integer or \code{NULL}. OpenMP threads.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' m   <- crossprod(matrix(rnorm(400), 20, 20))
#' X   <- hdf5_create_matrix(tmp, "data/M", data = m)
#' ev  <- eigen(X, symmetric = TRUE, k = 5L)
#' ev$values
#' close(ev$vectors)
#' close(X)
#' unlink(tmp)
#' }
#' @export
eigen.HDF5Matrix <- function(x, symmetric = TRUE, ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")
    x$eigen(symmetric = symmetric, ...)
}
