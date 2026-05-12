# S3_factorizations.R
#
# S3 generic overrides and methods for QR, Cholesky and matrix inverse of
# HDF5Matrix objects (Phase 9).
#
# Design notes:
#   - qr()    : base::qr(x, tol, LAPACK, ...) has NO "..." in the formals →
#               must override the generic exactly as done for svd() in Phase 8.
#   - chol()  : base::chol(x, pivot, ...) has "..." → no override needed.
#   - solve() : base::solve(a, b, ...) has "..." → no override needed.
#               For HDF5Matrix we only implement solve(a) (no rhs b).


# ── qr() — must override because base::qr has no ... ───────────────────────

#' QR decomposition of an HDF5Matrix
#'
#' @description
#' Overrides \code{base::qr()} to dispatch on \code{HDF5Matrix} objects.
#' For plain R matrices the call is forwarded to \code{base::qr()}.
#'
#' @param x   An \code{HDF5Matrix} or a plain R matrix/vector.
#' @param ... Additional arguments passed to \code{qr.HDF5Matrix()} or
#'   ignored for \code{qr.default()}.
#' @return For \code{HDF5Matrix}: a named list with elements \code{Q}
#'   and \code{R} (both \code{HDF5Matrix}).  For plain R objects: the
#'   result of \code{base::qr()}.
#'
#' @seealso \code{\link{qr.HDF5Matrix}}, \code{\link{chol.HDF5Matrix}},
#'   \code{\link{solve.HDF5Matrix}}
#'
#' @export
qr <- function(x, ...) UseMethod("qr")

#' @exportS3Method
qr.default <- function(x, tol = 1e-07, LAPACK = FALSE, ...) {
    base::qr.default(x, tol = tol, LAPACK = LAPACK)
}
# qr.default <- function(x, tol = 1e-07, LAPACK = FALSE, ...) {
#     # Avoid dispatch loop: call the C-level qr directly
#     stats::qr(x, tol = tol, LAPACK = LAPACK)  # si está en stats
# }


#' QR decomposition of an HDF5Matrix
#'
#' @description
#' Computes A = Q R block-wise on disk and returns Q and R as
#' \code{HDF5Matrix} objects.
#'
#' @param x          An \code{HDF5Matrix}.
#' @param thin       Logical. Compute thin (economy) QR.  Default \code{FALSE}.
#'   For \code{method = "tsqr"} this is strongly recommended; full Q via TSQR
#'   requires an expensive basis-completion step (O(m² k)).
#' @param block_size Integer or NULL. Row-block size hint for TSQR; ignored by
#'   \code{method = "lapack"}.  NULL = auto.
#' @param overwrite  Logical. Overwrite existing results.  Default \code{FALSE}.
#' @param threads    Integer. OpenMP threads (-1 = auto, CRAN-compliant).
#'   For \code{method = "lapack"} controls Eigen BLAS threads.
#'   For \code{method = "tsqr"} controls the OpenMP block loop.
#' @param method     Character. Algorithm selection:
#'   \describe{
#'     \item{"auto"}{Default. Selects TSQR when m > 5n AND m > 1000
#'       (tall-skinny), LAPACK otherwise.}
#'     \item{"lapack"}{Eigen HouseholderQR. Reliable for any matrix shape.}
#'     \item{"tsqr"}{Parallel Tall-Skinny QR. Requires m >= n.
#'       Efficient for big-omics tall matrices (e.g. samples × features).
#'       The R factor is mathematically equivalent to LAPACK R (up to sign).
#'       The Q factor is a valid orthogonal matrix but differs from LAPACK Q.}
#'   }
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result datasets.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 6).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @param ...        Ignored (for S3 compatibility).
#' @return Named list: \code{Q} (\code{HDF5Matrix}), \code{R} (\code{HDF5Matrix}).
#'
#' @examples
#' \donttest{
#' 
#' tmp <- tempfile(fileext = ".h5")
#' 
#' X  <- hdf5_create_matrix(tmp, "data/A", data = matrix(rnorm(10000), 100, 100))
#' 
#' X   <- hdf5_matrix(tmp, "data/A")
#'
#' # Default (auto method)
#' res <- qr(X)
#' dim(res$Q)   # m x m  (or m x min(m,n) if thin = TRUE)
#' dim(res$R)   # m x n  (or min(m,n) x n)
#'
#' # Explicit TSQR for a tall-skinny matrix (recommended: thin = TRUE)
# Explicit TSQR for a tall-skinny matrix (recommended: thin = TRUE)
# X_tall <- hdf5_create_matrix(tmp, "data/B", data = matrix(rnorm(200000), 2000, 100))
# X_tall <- hdf5_matrix(tmp, "data/B")
# res_tsqr <- qr(X_tall, method = "tsqr", thin = TRUE, overwrite = TRUE)
# dim(res_tsqr$Q)   # 2000 x 10dim(res_tsqr$R)   # 100 x 100
#' 
#' hdf5_close_all()
#' unlink(tmp)
#' 
#' }
#'
#' @note 20260304: Added \code{method} parameter and TSQR support.
#'
#' @export
qr.HDF5Matrix <- function(x,
                           thin        = FALSE,
                           block_size  = NULL,
                           overwrite   = FALSE,
                           threads     = -1L,
                           method      = "auto",
                           compression = NULL,
                           ...) {
    compression <- .get_option("compression", default = NULL,
                                          override = compression)
    x$qr(thin        = thin,
         block_size  = block_size,
         overwrite   = overwrite,
         threads     = threads,
         method      = method,
         compression = compression)
}


# ── chol() ──────────────────────────────────────────────────────────────────

#' Cholesky decomposition of a symmetric positive-definite HDF5Matrix
#'
#' @description
#' Computes the lower-triangular Cholesky factor L such that A = L L'.
#' The input matrix must be square and symmetric positive-definite.
#'
#' @param x          An \code{HDF5Matrix}.
#' @param full_matrix Logical. Return full symmetric matrix (L + L').
#'   Default \code{FALSE}.
#' @param overwrite  Logical. Overwrite existing result.  Default \code{FALSE}.
#' @param threads    Integer. OpenMP threads (-1 = auto).
#' @param block_size Integer or NULL. Elements per block.  NULL = auto.
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result dataset.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 6).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @param ...        Ignored (for S3 compatibility).
#' @return \code{HDF5Matrix} containing the Cholesky factor L.
#'
#' @examples
#' \donttest{
#' 
#' tmp <- tempfile(fileext = ".h5")
#' 
#' X  <- hdf5_create_matrix(tmp, "data/X", data = matrix(rnorm(10000), 100, 100))
#' 
#' # Create a symmetric positive-definite matrix: A = t(X) %*% X
#' X  <- hdf5_matrix(tmp, "data/X")
#' AtA <- crossprod(X)              # HDF5Matrix, square SPD
#' L   <- chol(AtA)
#' 
#' hdf5_close_all()
#' unlink(tmp)
#' 
#' }
#'
#' @exportS3Method base::chol HDF5Matrix
chol.HDF5Matrix <- function(x,
                             full_matrix = FALSE,
                             overwrite   = FALSE,
                             threads     = -1L,
                             block_size  = NULL,
                             compression = NULL,
                             ...) {
    compression <- .get_option("compression", default = NULL,
                                          override = compression)
    x$chol(full_matrix = full_matrix,
           overwrite   = overwrite,
           threads     = threads,
           block_size  = block_size,
           compression = compression)
}


# ── solve() ─────────────────────────────────────────────────────────────────

#' Matrix inverse of a symmetric positive-definite HDF5Matrix via Cholesky
#'
#' @description
#' Computes the matrix inverse of a symmetric positive-definite HDF5Matrix
#' using Cholesky decomposition + back-substitution.  Equivalent to
#' \code{base::solve(A)} for SPD matrices.
#'
#' @param a          An \code{HDF5Matrix} (square, symmetric positive-definite).
#' @param b          Not supported for \code{HDF5Matrix}; must be missing.
#' @param full_matrix Logical. Return full symmetric inverse.  Default \code{TRUE}.
#' @param overwrite  Logical. Overwrite existing result.  Default \code{FALSE}.
#' @param threads    Integer. OpenMP threads (-1 = auto).
#' @param block_size Integer or NULL. Elements per block.  NULL = auto.
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result dataset.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 6).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @param ...        Ignored (for S3 compatibility).
#' @return \code{HDF5Matrix} containing the matrix inverse.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' 
#' X  <- hdf5_create_matrix(tmp, "data/X", data = matrix(rnorm(10000), 100, 100))
#' 
#' X   <- hdf5_matrix(tmp, "data/X")
#' AtA <- crossprod(X)              # HDF5Matrix, square SPD
#' inv <- solve(AtA)               # inverse of AtA
#' 
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @exportS3Method base::solve HDF5Matrix
solve.HDF5Matrix <- function(a,
                              b,
                              full_matrix = TRUE,
                              overwrite   = FALSE,
                              threads     = -1L,
                              block_size  = NULL,
                              compression = NULL,
                              ...) {
    if (!missing(b))
        stop("solve.HDF5Matrix: argument 'b' is not supported; ",
             "only matrix inversion (solve(A)) is implemented for HDF5Matrix.")
    compression <- .get_option("compression", default = NULL,
                                          override = compression)
    a$solve(full_matrix = full_matrix,
            overwrite   = overwrite,
            threads     = threads,
            block_size  = block_size,
            compression = compression)
}
