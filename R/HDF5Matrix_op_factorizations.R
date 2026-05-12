# HDF5Matrix_op_factorizations.R
#
# Adds $qr(), $chol() and $solve() methods to the HDF5Matrix R6 class.
#
# All methods follow the established pattern (Rule 1.7):
#   - Close own handle before calling C++ (string-based, writes to same file)
#   - C++ opens all handles internally (no conflicts)
#   - Reopen own handle after C++ finishes
#   - When overwrite=TRUE, also close any live HDF5Matrix handles pointing
#     to the output paths (rcpp_hdf5_close_at_paths) so stale result
#     objects are safely invalidated before their datasets are overwritten.
#
# $qr()    → list(Q = HDF5Matrix, R = HDF5Matrix)
# $chol()  → HDF5Matrix  (lower-triangular Cholesky factor L, A = L L')
# $solve() → HDF5Matrix  (matrix inverse via Cholesky, for SPD matrices)


# ── $qr() ───────────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "qr",
# @description
# Block-wise QR decomposition of an HDF5Matrix.
#'
# Computes A = Q R entirely on disk.  Both Q and R are returned as
# \code{HDF5Matrix} objects pointing to the result datasets on disk.
#'
# @param thin       Logical. If \code{TRUE}, compute thin (economy) QR:
#   Q is m × min(m,n), R is min(m,n) × n.  Default \code{FALSE} (full QR).
#   For method = "tsqr" thin = TRUE is strongly recommended; full Q via
#   TSQR requires an expensive basis-completion step.
# @param block_size Integer or NULL. Row-block size hint for TSQR; ignored
#   by method = "lapack".  NULL = auto.
# @param overwrite  Logical. Overwrite existing QR results.  Default \code{FALSE}.
#   When TRUE, any open HDF5Matrix handles pointing to the output paths
#   (<group>/QRDec/Q.<dataset>, <group>/QRDec/R.<dataset>) are automatically
#   closed and invalidated.
# @param threads    Integer. OpenMP threads (-1 = auto, CRAN-compliant).
#   For method = "lapack" controls Eigen::setNbThreads (BLAS internal).
#   For method = "tsqr"   controls the OpenMP parallel block loop.
# @param method     Character. Algorithm selection:
#   \describe{
#     \item{"auto"}{Default. Selects TSQR when m > 5n AND m > 1000
#       (tall-skinny matrices), LAPACK otherwise.}
#     \item{"lapack"}{Eigen HouseholderQR. Reliable for any shape.
#       Single-threaded Householder; thread count applies to Eigen BLAS.}
#     \item{"tsqr"}{Tall-Skinny QR. Parallel row-block algorithm.
#       Requires m >= n. Efficient for big-omics tall matrices.
#       The R factor is equivalent to LAPACK R (up to row sign flips).
#       The Q factor is a valid orthogonal matrix but differs from LAPACK Q.}
#   }
# @return Named list with elements:
#   \describe{
#     \item{Q}{HDF5Matrix. Orthogonal factor (m × k or m × m).}
#     \item{R}{HDF5Matrix. Upper-triangular factor (k × n or m × n).}
#   }
#
# @note 20260304: Added method parameter and TSQR support.
function(thin        = FALSE,
         block_size  = NULL,
         overwrite   = FALSE,
         threads     = -1L,
         method      = "auto",
         compression = NULL) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # Validate method
    method <- match.arg(method, c("auto", "lapack", "tsqr"))

    # When overwrite=TRUE, close any live handles pointing to output paths
    if (isTRUE(overwrite)) {
        out_group <- paste0(private$group, "/QRDec")
        rcpp_hdf5_close_at_paths(
            private$filename,
            c(paste0(out_group, "/Q.", private$dataset),
              paste0(out_group, "/R.", private$dataset))
        )
    }

    # ── Close own handle — C++ will open its own handles to the same file ──
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL

    res <- rcpp_hdf5dataset_qr(
        private$filename,
        private$group,
        private$dataset,
        thin       = isTRUE(thin),
        block_size  = if (is.null(block_size)) -1L else as.integer(block_size),
        overwrite   = isTRUE(overwrite),
        threads     = as.integer(threads),
        method      = method,
        compression = .get_option("compression", default = NULL, override = compression)
    )

    # ── Reopen after C++ has finished and released all its handles ──────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    list(
        Q = hdf5_matrix(res$file, res$path_Q),
        R = hdf5_matrix(res$file, res$path_R)
    )
})


# ── $chol() ─────────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "chol",
# @description
# Cholesky decomposition of a symmetric positive-definite HDF5Matrix.
#'
# Computes the lower-triangular factor L such that A = L L'.
# The input matrix must be square and symmetric positive-definite.
#'
# @param full_matrix Logical. If \code{TRUE}, return the full symmetric matrix
#   (L + L').  Default \code{FALSE} (lower-triangular L only).
# @param overwrite   Logical. Overwrite existing result.  Default \code{FALSE}.
#   When TRUE, any open HDF5Matrix handle pointing to the output path
#   (<group>/chol_<dataset>) is automatically closed and invalidated.
# @param threads     Integer. OpenMP threads (-1 = auto).
# @param block_size  Integer or NULL. Elements per block.  NULL = auto.
# @return \code{HDF5Matrix} containing the Cholesky factor L.
function(full_matrix = FALSE,
         overwrite   = FALSE,
         threads     = -1L,
         block_size  = NULL,
         compression = NULL) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # When overwrite=TRUE, close any live handles pointing to the output path
    if (isTRUE(overwrite)) {
        rcpp_hdf5_close_at_paths(
            private$filename,
            paste0(private$group, "/chol_", private$dataset)
        )
    }

    # ── Close own handle — C++ will open its own handles to the same file ──
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL

    res <- rcpp_hdf5dataset_chol(
        private$filename,
        private$group,
        private$dataset,
        full_matrix = isTRUE(full_matrix),
        overwrite   = isTRUE(overwrite),
        threads     = as.integer(threads),
        block_size  = if (is.null(block_size)) -1L else as.integer(block_size),
        compression = .get_option("compression", default = NULL, override = compression)
    )

    # ── Reopen after C++ has finished and released all its handles ──────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    hdf5_matrix(res$file, res$path_L)
})


# ── $solve() ────────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "solve",
# @description
# Matrix inverse of a symmetric positive-definite HDF5Matrix via Cholesky.
#'
# Equivalent to \code{base::solve(A)} for SPD matrices.  Uses Cholesky
# decomposition + back-substitution for numerical stability.
# The input matrix must be square and symmetric positive-definite.
#'
# @param full_matrix Logical. If \code{TRUE}, return the full symmetric inverse.
#   Default \code{TRUE}.
# @param overwrite   Logical. Overwrite existing result.  Default \code{FALSE}.
#   When TRUE, any open HDF5Matrix handle pointing to the output path
#   (<group>/inv_<dataset>) is automatically closed and invalidated.
# @param threads     Integer. OpenMP threads (-1 = auto).
# @param block_size  Integer or NULL. Elements per block.  NULL = auto.
# @return \code{HDF5Matrix} containing the matrix inverse.
function(full_matrix = TRUE,
         overwrite   = FALSE,
         threads     = -1L,
         block_size  = NULL,
         compression = NULL) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # When overwrite=TRUE, close any live handles pointing to the output path
    if (isTRUE(overwrite)) {
        rcpp_hdf5_close_at_paths(
            private$filename,
            paste0(private$group, "/inv_", private$dataset)
        )
    }

    # ── Close own handle — C++ will open its own handles to the same file ──
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL

    res <- rcpp_hdf5dataset_solve(
        private$filename,
        private$group,
        private$dataset,
        full_matrix = isTRUE(full_matrix),
        overwrite   = isTRUE(overwrite),
        threads     = as.integer(threads),
        block_size  = if (is.null(block_size)) -1L else as.integer(block_size),
        compression = .get_option("compression", default = NULL, override = compression)
    )

    # ── Reopen after C++ has finished and released all its handles ──────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    hdf5_matrix(res$file, res$path_inv)
})
