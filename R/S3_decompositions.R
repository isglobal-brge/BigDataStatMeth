# S3_decompositions.R
#
# S3 dispatch for svd() and prcomp() on HDF5Matrix objects.
#
# Both svd() and prcomp() are proper S3 generics in base R / stats,
# so HDF5Matrix methods are registered through the standard mechanism
# without any masking of the base generic.
#
#   svd()    → base::svd    (UseMethod in base)
#   prcomp() → stats::prcomp (UseMethod in stats)
#
# User-facing interface follows base R conventions:
#   svd(X, nu, nv, ...)  → list(d, u, v) where d is numeric, u/v are HDF5Matrix
#   prcomp(X, retx, center, scale., ...) → HDF5PCA object
#
# NOTE: nu and nv (number of left/right vectors to return) are passed
# as nev = max(nu, nv) to the underlying block-wise algorithm.  The
# caller can then close unwanted matrices if needed.


# ---------------------------------------------------------------------------
# Generic override  (base::cor has no UseMethod → must mask)
# ---------------------------------------------------------------------------


# svd() override — base::svd() has no '...' so extra args fail before S3 dispatch.
# We create our own generic that routes HDF5Matrix to svd.HDF5Matrix and
# everything else to base::svd().
#' @export
svd <- function(x, nu = min(dim(x)), nv = min(dim(x)), ...) UseMethod("svd")

#' @exportS3Method
svd.default <- function(x, nu = min(dim(x)), nv = min(dim(x)), ...) {
    base::svd(x, nu = nu, nv = nv)
}


# ── svd.HDF5Matrix ───────────────────────────────────────────────────────────

#' Singular Value Decomposition of an HDF5Matrix
#'
#' Block-wise SVD entirely on disk.  The matrix \code{x} is decomposed
#' into \code{x = u \%*\% diag(d) \%*\% t(v)}.
#'
#' Singular values \code{d} are loaded into a plain numeric vector
#' (they are always small: at most \code{min(nrow(x), ncol(x))} values).
#' \code{u} and \code{v} are returned as \code{HDF5Matrix} objects.
#'
#' @param x   An \code{HDF5Matrix} object.
#' @param nu  Number of left  singular vectors to compute (default = \code{min(dim(x))}).
#' @param nv  Number of right singular vectors to compute (default = \code{min(dim(x))}).
#' @param center Logical.  Center columns before decomposition (default \code{TRUE}).
#' @param scale  Logical.  Scale  columns before decomposition (default \code{TRUE}).
#' @param k      Number of local SVDs per incremental level (default 2).
#' @param q      Number of incremental levels (default 1).
#' @param method Computation method: \code{"auto"} (default), \code{"blocks"},
#'   or \code{"full"}.
#' @param rankthreshold Numeric in \code{[0, 0.1]}.  Rank approximation
#'   threshold (default 0).
#' @param overwrite Logical.  Overwrite existing SVD results (default \code{FALSE}).
#' @param threads Integer.  OpenMP threads (\code{-1} = auto-detect).
#' @param ... Ignored (S3 compatibility).
#'
#' @return Named list with:
#'   \describe{
#'     \item{\code{d}}{Numeric vector of non-negative singular values, decreasing.}
#'     \item{\code{u}}{HDF5Matrix of left  singular vectors, \code{nrow(x) × nu}.}
#'     \item{\code{v}}{HDF5Matrix of right singular vectors, \code{ncol(x) × nv}.}
#'   }
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' X   <- hdf5_create_matrix(tmp, "data/M", data = matrix(rnorm(500), 50, 10))
#' res <- svd(X)
#' length(res$d)   # 10  (min(50,10))
#' dim(res$u)      # 50 x 10
#' dim(res$v)      # 10 x 10
#' X$close()
#' res$u$close(); res$v$close()
#' unlink(tmp)
#' }
#'
#' @export
svd.HDF5Matrix <- function(x,
                            nu     = min(dim(x)),
                            nv     = min(dim(x)),
                            center = TRUE,
                            scale  = TRUE,
                            k      = 2L,
                            q      = 1L,
                            method = "auto",
                            rankthreshold = 0.0,
                            overwrite     = FALSE,
                            threads       = -1L,
                            ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")

    nu_int <- as.integer(nu)
    nv_int <- as.integer(nv)

    # NOTE: nev in the C++ block-SVD algorithm controls eigenvalues *per block*,
    # not the output rank.  The algorithm always writes min(m,n) singular triplets
    # to disk.  We compute the full decomposition and then truncate d/u/v here
    # in R so the caller gets exactly nu left vectors and nv right vectors,
    # matching base::svd() semantics.
    res <- x$svd(
        k             = k,
        q             = q,
        nev           = 0L,   # 0 = all; nev as block param is unused in C++
        center        = center,
        scale         = scale,
        method        = method,
        rankthreshold = rankthreshold,
        overwrite     = overwrite,
        threads       = threads
    )

    # ── Truncate d ───────────────────────────────────────────────────────────
    k_out <- max(nu_int, nv_int)
    d_full <- res$d
    if (k_out > 0L && k_out < length(d_full))
        res$d <- d_full[seq_len(k_out)]

    # ── Truncate U to nu columns ─────────────────────────────────────────────
    # u is a HDF5Matrix of shape nrow(x) × r where r = min(dim(x)).
    # When nu < r, select only the first nu columns via the existing [r,c] method.
    if (nu_int > 0L && nu_int < ncol(res$u)) {
        u_mat  <- res$u[, seq_len(nu_int)]    # reads into memory (nu cols only)
        u_path <- paste0(res$u$get_group(), "/", res$u$get_dataset(), "_nu")
        try(close(res$u), silent = TRUE)
        res$u  <- hdf5_create_matrix(x$get_filename(), u_path,
                                     data = u_mat, overwrite = TRUE)
    }

    # ── Truncate V to nv columns ─────────────────────────────────────────────
    if (nv_int > 0L && nv_int < ncol(res$v)) {
        v_mat  <- res$v[, seq_len(nv_int)]
        v_path <- paste0(res$v$get_group(), "/", res$v$get_dataset(), "_nv")
        try(close(res$v), silent = TRUE)
        res$v  <- hdf5_create_matrix(x$get_filename(), v_path,
                                     data = v_mat, overwrite = TRUE)
    }

    res
}


# ── prcomp.HDF5Matrix ────────────────────────────────────────────────────────

#' Principal Component Analysis of an HDF5Matrix
#'
#' Block-wise PCA entirely on disk, equivalent to \code{\link{prcomp}()}.
#' Implements the same interface as \code{stats::prcomp()} but operates on
#' data stored in an HDF5 file without loading it into RAM.
#'
#' @param x        An \code{HDF5Matrix} object.
#' @param retx     Logical.  If \code{TRUE} (default) return the individual
#'   coordinates (\code{x} slot).  If \code{FALSE} the \code{x} slot is
#'   \code{NULL} in the returned object.
#' @param center   Logical.  Subtract column means before PCA (default \code{TRUE}).
#' @param scale.   Logical.  Divide by column SDs before PCA (default \code{FALSE}).
#' @param tol      Ignored (present for interface compatibility with \code{prcomp()}).
#' @param ncomponents Integer.  Number of PCs to compute (0 = all, default).
#' @param k        Number of local SVDs per incremental level (default 2).
#' @param q        Number of incremental levels (default 1).
#' @param method   Computation method: \code{"auto"} (default), \code{"blocks"},
#'   or \code{"full"}.
#' @param rankthreshold Numeric in \code{[0, 0.1]}.  Rank approximation threshold.
#' @param svdgroup HDF5 group for intermediate SVD storage (default \code{"SVD/"}).
#' @param overwrite Logical.  Recompute even if PCA results exist (default \code{FALSE}).
#' @param threads  Integer.  OpenMP threads (\code{-1} = auto-detect).
#' @param ... Ignored (S3 compatibility).
#'
#' @return An object of class \code{c("HDF5PCA", "list")} with elements:
#'   \describe{
#'     \item{\code{sdev}}{Numeric vector.  Standard deviations of the PCs.}
#'     \item{\code{rotation}}{HDF5Matrix.  Variable loadings (rotation matrix).}
#'     \item{\code{x}}{HDF5Matrix or \code{NULL}.  Individual coordinates.}
#'     \item{\code{center}}{Logical.  Whether columns were centered.}
#'     \item{\code{scale}}{Logical.  Whether columns were scaled.}
#'     \item{\code{cumvar}}{Numeric vector.  Cumulative variance explained (percent).}
#'     \item{\code{lambda}}{Numeric vector.  Eigenvalues.}
#'     \item{\code{var.cos2}}{HDF5Matrix.  Squared cosines for variables.}
#'     \item{\code{ind.cos2}}{HDF5Matrix.  Squared cosines for individuals.}
#'     \item{\code{ind.contrib}}{HDF5Matrix.  Contributions of individuals to PCs.}
#'     \item{\code{file}}{Character.  Path to the HDF5 file with all results.}
#'   }
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' X   <- hdf5_create_matrix(tmp, "data/M", data = matrix(rnorm(1000), 100, 10))
#' pca <- prcomp(X, center = TRUE, scale. = FALSE)
#' cat("Variance explained (PC1-3):", pca$cumvar[1:3], "\n")
#' dim(pca$rotation)   # 10 x nPC
#' dim(pca$x)          # 100 x nPC
#' X$close()
#' pca$rotation$close(); pca$x$close()
#' unlink(tmp)
#' }
#'
#' @exportS3Method stats::prcomp HDF5Matrix
prcomp.HDF5Matrix <- function(x,
                               retx          = TRUE,
                               center        = TRUE,
                               scale.        = FALSE,
                               tol           = NULL,
                               rank.         = NULL,       # base R alias → ncomponents
                               ncomponents   = 0L,
                               k             = 2L,
                               q             = 1L,
                               method        = "auto",
                               rankthreshold = 0.0,
                               svdgroup      = "SVD/",
                               overwrite     = FALSE,
                               threads       = -1L,
                               ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")

    # rank. (base R prcomp convention) takes precedence over ncomponents
    if (!is.null(rank.))
        ncomponents <- as.integer(rank.)

    result <- x$pca(
        ncomponents   = ncomponents,
        center        = center,
        scale         = scale.,
        k             = k,
        q             = q,
        method        = method,
        rankthreshold = rankthreshold,
        svdgroup      = svdgroup,
        overwrite     = overwrite,
        threads       = threads
    )

    if (!isTRUE(retx)) result$x <- NULL

    result
}


# ── print method for HDF5PCA ─────────────────────────────────────────────────

#' Print method for HDF5PCA objects
#'
#' @param x   An \code{HDF5PCA} object returned by \code{prcomp()} or \code{$pca()}.
#' @param ... Ignored.
#' @exportS3Method base::print HDF5PCA
print.HDF5PCA <- function(x, ...) {
    cat("HDF5PCA — Principal Component Analysis (on-disk)\n")
    cat(sprintf("  File     : %s\n", x$file))
    cat(sprintf("  PCs      : %d\n", length(x$sdev)))
    cat(sprintf("  center   : %s  |  scale: %s\n",
                x$center, x$scale))
    cat(sprintf("  sdev[1:3]: %s\n",
                paste(round(head(x$sdev, 3), 4), collapse = ", ")))
    if (!is.null(x$cumvar))
        cat(sprintf("  cumvar%%  : %s ...\n",
                    paste(round(head(x$cumvar, 3), 2), collapse = ", ")))
    if (!is.null(x$rotation) && inherits(x$rotation, "HDF5Matrix"))
        cat(sprintf("  rotation : %d × %d  [HDF5Matrix]\n",
                    dim(x$rotation)[1], dim(x$rotation)[2]))
    if (!is.null(x$x) && inherits(x$x, "HDF5Matrix"))
        cat(sprintf("  x (ind.) : %d × %d  [HDF5Matrix]\n",
                    dim(x$x)[1], dim(x$x)[2]))
    invisible(x)
}
