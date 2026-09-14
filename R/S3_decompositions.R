# S3_decompositions.R
#
# S3 dispatch for svd() and prcomp() on HDF5Matrix objects.
#
# Both svd() and prcomp() are proper S3 generics in base R / stats,
# so HDF5Matrix methods are registered through the standard mechanism
# without any masking of the base generic.
#
#   svd() -> base::svd    (UseMethod in base)
#   prcomp() -> stats::prcomp (UseMethod in stats)
#
# User-facing interface follows base R conventions:
#   svd(X, nu, nv, ...)  -> list(d, u, v) where d is numeric, u/v are HDF5Matrix
#   prcomp(X, retx, center, scale., ...) -> HDF5PCA object
#
# NOTE: nu and nv (number of left/right vectors to return) are passed
# as nev = max(nu, nv) to the underlying block-wise algorithm.  The
# caller can then close unwanted matrices if needed.


# ---------------------------------------------------------------------------
# Generic override  (base::cor has no UseMethod -> must mask)
# ---------------------------------------------------------------------------


# svd() override - base::svd() has no '...' so extra args fail before S3 dispatch.
# We create our own generic that routes HDF5Matrix to svd.HDF5Matrix and
# everything else to base::svd().
# #' @export
# svd <- function(x, nu = min(dim(x)), nv = min(dim(x)), ...) UseMethod("svd")

#' Singular Value Decomposition (generic)
#'
#' @description
#' S3 generic for \code{svd()}. Dispatches to \code{\link{svd.HDF5Matrix}}
#' for \code{HDF5Matrix} objects, and to \code{base::svd()} for all others.
#'
#' @param x A matrix or \code{HDF5Matrix} object.
#' @param nu Number of left singular vectors to compute.
#' @param nv Number of right singular vectors to compute.
#' @param ... Additional arguments passed to the method.
#' @return Named list with components \code{d}, \code{u}, \code{v}.
#' @name svd
#' @rdname svd
#' @export
svd <- function(x, nu = min(dim(x)), nv = min(dim(x)), ...) UseMethod("svd")


#' @exportS3Method
svd.default <- function(x, nu = min(dim(x)), nv = min(dim(x)), ...) {
    base::svd(x, nu = nu, nv = nv)
}


# -- svd.HDF5Matrix ------------------------------------------------------------

#' Singular Value Decomposition of an HDF5Matrix
#'
#' Block-wise SVD entirely on disk.  The matrix \code{x} is decomposed
#' into \code{x = u \%*\% diag(d) \%*\% t(v)}.
#'
#' Singular values \code{d} are loaded into a plain numeric vector
#' (they are always small: at most \code{min(nrow(x), ncol(x))} values).
#' \code{u} and \code{v} are returned as \code{HDF5Matrix} objects.
#'
#' @details
#' \strong{Constant columns and \code{scale = TRUE}.} A column with zero
#' variance cannot be rescaled to unit variance, so the call stops with an
#' error reporting how many columns are affected and a few of their positions
#' -- the same contract as \code{stats::prcomp(x, scale. = TRUE)}. Previously
#' the division by zero propagated silently and every singular value came back
#' as \code{0}. Either pass \code{scale = FALSE} or remove the constant columns
#' first. \code{\link{scale.HDF5Matrix}} is deliberately unaffected: like
#' \code{base::scale()} it still returns \code{NaN} for such a column.
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
#' @param threads Integer.  OpenMP threads (\code{-1} = auto-detect).  The
#'   effective number may be lower than requested: the system ceiling
#'   (\code{OMP_NUM_THREADS}, \code{OMP_THREAD_LIMIT} and, by default, half of
#'   the detected cores) is applied silently by the OpenMP runtime.  Requesting
#'   more than that is not an error, but it now emits a warning stating how
#'   many threads can actually be used.
#' @param ... Ignored (S3 compatibility).
#'
#' @return Named list with:
#'   \describe{
#'     \item{\code{d}}{Numeric vector of non-negative singular values, decreasing.}
#'     \item{\code{u}}{HDF5Matrix of left  singular vectors, \code{nrow(x) x nu}.}
#'     \item{\code{v}}{HDF5Matrix of right singular vectors, \code{ncol(x) x nv}.}
#'   }
#'   The list carries attributes describing what was actually computed:
#'   \describe{
#'     \item{\code{method}}{\code{"full"} (exact LAPACK) or \code{"blocks"}
#'       (hierarchical block algorithm).}
#'     \item{\code{exact}}{\code{TRUE} when the exact LAPACK algorithm was
#'       used, i.e. \code{method == "full"} and no per-block truncation. It
#'       reports the \emph{algorithm}, not a measured accuracy: the block path
#'       at full rank is not exact by construction, though in our tests it
#'       agreed with LAPACK to ~1e-15 (see below).}
#'     \item{\code{nev}}{Per-block truncation rank actually applied
#'       (\code{0} = none). Requesting fewer vectors below the boundary does
#'       not truncate anything: the exact path runs and the reduction happens
#'       afterwards in R, losslessly.}
#'     \item{\code{truncated}}{\code{TRUE} when per-block truncation was
#'       applied, which is the dominant source of error (see below).}
#'     \item{\code{elements}}{\code{nrow(x) * ncol(x)}, the quantity compared
#'       against the boundary.}
#'     \item{\code{auto_threshold}}{The boundary itself, see
#'       \code{\link{svd_auto_threshold}}.}
#'     \item{\code{blocking}}{The \code{k} and \code{q} actually used.}
#'     \item{\code{rank}}{Number of singular values returned.}
#'   }
#'
#' @section Exact versus approximate decomposition:
#' This is the most important behavioural difference from \code{base::svd()}.
#' \strong{Two independent things} can make the result approximate, and they
#' are not equally dangerous.
#'
#' \strong{1. The algorithm (minor).} \code{method = "full"} reads the whole
#' matrix into RAM and computes an exact LAPACK SVD, needing
#' \code{nrow * ncol * 8} bytes. \code{method = "blocks"} computes a
#' hierarchical block SVD entirely on disk. \code{method = "auto"} (the
#' default) picks \code{"full"} while
#' \code{nrow * ncol < svd_auto_threshold()} (currently 53,687,091 elements,
#' about 410 MB of doubles) and \code{"blocks"} at or above it. When the full
#' rank is requested, the block path is essentially exact: relative errors of
#' order 1e-15 across the whole spectrum, independent of \code{k} and
#' \code{q}, in our measurements. Crossing the boundary is, on its own, not a
#' numerical event.
#'
#' \strong{2. Per-block rank truncation (major).} Asking for fewer singular
#' triplets than \code{min(dim(x))} -- \code{nu} or \code{nv} below full rank,
#' or \code{ncomponents}/\code{rank.} in \code{\link{prcomp.HDF5Matrix}} --
#' switches on truncation of every \emph{local} block SVD to that rank before
#' the blocks are merged. That is what actually costs accuracy, and it
#' compounds with blocking depth. On a 1600 x 200 matrix with a decaying
#' spectrum, relative error in the singular values:
#'
#' \tabular{lrrr}{
#'   \strong{requested rank} \tab \strong{k, q} \tab \strong{leading} \tab \strong{trailing} \cr
#'   200 (full) \tab 2, 1 \tab 3e-16 \tab 7e-15 \cr
#'   200 (full) \tab 4, 3 \tab 1e-15 \tab 4e-15 \cr
#'    50        \tab 2, 1 \tab 7e-06 \tab 8e-02 \cr
#'    50        \tab 2, 3 \tab 5e-05 \tab 1e-01 \cr
#'    10        \tab 2, 1 \tab 2e-03 \tab 1e-01 \cr
#'    10        \tab 2, 3 \tab 1e-02 \tab 2e-01
#' }
#'
#' Note that the error is not confined to the tail: at rank 10 even the
#' \emph{leading} singular value was wrong by 1-3\%. If you need a small number
#' of accurate components, ask for a generously oversampled rank and discard
#' the extra ones -- in the table above, computing 50 and keeping 10 is about
#' 250 times more accurate in the leading value than computing 10 directly --
#' or request the full rank when you can afford it.
#'
#' \strong{Detecting it.} Both switches are silent by default, so the result
#' records them: check \code{attr(res, "exact")}, \code{attr(res, "truncated")}
#' and \code{attr(res, "method")}. \code{svd()} also emits a \code{message()}
#' when \code{"auto"} falls through to the block path and when per-block
#' truncation is applied (silence with \code{suppressMessages()}).
#'
#' There is no cheap residual diagnostic: verifying a decomposition means
#' reconstructing \code{u \%*\% diag(d) \%*\% t(v)} and comparing it with
#' \code{x}, another full pass over the data. The practical checks are to
#' recompute with a larger requested rank, or with \code{method = "full"} on a
#' subset, and see whether the leading components move.
#'

#' There is no cheap residual diagnostic: verifying the decomposition means
#' reconstructing \code{u \%*\% diag(d) \%*\% t(v)} and comparing it with
#' \code{x}, which costs another full pass over the data. When accuracy
#' matters, the practical check is to recompute a small trailing portion with
#' \code{method = "full"} on a subset, or to compare across two values of
#' \code{q}.
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

    # NOTE: nev controls per-block truncation in the block-SVD algorithm.
    # Passing max(nu,nv) when the user requests fewer than full rank keeps
    # B_i = U_i * d_i compact, reducing the cost of the final joined SVD.
    # Pass nev to the block algorithm so per-block U matrices are truncated
    # to the requested rank, reducing the cost of the final joined SVD.
    # nev=0 means full rank (no truncation) — only truncate when the user
    # explicitly requests fewer vectors than the full decomposition.
    nev_blk <- if (max(nu_int, nv_int) < min(dim(x))) as.integer(max(nu_int, nv_int)) else 0L
    res <- x$svd(
        k             = k,
        q             = q,
        nev           = nev_blk,
        ##..## nev           = 0L,   # 0 = all; nev as block param is unused in C++
        center        = center,
        scale         = scale,
        method        = method,
        rankthreshold = rankthreshold,
        overwrite     = overwrite,
        threads       = threads
    )

    # -- Truncate d ------------------------------------------------------------
    k_out <- max(nu_int, nv_int)
    d_full <- res$d
    if (k_out > 0L && k_out < length(d_full))
        res$d <- d_full[seq_len(k_out)]

    # -- Truncate U to nu columns ----------------------------------------------
    # u is a HDF5Matrix of shape nrow(x) x r where r = min(dim(x)).
    # When nu < r, select only the first nu columns via the existing [r,c] method.
    if (nu_int > 0L && nu_int < ncol(res$u)) {
        u_mat  <- res$u[, seq_len(nu_int)]    # reads into memory (nu cols only)
        u_path <- paste0(res$u$get_group(), "/", res$u$get_dataset(), "_nu")
        try(close(res$u), silent = TRUE)
        res$u  <- hdf5_create_matrix(x$get_filename(), u_path,
                                     data = u_mat, overwrite = TRUE)
    }

    # -- Truncate V to nv columns ----------------------------------------------
    if (nv_int > 0L && nv_int < ncol(res$v)) {
        v_mat  <- res$v[, seq_len(nv_int)]
        v_path <- paste0(res$v$get_group(), "/", res$v$get_dataset(), "_nv")
        try(close(res$v), silent = TRUE)
        res$v  <- hdf5_create_matrix(x$get_filename(), v_path,
                                     data = v_mat, overwrite = TRUE)
    }

    # -- Path metadata ---------------------------------------------------------
    # $svd() already attached it; only 'rank' can have changed above.
    attr(res, "rank") <- length(res$d)

    res
}


# -- prcomp.HDF5Matrix ---------------------------------------------------------

#' Principal Component Analysis of an HDF5Matrix
#'
#' Block-wise PCA entirely on disk, equivalent to \code{\link{prcomp}()}.
#' Implements the same interface as \code{stats::prcomp()} but operates on
#' data stored in an HDF5 file without loading it into RAM.
#'
#' @details
#' \strong{Constant columns and \code{scale. = TRUE}.} As in
#' \code{stats::prcomp()}, a column with zero variance cannot be rescaled to
#' unit variance and the call stops with an error reporting how many columns
#' are affected and a few of their positions. Previously the division by zero
#' propagated silently and the whole decomposition degenerated to zeros. Either
#' leave \code{scale. = FALSE} (the default) or remove the constant columns
#' first.
#'
#' @param x        An \code{HDF5Matrix} object.
#' @param retx     Logical.  If \code{TRUE} (default) return the individual
#'   coordinates (\code{x} slot).  If \code{FALSE} the \code{x} slot is
#'   \code{NULL} in the returned object.
#' @param center   Logical.  Subtract column means before PCA (default \code{TRUE}).
#' @param scale.   Logical.  Divide by column SDs before PCA (default \code{FALSE}).
#' @param tol      Ignored (present for interface compatibility with \code{prcomp()}).
#' @param rank. Integer. Number of principal components to compute. When
#'   supplied (non-\code{NULL}), takes precedence over \code{ncomponents}.
#'   Present for compatibility with \code{stats::prcomp()}; unlike the base R
#'   method, it is not ignored here.
#' @param ncomponents Integer.  Number of PCs to compute (0 = all, default).
#' @param k        Number of local SVDs per incremental level (default 2).
#' @param q        Number of incremental levels (default 1).
#' @param method   Computation method: \code{"auto"} (default), \code{"blocks"},
#'   or \code{"full"}.
#' @param rankthreshold Numeric in \code{[0, 0.1]}.  Rank approximation threshold.
#' @param svdgroup HDF5 group for intermediate SVD storage (default \code{"SVD/"}).
#' @param overwrite Logical.  Recompute even if PCA results exist (default \code{FALSE}).
#' @param threads  Integer.  OpenMP threads (\code{-1} = auto-detect).  The
#'   effective number may be lower than requested: the system ceiling
#'   (\code{OMP_NUM_THREADS}, \code{OMP_THREAD_LIMIT} and, by default, half of
#'   the detected cores) is applied silently by the OpenMP runtime.  Requesting
#'   more than that is not an error, but it now emits a warning stating how
#'   many threads can actually be used.
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
#'   The object also carries the attributes \code{method}, \code{exact},
#'   \code{elements}, \code{auto_threshold}, \code{blocking} and \code{rank},
#'   recording which decomposition path was actually taken; \code{print()}
#'   shows them.
#'
#' @section Exact versus approximate decomposition:
#' PCA runs on the same machinery as \code{\link{svd.HDF5Matrix}} and inherits
#' both of its silent approximations. Read that help page before interpreting a
#' large PCA. In particular, supplying \code{ncomponents} or \code{rank.}
#' below \code{min(dim(x))} switches on per-block rank truncation, which is
#' the larger error source and can perturb even the first principal component
#' by a few percent; asking for a generously oversampled number of components
#' and discarding the extra ones is much more accurate than asking for exactly
#' the few you want. Check \code{attr(pca, "exact")} and
#' \code{attr(pca, "truncated")}; \code{print()} shows both.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' X   <- hdf5_create_matrix(tmp, "data/M", data = matrix(rnorm(1000), 100, 10))
#' pca <- prcomp(X, center = TRUE, scale. = FALSE)
#' cat("Variance explained (PC1-3):", pca$cumvar[1:3], "\n")
#' dim(pca$rotation)   # 10 x nPC
#' dim(pca$x)          # 100 x nPC
#'
#' # rank. takes precedence over ncomponents when both are supplied.
#' # NOTE: this reuses the same output location as the call above, so it
#' # must come after pca's results have already been read/used -- once
#' # overwrite = TRUE recreates the output datasets, any earlier
#' # HDF5Matrix handle pointing at that location (here, pca$rotation and
#' # pca$x) becomes invalid by design (see HDF5Matrix lifecycle, Section
#' # 5.3.3): the package proactively clears stale external pointers when
#' # overwriting to convert what would otherwise be a memory-unsafe crash
#' # into a clean R-level error.
#' 
#' pca5 <- prcomp(X, rank. = 5, ncomponents = 10, overwrite = TRUE)
#' dim(pca5$rotation)   # 10 x 5 -- rank. (5) was used, not ncomponents (10)
#' 
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @exportS3Method stats::prcomp HDF5Matrix
prcomp.HDF5Matrix <- function(x,
                               retx          = TRUE,
                               center        = TRUE,
                               scale.        = FALSE,
                               tol           = NULL,
                               rank.         = NULL,       # base R alias -> ncomponents
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


# -- print method for HDF5PCA --------------------------------------------------

#' Print method for HDF5PCA objects
#'
#' @param x   An \code{HDF5PCA} object returned by \code{prcomp()} or \code{$pca()}.
#' @param ... Ignored.
#' @exportS3Method base::print HDF5PCA
#' @importFrom utils head
print.HDF5PCA <- function(x, ...) {
    cat("HDF5PCA - Principal Component Analysis (on-disk)\n")
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
        cat(sprintf("  rotation : %d x %d  [HDF5Matrix]\n",
                    dim(x$rotation)[1], dim(x$rotation)[2]))
    if (!is.null(x$x) && inherits(x$x, "HDF5Matrix"))
        cat(sprintf("  x (ind.) : %d x %d  [HDF5Matrix]\n",
                    dim(x$x)[1], dim(x$x)[2]))
    if (!is.null(attr(x, "method"))) {
        cat(sprintf("  method   : %s (%s)\n",
                    attr(x, "method"),
                    if (isTRUE(attr(x, "exact"))) "exact" else "APPROXIMATE"))
        if (!isTRUE(attr(x, "exact")))
            cat(sprintf("             k = %d, q = %d%s\n",
                        attr(x, "blocking")[["k"]],
                        attr(x, "blocking")[["q"]],
                        if (isTRUE(attr(x, "truncated")))
                            sprintf(", per-block rank truncated to %d - see ?prcomp.HDF5Matrix",
                                    attr(x, "nev"))
                        else ""))
    }
    invisible(x)
}
