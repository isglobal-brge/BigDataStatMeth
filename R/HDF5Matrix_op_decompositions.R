# HDF5Matrix_op_decompositions.R
#
# Adds $svd() and $pca() methods to the HDF5Matrix R6 class.
#
# Both methods follow the established pattern:
#   - Pass filename/group/dataset as plain strings to C++ (never SEXP ptr)
#   - C++ opens all datasets internally with unique_ptr (no handle conflicts)
#   - Results are returned as named lists of HDF5Matrix objects
#   - When overwrite=TRUE, close any live HDF5Matrix handles pointing to
#     the output paths (rcpp_hdf5_close_at_paths) so stale result objects
#     are safely invalidated before their datasets are overwritten.
#
# $svd() → list(d = numeric_vector, u = HDF5Matrix, v = HDF5Matrix)
#   d is loaded into memory: it is always small (min(m,n) values).
#   u and v remain on disk as HDF5Matrix objects.
#
# $pca() → object of class c("HDF5PCA", "list") with fields matching
#   base R prcomp() as closely as possible:
#     sdev, rotation (HDF5Matrix), x (HDF5Matrix),
#     center, scale, cumvar, lambda,
#     var.cos2 (HDF5Matrix), ind.cos2 (HDF5Matrix), ind.contrib (HDF5Matrix)


# ── SVD path metadata ────────────────────────────────────────────────────────
#
# svd()/prcomp() on an HDF5Matrix silently switch between two algorithms:
#
#   "full"   exact LAPACK SVD of the whole matrix held in RAM
#   "blocks" hierarchical (incremental) block SVD — an APPROXIMATION whose
#            accuracy degrades with blocking depth, most visibly in the
#            trailing singular values
#
# With method = "auto" the switch happens at an internal element-count
# boundary, so the same call can be exact or approximate depending only on the
# size of the input.  Downstream code cannot detect this from the result alone,
# which is what these helpers fix: the boundary is read from C++ (single source
# of truth) and the path actually taken is attached to the returned object.

#' Element-count boundary of the automatic SVD/PCA path selection
#'
#' @description
#' Number of matrix elements (\code{nrow * ncol}) at or above which
#' \code{method = "auto"} stops computing an exact decomposition and switches
#' to the hierarchical block approximation.
#'
#' @return Numeric scalar: the element-count boundary.
#'
#' @details
#' The value is read from the C++ constant \code{SVD_EXACT_MAX_ELEMENTS} so
#' that documentation and implementation cannot drift apart. It currently
#' corresponds to 53,687,091 elements, about 410 MB of double-precision data.
#'
#' @examples
#' svd_auto_threshold()
#'
#' @seealso \code{\link{svd.HDF5Matrix}}, \code{\link{prcomp.HDF5Matrix}}
#' @export
svd_auto_threshold <- function() rcpp_svd_auto_threshold()


# ── Effective thread count ──────────────────────────────────────────────────
#
# The C++ thread ceiling silently discards a requested count it cannot honour:
# OMP_THREAD_LIMIT and OMP_NUM_THREADS are respected, and by default only half
# of the detected cores are usable.  Asking for more than that used to be a
# no-op with no indication whatsoever, so `threads = 32` on a 4-core machine
# looked like it had been applied.  These two helpers make the requested value
# explicit and report the gap once, per call.

# Bring `threads` to the single integer the C++ layer expects.  A negative
# value (the documented -1 default) and NULL both mean "let the system decide".
.normalize_threads <- function(threads) {
    if (is.null(threads) || length(threads) != 1L) return(-1L)
    n <- suppressWarnings(as.integer(threads))
    if (is.na(n)) return(-1L)
    n
}

# Warn once when fewer threads than requested can actually be used.  Says
# nothing when no count was requested, so the default path stays silent.
.warn_threads_capped <- function(threads) {
    req <- .normalize_threads(threads)
    if (req <= 0L) return(invisible(NULL))

    eff <- rcpp_effective_threads(req)
    if (!is.na(eff) && eff < req)
        warning(sprintf(
            paste0("requested threads = %d, but only %d can be used on this ",
                   "system (limited by OMP_NUM_THREADS / available cores)"),
            req, eff), call. = FALSE)

    invisible(NULL)
}

# Attach the path metadata to a decomposition result, and warn once about each
# of the two things that can make it approximate without the user asking:
#
#   (1) method = "auto" falling through to the block algorithm above the
#       element-count boundary;
#   (2) per-block rank truncation, which is switched on by requesting fewer
#       singular triplets than min(dim(x)).  This one is by far the larger
#       error source in practice -- see ?svd.HDF5Matrix.
.attach_svd_path_info <- function(res, x, used_method, exact,
                                  method_requested, k, q,
                                  rank = NULL, nev = 0L) {

    n_elem     <- prod(as.numeric(dim(x)))
    thr        <- rcpp_svd_auto_threshold()
    nev        <- as.integer(nev)
    truncated  <- identical(used_method, "blocks") && nev > 0L

    attr(res, "method")         <- used_method
    attr(res, "exact")          <- isTRUE(exact) && !truncated
    attr(res, "path")           <- used_method
    attr(res, "elements")       <- n_elem
    attr(res, "auto_threshold") <- thr
    attr(res, "blocking")       <- c(k = as.integer(k), q = as.integer(q))
    # Report the truncation that was APPLIED, not the one that was requested:
    # below the auto boundary the exact LAPACK path runs and the reduction to
    # `nev` columns happens afterwards, losslessly, in R.
    attr(res, "nev")            <- if (truncated) nev else 0L
    attr(res, "truncated")      <- truncated
    if (!is.null(rank)) attr(res, "rank") <- as.integer(rank)

    if (!isTRUE(exact) && identical(method_requested, "auto"))
        message(sprintf(
            paste0("BigDataStatMeth: %s x %s (%.3g elements) is at or above the ",
                   "exact/approximate boundary (%.3g elements), so an APPROXIMATE ",
                   "hierarchical block decomposition was used (k = %d, q = %d). ",
                   "See attr(., \"method\") and ?svd.HDF5Matrix; ",
                   "method = \"full\" forces the exact path (needs the whole ",
                   "matrix in RAM)."),
            nrow(x), ncol(x), n_elem, thr, as.integer(k), as.integer(q)))

    if (truncated)
        message(sprintf(
            paste0("BigDataStatMeth: %d of %d singular triplets requested, so each ",
                   "local block was truncated to rank %d before merging. The result ",
                   "is an APPROXIMATION -- relative errors of 1e-3 to 1e-1 have been ",
                   "measured in this regime, including in the LEADING singular value. ",
                   "attr(., \"truncated\") records it. Request the full ",
                   "min(dim(x)) = %d triplets (or a generously oversampled rank) for ",
                   "an accurate decomposition. See ?svd.HDF5Matrix."),
            nev, min(dim(x)), nev, min(dim(x))))

    res
}


# ── $svd() ──────────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "svd",
# @description
# Block-wise Singular Value Decomposition of an HDF5Matrix.
#
# Computes A = U diag(d) t(V) entirely on disk.  The singular values
# \code{d} are loaded into memory (they are always a small vector of
# length \code{min(nrow, ncol)}).  \code{u} and \code{v} are returned
# as \code{HDF5Matrix} objects pointing to the result datasets on disk.
#
# @param k          Number of local SVDs per level (default 2).
# @param q          Number of incremental levels (default 1).
# @param nev        Number of singular triplets to compute (0 = all).
# @param center     Logical. Center columns before decomposition.
# @param scale      Logical. Scale columns before decomposition.
# @param method     "auto" (default), "blocks", or "full".
# @param rankthreshold Threshold in [0, 0.1] for rank approximation.
# @param overwrite  Logical. Overwrite existing SVD results.
#   When TRUE, any open HDF5Matrix handles pointing to the output paths
#   (SVD/<dataset>/d, SVD/<dataset>/u, SVD/<dataset>/v) are automatically
#   closed and invalidated before the datasets are overwritten.
# @param threads    Integer. OpenMP threads (-1 = auto).
# @return Named list with elements:
#   \describe{
#     \item{d}{Numeric vector of singular values (non-negative, decreasing).}
#     \item{u}{HDF5Matrix of left  singular vectors (m × r).}
#     \item{v}{HDF5Matrix of right singular vectors (n × r).}
#   }
function(k             = 2,
         q             = 1,
         nev           = 0,
         center        = TRUE,
         scale         = TRUE,
         method        = "auto",
         rankthreshold = 0.0,
         overwrite     = FALSE,
         threads       = -1L) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # Report a requested thread count this system cannot honour, then carry on
    # with whatever it can give.
    threads <- .normalize_threads(threads)
    .warn_threads_capped(threads)

    # When overwrite=TRUE, close any live handles pointing to output paths
    if (isTRUE(overwrite)) {
        svd_root <- paste0("SVD/", private$dataset, "/")
        rcpp_hdf5_close_at_paths(
            private$filename,
            c(paste0(svd_root, "d"),
              paste0(svd_root, "u"),
              paste0(svd_root, "v"))
        )
    }

    # ── Close own handle — C++ will open its own handles to the same file ──
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    on.exit({
        if (is.null(private$ptr))
            private$ptr <- rcpp_hdf5dataset_open(private$filename, private$group, private$dataset)
    }, add = TRUE)

    res <- rcpp_hdf5dataset_svd(
        private$filename,
        private$group,
        private$dataset,
        k             = as.integer(k),
        q             = as.integer(q),
        nev           = as.integer(nev),
        bcenter       = isTRUE(center),
        bscale        = isTRUE(scale),
        rankthreshold = as.double(rankthreshold),
        overwrite     = isTRUE(overwrite),
        method        = as.character(method),
        threads       = as.integer(threads)
    )

    # ── Reopen after C++ has finished and released all its handles ──────────
    private$ptr <- rcpp_hdf5dataset_open( private$filename, private$group, private$dataset )

    # d is always a 1 × r matrix in HDF5 — load into a plain numeric vector
    d_tmp <- hdf5_matrix(res$file, res$path_d)
    d_vec <- as.numeric(as.matrix(d_tmp))
    close(d_tmp)

    out <- list(
        d = d_vec,
        u = hdf5_matrix(res$file, res$path_u),
        v = hdf5_matrix(res$file, res$path_v)
    )

    .attach_svd_path_info(out, self,
                          used_method      = res$used_method,
                          exact            = res$exact,
                          method_requested = as.character(method),
                          k = k, q = q, rank = length(d_vec),
                          nev = nev)

})


# ── $pca() ──────────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "pca",
# @description
# Block-wise Principal Component Analysis of an HDF5Matrix.
#
# Uses SVD internally (A = UDV') and computes all PCA statistics on disk.
# Follows the \code{prcomp()} contract as closely as possible: \code{sdev}
# and \code{cumvar} are loaded into memory (small vectors); heavy matrices
# (\code{rotation}, \code{x}, etc.) remain as \code{HDF5Matrix} objects.
#
# @param ncomponents Number of PCs to compute (0 = all).
# @param center      Logical. Center columns before PCA.
# @param scale       Logical. Scale columns before PCA.
# @param k           Number of local SVDs per level (default 2).
# @param q           Number of incremental levels (default 1).
# @param method      "auto" (default), "blocks", or "full".
# @param rankthreshold Threshold in [0, 0.1] for rank approximation.
# @param svdgroup    HDF5 group for intermediate SVD storage ("SVD/" default).
# @param overwrite   Logical. Recompute even if PCA results exist.
#   When TRUE, any open HDF5Matrix handles pointing to PCA output paths
#   (PCA/<dataset>/*) are automatically closed and invalidated.
# @param threads     Integer. OpenMP threads (-1 = auto).
# @return Object of class \code{c("HDF5PCA", "list")} with elements:
#   \describe{
#     \item{sdev}{Numeric vector. Standard deviations of the PCs (sqrt of eigenvalues).}
#     \item{rotation}{HDF5Matrix. Variable loadings — equivalent to \code{prcomp()$rotation}.}
#     \item{x}{HDF5Matrix. Individual coordinates — equivalent to \code{prcomp()$x}.}
#     \item{center}{Logical. Whether columns were centered.}
#     \item{scale}{Logical. Whether columns were scaled.}
#     \item{cumvar}{Numeric vector. Cumulative variance explained (0-100).}
#     \item{lambda}{Numeric vector. Eigenvalues.}
#     \item{var.cos2}{HDF5Matrix. Squared cosines for variables.}
#     \item{ind.cos2}{HDF5Matrix. Squared cosines for individuals.}
#     \item{ind.contrib}{HDF5Matrix. Contributions of individuals to PCs.}
#     \item{file}{Character. Path to the HDF5 file containing all results.}
#   }
#   The object also carries the exact/approximate path metadata described in
#   \code{?svd.HDF5Matrix} (attributes \code{method}, \code{exact},
#   \code{elements}, \code{auto_threshold}, \code{blocking}, \code{rank}),
#   because PCA is computed through the same SVD machinery.
function(ncomponents   = 0L,
         center        = FALSE,
         scale         = FALSE,
         k             = 2L,
         q             = 1L,
         method        = "auto",
         rankthreshold = 0.0,
         svdgroup      = "SVD/",
         overwrite     = FALSE,
         threads       = -1L) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # Report a requested thread count this system cannot honour, then carry on
    # with whatever it can give.
    threads <- .normalize_threads(threads)
    .warn_threads_capped(threads)

    # Blocking parameters, kept aside: the truncation block below reuses `k`
    # for the component count.
    k_blocking <- as.integer(k)
    q_blocking <- as.integer(q)

    # When overwrite=TRUE, close any live handles pointing to PCA output paths.
    # PCA also uses SVD internally, so close those paths too.
    if (isTRUE(overwrite)) {
        pca_root <- paste0("PCA/", private$dataset, "/")
        svd_root <- paste0(if (nzchar(svdgroup)) svdgroup else "SVD/",
                           private$dataset, "/")
        rcpp_hdf5_close_at_paths(
            private$filename,
            c(paste0(pca_root, "lambda"),
              paste0(pca_root, "variance"),
              paste0(pca_root, "cumvar"),
              paste0(pca_root, "var.coord"),
              paste0(pca_root, "var.cos2"),
              paste0(pca_root, "ind.dist"),
              paste0(pca_root, "components"),
              paste0(pca_root, "ind.coord"),
              paste0(pca_root, "ind.cos2"),
              paste0(pca_root, "ind.contrib"),
              paste0(svd_root, "d"),
              paste0(svd_root, "u"),
              paste0(svd_root, "v"))
        )
    }

    # ── Close own handle — C++ will open its own handles to the same file ──
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    on.exit({
        if (is.null(private$ptr))
            private$ptr <- rcpp_hdf5dataset_open(private$filename, private$group, private$dataset)
    }, add = TRUE)

    res <- rcpp_hdf5dataset_pca(
        private$filename,
        private$group,
        private$dataset,
        ncomponents   = as.integer(ncomponents),
        bcenter       = isTRUE(center),
        bscale        = isTRUE(scale),
        k             = as.integer(k),
        q             = as.integer(q),
        rankthreshold = as.double(rankthreshold),
        svdgroup      = as.character(svdgroup),
        overwrite     = isTRUE(overwrite),
        method        = as.character(method),
        threads       = as.integer(threads)
    )

    # ── Reopen after C++ has finished and released all its handles ──────────
    private$ptr <- rcpp_hdf5dataset_open( private$filename, private$group, private$dataset )

    # Small vectors — load into memory
    lambda_tmp <- hdf5_matrix(res$file, res$path_lambda)
    lambda     <- as.numeric(as.matrix(lambda_tmp))
    close(lambda_tmp)
    cumvar_tmp <- hdf5_matrix(res$file, res$path_cumvar)
    cumvar     <- as.numeric(as.matrix(cumvar_tmp))
    close(cumvar_tmp)
    sdev       <- sqrt(abs(lambda))   # standard deviations of PCs

    # Truncate to requested number of components
    if (ncomponents > 0L && ncomponents < length(sdev)) {
        k      <- as.integer(ncomponents)
        sdev   <- sdev[seq_len(k)]
        lambda <- lambda[seq_len(k)]
        cumvar <- cumvar[seq_len(k)]
    }

    result <- list(
        sdev        = sdev,
        rotation    = hdf5_matrix(res$file, res$path_var_coord),
        x           = hdf5_matrix(res$file, res$path_ind_coord),
        center      = isTRUE(center),
        scale       = isTRUE(scale),
        cumvar      = cumvar,
        lambda      = lambda,
        var.cos2    = hdf5_matrix(res$file, res$path_var_cos2),
        ind.cos2    = hdf5_matrix(res$file, res$path_ind_cos2),
        ind.contrib = hdf5_matrix(res$file, res$path_ind_contrib),
        file        = res$file
    )

    class(result) <- c("HDF5PCA", "list")

    # PCA runs an SVD underneath and inherits its exact/approximate switch —
    # record which path was taken (see .attach_svd_path_info above).
    .attach_svd_path_info(result, self,
                          used_method      = res$used_method,
                          exact            = res$exact,
                          method_requested = as.character(method),
                          k = k_blocking, q = q_blocking,
                          rank = length(sdev),
                          nev = ncomponents)
})
