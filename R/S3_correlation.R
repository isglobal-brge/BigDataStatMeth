# S3_correlation.R
#
# S3 dispatch for cor() on HDF5Matrix objects.
#
# base::cor() is NOT a proper S3 generic (it does not call UseMethod internally),
# so we must mask it with a UseMethod-based generic — exactly the same pattern
# used for scale(), var(), sd(), colSums(), rowSums(), etc. in S3_aggregations.R.
#
# cor.HDF5Matrix follows the base R cor() contract as closely as possible:
#   - x       : HDF5Matrix (required)
#   - y       : HDF5Matrix or NULL (cross-correlation)
#   - use     : character — "everything" / "complete.obs" / ...
#   - method  : "pearson" (default) or "spearman"
#
# Additional HDF5-specific parameters (not in base::cor):
#   - trans_x, trans_y   : correlate rows instead of columns
#   - compute_pvalues    : also store p-values on disk
#   - block_size, threads: performance tuning
#   - result_path        : output location (NULL = auto)
#
# The result is an HDF5Matrix pointing to the correlation dataset on disk,
# with metadata attributes:  cor.method, cor.type, cor.n.vars, cor.n.obs,
# and (if computed) cor.pvalues.path.


# ---------------------------------------------------------------------------
# Generic override  (base::cor has no UseMethod → must mask)
# ---------------------------------------------------------------------------

# #' @export
# cor <- function(x, y = NULL, use = "everything", method = "pearson", ...)
#     UseMethod("cor")

#' Correlation (generic)
#'
#' @description
#' S3 generic for \code{cor()}. Dispatches to \code{\link{cor.HDF5Matrix}}
#' for \code{HDF5Matrix} objects, and to \code{stats::cor()} for all others.
#'
#' @param x A matrix or \code{HDF5Matrix} object.
#' @param y Optional second matrix.
#' @param use Character. Method for handling missing values.
#' @param method Character. Correlation method: \code{"pearson"} or \code{"spearman"}.
#' @param ... Additional arguments passed to the method.
#' @return Correlation matrix.
#' @name cor
#' @rdname cor
#' @export
cor <- function(x, y = NULL, use = "everything", method = "pearson", ...)
    UseMethod("cor")

#' @exportS3Method
cor.default <- function(x, y = NULL, use = "everything", method = "pearson", ...)
    stats::cor(x, y = y, use = use, method = method)


# ---------------------------------------------------------------------------
# S3 method for HDF5Matrix
# ---------------------------------------------------------------------------

#' Correlation matrix for HDF5Matrix objects
#'
#' Block-wise computation of Pearson or Spearman correlation, running
#' entirely on disk without loading the full matrix into RAM.
#' Supports both auto-correlation \code{cor(X)} and cross-correlation
#' \code{cor(X, Y)}.
#'
#' @param x             An \code{HDF5Matrix} object.
#' @param y             An \code{HDF5Matrix} for cross-correlation, or
#'                      \code{NULL} (default) to compute \code{cor(x, x)}.
#' @param use           Character string.  Only \code{"everything"} (default)
#'                      and \code{"complete.obs"} are currently supported.
#' @param method        \code{"pearson"} (default) or \code{"spearman"}.
#' @param trans_x       Logical.  If \code{TRUE}, correlate rows of \code{x}
#'                      instead of columns (useful for sample-sample
#'                      correlations in omics data). Default \code{FALSE}.
#' @param trans_y       Logical.  Same for \code{y}. Default \code{FALSE}.
#' @param compute_pvalues Logical. Also compute and store p-values on disk.
#'                      Default \code{TRUE}.
#' @param block_size    Integer or NULL. Block size for HDF5 reads (NULL = auto).
#' @param threads       Integer or NULL. Number of OpenMP threads (NULL = auto).
#' @param result_path   Output location:
#'   \code{NULL} (default) writes to \code{"CORR/<dataset>/correlation"} in the
#'   same file as \code{x}. A character string specifies a custom output group
#'   in the same file.  A named list \code{list(file=, group=)} writes to a
#'   different file.
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result datasets.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 6).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @param ...           Ignored (for S3 compatibility).
#'
#' @return An \code{HDF5Matrix} pointing to the correlation matrix on disk.
#'   Attributes attached to the result:
#'   \describe{
#'     \item{\code{cor.method}}{The correlation method used.}
#'     \item{\code{cor.type}}{\code{"single"} or \code{"cross"}.}
#'     \item{\code{cor.n.vars}}{Number of variables (columns/rows correlated).}
#'     \item{\code{cor.n.obs}}{Number of observations used.}
#'     \item{\code{cor.pvalues.path}}{HDF5 path to the p-values dataset
#'       (present only when \code{compute_pvalues = TRUE}).}
#'   }
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' X   <- hdf5_create_matrix(tmp, "data/X",
#'                            data = matrix(rnorm(500), 50, 10))
#'
#' # Auto-correlation: cor(X) — 10 x 10 matrix
#' C <- cor(X)
#' dim(C)
#' cat("method:", attr(C, "cor.method"), "\n")
#'
#' # Spearman
#' Cs <- cor(X, method = "spearman")
#' dim(Cs)
#'
#' # Sample-sample correlation (rows)
#' Sr <- cor(X, trans_x = TRUE)   # 50 x 50
#' dim(Sr)
#'
#' X$close(); C$close(); Cs$close(); Sr$close()
#' unlink(tmp)
#' }
#'
#' @export
cor.HDF5Matrix <- function(x,
                           y               = NULL,
                           use             = "everything",
                           method          = "pearson",
                           trans_x         = FALSE,
                           trans_y         = FALSE,
                           compute_pvalues = TRUE,
                           block_size      = NULL,
                           threads         = NULL,
                           result_path     = NULL,
                           compression     = NULL,
                           ...) {
    
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")
    
    x$cor(y               = y,
          method          = method,
          use             = use,
          trans_x         = isTRUE(trans_x),
          trans_y         = isTRUE(trans_y),
          compute_pvalues = isTRUE(compute_pvalues),
          block_size      = block_size,
          threads         = threads,
          result_path     = result_path,
          compression     = compression)
}