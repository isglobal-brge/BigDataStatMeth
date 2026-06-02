# S3_normalize.R
#
# S3 dispatch for scale() on HDF5Matrix objects.
#
# base::scale() is a proper S3 generic (UseMethod("scale")), so dispatch
# works without overriding it — no masking needed.
#
# scale.HDF5Matrix follows the base R scale() contract:
#   - center = TRUE  or a numeric vector → center each column (or row)
#   - scale  = TRUE  or a numeric vector → scale each column (or row)
#   - The returned HDF5Matrix carries 'scaled:center' and 'scaled:scale'
#     attributes, exactly as base::scale() does for ordinary matrices.
#
# Passing pre-computed center/scale vectors is NOT supported in this version
# (the C++ backend computes them internally).  If center or scale is a
# numeric vector a warning is issued and the boolean interpretation is used.


#' Scale / normalize an HDF5Matrix
#'
#' Block-wise centering and scaling equivalent to base R \code{\link{scale}()}.
#' The computation runs entirely on disk — the full matrix is never loaded into
#' RAM.
#'
#' @param x           An \code{HDF5Matrix} object.
#' @param center      Logical (or numeric vector, see Details).  If \code{TRUE}
#'                    (default) subtract column means before scaling.
#' @param scale       Logical (or numeric vector, see Details).  If \code{TRUE}
#'                    (default) divide by column standard deviations.
#' @param byrows      Logical.  If \code{TRUE} normalize row-wise instead of
#'                    column-wise. Default \code{FALSE}.
#' @param wsize       Integer or NULL.  Block size for HDF5 reads (NULL = auto).
#' @param result_path Output location.  \code{NULL} (default) writes to
#'                    \code{"NORMALIZED/<group>/<dataset>"} in the same file.
#'                    A character string writes to that path in the same file.
#'                    A named list \code{list(file=, path=)} writes to a
#'                    different file.
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result datasets.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 6).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @param paral   Logical or NULL. Enable OpenMP parallelism. \code{TRUE}
#'   forces block-wise streaming (PATH 2) regardless of matrix size, so the
#'   thread count is respected. \code{NULL} or \code{FALSE} uses the preload
#'   path (PATH 1) when the matrix fits in RAM. Overrides the global option
#'   set by \code{\link{hdf5matrix_options}}. Default \code{NULL}.
#' @param threads Integer or NULL. Number of OpenMP threads when
#'   \code{paral = TRUE}. Always capped by \code{OMP_THREAD_LIMIT}
#'   (CRAN compliance). \code{NULL} uses the system default. Overrides
#'   the global option set by \code{\link{hdf5matrix_options}}.
#' @param ...         Ignored (for S3 compatibility).
#'
#' @details
#' Passing a pre-computed numeric vector as \code{center} or \code{scale} is
#' not currently supported.  If a vector is supplied it is coerced to a logical
#' (\code{TRUE} if \code{length(x) > 0}) and a warning is issued.
#'
#' The returned \code{HDF5Matrix} carries \code{scaled:center} and
#' \code{scaled:scale} attributes (numeric vectors), mirroring the behavior of
#' \code{base::scale()}.
#' 
#' **Performance settings:**
#'
#' Parallelization and thread count can be set globally via
#' \code{\link{hdf5matrix_options}} or passed explicitly via \code{paral}
#' and \code{threads}. Explicit parameters take priority over global options.
#' \preformatted{
#' # Global configuration
#' hdf5matrix_options(paral = TRUE, threads = 4)
#' Xs <- scale(X)           # uses 4 threads
#'
#' # Explicit per-call override
#' Xs <- scale(X, paral = TRUE, threads = 8)
#' }
#'
#' @return An \code{HDF5Matrix} pointing to the normalized dataset on disk.
#'
#' @seealso
#' \code{\link{hdf5matrix_options}} for global performance settings.
#' 
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' X   <- hdf5_create_matrix(tmp, "data/M",
#'                            data = matrix(rnorm(500), 50, 10))
#' Xs  <- scale(X)                         # center=TRUE, scale=TRUE by cols
#' cat("scaled:center[1]:", attr(Xs, "scaled:center")[1], "\n")
#' X$close(); Xs$close(); unlink(tmp)
#' }
#'
#' @name scale
#' @rdname scale
#' 
#' @export
scale.HDF5Matrix <- function(x,
                             center      = TRUE,
                             scale       = TRUE,
                             byrows      = FALSE,
                             wsize       = NULL,
                             result_path = NULL,
                             compression = NULL,
                             paral       = NULL,
                             threads     = NULL,
                              ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")

    # Coerce numeric vectors to logical with a warning
    if (is.numeric(center)) {
        warning("Passing a numeric vector as 'center' is not supported for ",
                "HDF5Matrix. Using center = TRUE.")
        center <- TRUE
    }
    if (is.numeric(scale)) {
        warning("Passing a numeric vector as 'scale' is not supported for ",
                "HDF5Matrix. Using scale = TRUE.")
        scale <- TRUE
    }

    x$normalize(center      = isTRUE(center),
                scale       = isTRUE(scale),
                byrows      = isTRUE(byrows),
                wsize       = wsize,
                result_path = result_path,
                compression = compression,
                paral       = paral,
                threads     = threads)
}
