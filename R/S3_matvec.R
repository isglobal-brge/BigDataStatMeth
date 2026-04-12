# S3_matvec.R
#
# S3 methods for sweep, diag/diag<-, diag_op, diag_scale for HDF5Matrix.
#
# Design notes:
#   - sweep():   base::sweep(x, MARGIN, STATS, FUN, check.margin, ...) has "..."
#                -> no generic override needed.
#                -> Use @exportS3Method base::sweep HDF5Matrix (same as chol).
#   - diag():    base::diag(x, nrow, ncol, names) has NO "..."
#                -> must override the generic (same pattern as qr()/eigen()).
#   - diag<-():  base::`diag<-`(x, value) already exists as generic.
#                -> define only the HDF5Matrix method with @export.
#   - diag_op():    new generic, no base R equivalent.
#   - diag_scale(): new generic, no base R equivalent.


# -- sweep.HDF5Matrix ---------------------------------------------------------

#' @export
sweep <- function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)
    UseMethod("sweep")

#' @exportS3Method
sweep.default <- function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)
    base::sweep(x, MARGIN = MARGIN, STATS = STATS, FUN = FUN,
                check.margin = check.margin, ...)



#' Broadcast a vector over an HDF5Matrix (sweep)
#'
#' @description
#' S3 method of \code{base::sweep()} for \code{HDF5Matrix} objects.
#' Broadcasts a 1-row \code{HDF5Matrix} (acting as the STATS vector)
#' across every row or column of the matrix, element-wise.
#'
#' @param x           An \code{HDF5Matrix} (the matrix to sweep over).
#' @param MARGIN      Integer. \code{2} (default, sweep over columns) or
#'   \code{1} (sweep over rows). Corresponds to \code{byrows}:
#'   \code{MARGIN = 1} means \code{byrows = TRUE}.
#' @param STATS       An \code{HDF5Matrix} with one row or one column
#'   (the vector to broadcast). Must be in the same HDF5 file as \code{x}.
#' @param FUN         Character. Operation: \code{"+"}, \code{"-"},
#'   \code{"*"} (default), \code{"/"}, \code{"pow"}.
#' @param check.margin Ignored (kept for S3 signature compatibility).
#' @param paral       Logical or NULL.
#' @param threads     Integer or NULL.
#' @param compression Integer (0-9) or NULL.
#' @param ...         Ignored.
#' @return A new \code{HDF5Matrix}.
#'
#' @examples
#' \donttest{
#' X         <- hdf5_matrix("data.h5", "data/X")   # 5000 x 800
#' col_means <- colMeans(X)                         # HDF5Matrix (1 x 800)
#' X_c       <- sweep(X, 2, col_means, "-")         # column-center X
#' }
#'
#' @export 
sweep.HDF5Matrix <- function(x,
                              MARGIN       = 2L,
                              STATS,
                              FUN          = "*",
                              check.margin = TRUE,
                              paral        = NULL,
                              threads      = NULL,
                              compression  = NULL,
                              ...) {
    if (!inherits(STATS, "HDF5Matrix"))
        stop("sweep.HDF5Matrix: STATS must be an HDF5Matrix")
    byrows <- (as.integer(MARGIN) == 1L)
    x$sweep(STATS, FUN = FUN, byrows = byrows,
            paral = paral, threads = threads, compression = compression)
}


# -- diag() -- must override because base::diag has no ... -------------------

#' Extract or construct a diagonal for HDF5Matrix
#'
#' @description
#' Overrides \code{base::diag()} to dispatch on \code{HDF5Matrix} objects.
#' For plain R matrices/vectors the call is forwarded to \code{base::diag()}.
#'
#' When \code{x} is an \code{HDF5Matrix}, extracts the diagonal as an
#' in-memory numeric vector (length = min(nrow, ncol)).
#'
#' @param x     An \code{HDF5Matrix} or any object accepted by \code{base::diag()}.
#' @param nrow  Passed to \code{base::diag()} for non-HDF5Matrix objects.
#' @param ncol  Passed to \code{base::diag()} for non-HDF5Matrix objects.
#' @param names Passed to \code{base::diag()} for non-HDF5Matrix objects.
#' @param ...   Ignored.
#' @return For \code{HDF5Matrix}: numeric vector of diagonal elements.
#'   For plain R objects: result of \code{base::diag()}.
#'
#' @seealso \code{\link[base]{diag}}
#'
#' @export
diag <- function(x, ...) UseMethod("diag")

#' @rdname diag
#' @exportS3Method
diag.default <- function(x, nrow, ncol, names = TRUE, ...) {
    if (missing(nrow) && missing(ncol))
        base::diag(x)
    else if (missing(ncol))
        base::diag(x, nrow = nrow)
    else if (missing(nrow))
        base::diag(x, ncol = ncol)
    else
        base::diag(x, nrow = nrow, ncol = ncol)
}
# diag.default <- function(x, nrow, ncol, names = TRUE, ...) {
#     if (missing(nrow) && missing(ncol))
#         base::diag(x, names = names)
#     else if (missing(ncol))
#         base::diag(x, nrow = nrow, names = names)
#     else if (missing(nrow))
#         base::diag(x, ncol = ncol, names = names)
#     else
#         base::diag(x, nrow = nrow, ncol = ncol, names = names)
# }

#' @rdname diag
#' @exportS3Method
diag.HDF5Matrix <- function(x, ...) {
    x$diag()
}

# -- diag<-() ----------------------------------------------------------------
#' @rawNamespace export("diag<-")
#' @rawNamespace S3method("diag<-",HDF5Matrix)
#' @rawNamespace S3method("diag<-",default)
NULL

# generic — dispatches to HDF5Matrix or default
#' @noRd 
`diag<-` <- function(x, value) UseMethod("diag<-")

# fallback for plain R matrices/vectors
#' @noRd 
`diag<-.default` <- function(x, value) base::`diag<-`(x, value)

# HDF5Matrix method — writes diagonal in-place via C++
#' @noRd 
`diag<-.HDF5Matrix` <- function(x, value) {
    x$diag_assign(value)
    invisible(x)
}