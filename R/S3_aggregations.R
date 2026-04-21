# S3 methods - Aggregations
# Standard R aggregation generics dispatched to HDF5Matrix R6 methods.
#
# All functions accept an optional `save_to` argument that allows the result
# to be persisted to disk without loading the full matrix into RAM, which is
# useful when the result will be reused later and recomputation is expensive.
#
# save_to formats accepted:
#   NULL            -> result returned as an R vector/scalar (default)
#   "group/dataset" -> saved in the same HDF5 file as x, returns HDF5Matrix
#   list(file  = "path/to/other.h5",
#        path  = "group/dataset")
#                   -> saved in a different HDF5 file, returns HDF5Matrix
#
# When save_to is not NULL the function:
#   1. Computes the aggregation (block-wise on disk).
#   2. Writes the result vector/scalar as a 1-D HDF5 dataset.
#   3. Returns an HDF5Matrix handle pointing to that dataset so the caller
#      can retrieve or reuse it without recomputing.


# ---------------------------------------------------------------------------
# Generic overrides for base R functions that do NOT dispatch S3 properly
# ---------------------------------------------------------------------------
#
# base::colSums / rowSums / colMeans / rowMeans check !is.array(x) BEFORE
# calling UseMethod(), so they never reach *.HDF5Matrix methods.
# Masking them here with proper UseMethod generics fixes dispatch for all
# classes while delegating to base for normal matrix/array/data.frame objects.
#
# base::var and base::sd are not S3 generics at all (no UseMethod call),
# so they must also be overridden here.
#
# base::mean IS a proper S3 generic → no override needed.
# base::sum/min/max dispatch via the Summary group generic → no override needed.
#
# base::scale is defined as function(x, center, scale) with NO ..., so any
# extra argument (e.g. result_path) causes "unused argument" before dispatch.

#' @export
scale <- function(x, center = TRUE, scale = TRUE, ...)
    UseMethod("scale")

#' @exportS3Method
scale.default <- function(x, center = TRUE, scale = TRUE, ...)
    base::scale.default(x, center = center, scale = scale)
# scale.default <- function(x, center = TRUE, scale = TRUE, ...)
#     base::scale(x, center = center, scale = scale)



#' Column and row sums for HDF5Matrix
#'
#' Block-wise computation of column and row sums without loading the full
#' matrix into RAM.  Results can optionally be persisted to an HDF5 file.
#'
#' @param x        An \code{HDF5Matrix} object.
#' @param na.rm    Ignored (included for compatibility with the base generic).
#' @param dims     Ignored (included for compatibility with the base generic).
#' @param paral    Logical or NULL. Enable OpenMP parallelisation.
#' @param wsize    Integer or NULL. Block size for HDF5 reads (NULL = auto).
#' @param threads  Integer or NULL. Number of OpenMP threads (NULL = auto).
#' @param save_to  Where to save the result (see Details).  NULL returns a
#'                 plain R vector; a character string \code{"group/dataset"}
#'                 saves in the same file as \code{x}; a named list
#'                 \code{list(file = "f.h5", path = "group/dataset")} saves in
#'                 a different file.  In both non-NULL cases an
#'                 \code{HDF5Matrix} handle is returned.
#' @param overwrite Logical.  Overwrite an existing dataset at \code{save_to}?
#'                 (Currently always TRUE when \code{save_to} is not NULL.)
#' @param ...      Ignored.
#'
#' @return A numeric vector (when \code{save_to = NULL}) or an
#'   \code{HDF5Matrix} handle to the persisted result.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' X <- hdf5_create_matrix(tmp, "data/M", data = matrix(rnorm(200), 20, 10))
#' cs <- colSums(X)         # R vector
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @export
colSums <- function(x, na.rm = FALSE, dims = 1L, ...)
    UseMethod("colSums")

#' @exportS3Method
colSums.default <- function(x, na.rm = FALSE, dims = 1L, ...)
    base::colSums(x, na.rm = na.rm, dims = dims)

#' @rdname colSums
#' @export
rowSums <- function(x, na.rm = FALSE, dims = 1L, ...)
    UseMethod("rowSums")

#' @exportS3Method
rowSums.default <- function(x, na.rm = FALSE, dims = 1L, ...)
    base::rowSums(x, na.rm = na.rm, dims = dims)





#' Column and row means for HDF5Matrix
#'
#' Block-wise computation of column and row means without loading the full
#' matrix into RAM.
#'
#' @inheritParams colSums.HDF5Matrix
#' @return A numeric vector (when \code{save_to = NULL}) or an
#'   \code{HDF5Matrix} handle to the persisted result.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' X <- hdf5_create_matrix(tmp, "data/M", data = matrix(rnorm(200), 20, 10))
#' cm <- colMeans(X)
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @export
colMeans <- function(x, na.rm = FALSE, dims = 1L, ...)
    UseMethod("colMeans")

#' @exportS3Method
colMeans.default <- function(x, na.rm = FALSE, dims = 1L, ...)
    base::colMeans(x, na.rm = na.rm, dims = dims)

#' @rdname colMeans
#' @export
rowMeans <- function(x, na.rm = FALSE, dims = 1L, ...)
    UseMethod("rowMeans")

#' @exportS3Method
rowMeans.default <- function(x, na.rm = FALSE, dims = 1L, ...)
    base::rowMeans(x, na.rm = na.rm, dims = dims)

#' @rdname var.HDF5Matrix
#' @export
var <- function(x, y = NULL, na.rm = FALSE, use, ...)
    UseMethod("var")

#' @exportS3Method
var.default <- function(x, y = NULL, na.rm = FALSE, use, ...)
    stats::var(x, y = y, na.rm = na.rm,
               use = if (missing(use)) "everything" else use)

#' @rdname sd.HDF5Matrix
#' @export
sd <- function(x, na.rm = FALSE, ...)
    UseMethod("sd")

#' @exportS3Method
sd.default <- function(x, na.rm = FALSE, ...)
    stats::sd(x, na.rm = na.rm)


# ---------------------------------------------------------------------------
# Internal: persist a numeric vector/scalar result to HDF5
# ---------------------------------------------------------------------------

#' @noRd
.agg_save_result <- function(result, save_to, src_x) {
    # Resolve target file and full dataset path
    if (is.character(save_to)) {
        info     <- src_x$info()
        out_file <- info$filename
        out_path <- save_to
    } else if (is.list(save_to)) {
        if (is.null(save_to$file) || is.null(save_to$path))
            stop("save_to list must have 'file' and 'path' elements")
        out_file <- save_to$file
        out_path <- save_to$path
    } else {
        stop("save_to must be NULL, a character string, or a list(file, path)")
    }

    # For a vector treat it as a 1-column matrix so hdf5_create_matrix
    # receives a proper 2-D object.
    mat <- if (is.null(dim(result))) matrix(result, ncol = 1L) else result

    # hdf5_create_matrix(filename, dataset_full_path, data=, overwrite=)
    hdf5_create_matrix(out_file, out_path, data = mat, overwrite = TRUE)

    # Return an open HDF5Matrix handle pointing to the saved dataset
    hdf5_matrix(out_file, out_path)
}


# ---------------------------------------------------------------------------
# Helper: apply names from HDF5Matrix to a result vector when available
# ---------------------------------------------------------------------------

#' @noRd
.agg_apply_colnames <- function(vec, x) {
    cn <- tryCatch(base::colnames(x), error = function(e) NULL)
    if (!is.null(cn) && length(cn) == length(vec)) names(vec) <- cn
    vec
}

#' @noRd
.agg_apply_rownames <- function(vec, x) {
    rn <- tryCatch(base::rownames(x), error = function(e) NULL)
    if (!is.null(rn) && length(rn) == length(vec)) names(vec) <- rn
    vec
}


# ---------------------------------------------------------------------------
# colSums / rowSums
# ---------------------------------------------------------------------------

# Column and row sums for HDF5Matrix
#'
# Block-wise computation of column and row sums without loading the full
# matrix into RAM.  Results can optionally be persisted to an HDF5 file.
#'
# @param x        An \code{HDF5Matrix} object.
# @param na.rm    Ignored (included for compatibility with the base generic).
# @param dims     Ignored (included for compatibility with the base generic).
# @param paral    Logical or NULL. Enable OpenMP parallelisation.
# @param wsize    Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads  Integer or NULL. Number of OpenMP threads (NULL = auto).
# @param save_to  Where to save the result (see Details).  NULL returns a
#                 plain R vector; a character string \code{"group/dataset"}
#                 saves in the same file as \code{x}; a named list
#                 \code{list(file = "f.h5", path = "group/dataset")} saves in
#                 a different file.  In both non-NULL cases an
#                 \code{HDF5Matrix} handle is returned.
# @param overwrite Logical.  Overwrite an existing dataset at \code{save_to}?
#                 (Currently always TRUE when \code{save_to} is not NULL.)
# @param ...      Ignored.
#'
# @return A numeric vector (when \code{save_to = NULL}) or an
#   \code{HDF5Matrix} handle to the persisted result.
#'
# @examples
# \donttest{
# tmp <- tempfile(fileext = ".h5")
# X <- hdf5_create_matrix(tmp, "data/M",data = matrix(rnorm(200), 20, 10))
# cs <- colSums(X)                                     # R vector
# hdf5_close_all()
# unlink(tmp)
# }
#'
# @export
#' @rdname colSums
#' @export
colSums.HDF5Matrix <- function(x, na.rm = FALSE, dims = 1,
                                paral   = NULL,
                                wsize   = NULL,
                                threads = NULL,
                                save_to = NULL,
                                overwrite = TRUE,
                                ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_colnames(x$colSums(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}

#' @rdname colSums
#' @export
rowSums.HDF5Matrix <- function(x, na.rm = FALSE, dims = 1,
                                paral   = NULL,
                                wsize   = NULL,
                                threads = NULL,
                                save_to = NULL,
                                overwrite = TRUE,
                                ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_rownames(x$rowSums(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}


# ---------------------------------------------------------------------------
# colMeans / rowMeans
# ---------------------------------------------------------------------------

# Column and row means for HDF5Matrix
#
# Block-wise computation of column and row means without loading the full
# matrix into RAM.
#
# @inheritParams colSums.HDF5Matrix
# @return A numeric vector (when \code{save_to = NULL}) or an
#   \code{HDF5Matrix} handle to the persisted result.
#
# @examples
# \donttest{
# tmp <- tempfile(fileext = ".h5")
# X <- hdf5_create_matrix(tmp, "data/M",
#                          data = matrix(rnorm(200), 20, 10))
# cm <- colMeans(X)
# hdf5_close_all()
# unlink(tmp)
# }
#
# @export
#' @rdname colMeans
#' @export
colMeans.HDF5Matrix <- function(x, na.rm = FALSE, dims = 1,
                                 paral   = NULL,
                                 wsize   = NULL,
                                 threads = NULL,
                                 save_to = NULL,
                                 overwrite = TRUE,
                                 ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_colnames(x$colMeans(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}

#' @rdname colMeans
#' @export
rowMeans.HDF5Matrix <- function(x, na.rm = FALSE, dims = 1,
                                 paral   = NULL,
                                 wsize   = NULL,
                                 threads = NULL,
                                 save_to = NULL,
                                 overwrite = TRUE,
                                 ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_rownames(x$rowMeans(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}


# ---------------------------------------------------------------------------
# colMins / rowMins
# ---------------------------------------------------------------------------

#' Column and row minimums for HDF5Matrix
#'
#' Block-wise computation of column and row minimums.
#'
#' @inheritParams colSums.HDF5Matrix
#' @return A numeric vector (when \code{save_to = NULL}) or an
#'   \code{HDF5Matrix} handle to the persisted result.
#'
#' @export
colMins <- function(x, ...) UseMethod("colMins")

#' @rdname colMins
#' @export
colMins.HDF5Matrix <- function(x,
                                paral   = NULL,
                                wsize   = NULL,
                                threads = NULL,
                                save_to = NULL,
                                overwrite = TRUE,
                                ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_colnames(x$colMins(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}

#' @rdname colMins
#' @export
rowMins <- function(x, ...) UseMethod("rowMins")

#' @rdname colMins
#' @export
rowMins.HDF5Matrix <- function(x,
                                paral   = NULL,
                                wsize   = NULL,
                                threads = NULL,
                                save_to = NULL,
                                overwrite = TRUE,
                                ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_rownames(x$rowMins(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}


# ---------------------------------------------------------------------------
# colMaxs / rowMaxs
# ---------------------------------------------------------------------------

#' Column and row maximums for HDF5Matrix
#'
#' Block-wise computation of column and row maximums.
#'
#' @inheritParams colSums.HDF5Matrix
#' @return A numeric vector (when \code{save_to = NULL}) or an
#'   \code{HDF5Matrix} handle to the persisted result.
#'
#' @export
colMaxs <- function(x, ...) UseMethod("colMaxs")

#' @rdname colMaxs
#' @export
colMaxs.HDF5Matrix <- function(x,
                                paral   = NULL,
                                wsize   = NULL,
                                threads = NULL,
                                save_to = NULL,
                                overwrite = TRUE,
                                ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_colnames(x$colMaxs(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}

#' @rdname colMaxs
#' @export
rowMaxs <- function(x, ...) UseMethod("rowMaxs")

#' @rdname colMaxs
#' @export
rowMaxs.HDF5Matrix <- function(x,
                                paral   = NULL,
                                wsize   = NULL,
                                threads = NULL,
                                save_to = NULL,
                                overwrite = TRUE,
                                ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_rownames(x$rowMaxs(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}


# ---------------------------------------------------------------------------
# colVars / rowVars
# ---------------------------------------------------------------------------

#' Column and row variances for HDF5Matrix
#'
#' Block-wise computation of column and row variances (Bessel's correction,
#' n-1).  Returns NaN for columns/rows with fewer than 2 observations,
#' matching base R behaviour.
#'
#' @inheritParams colSums.HDF5Matrix
#' @return A numeric vector (when \code{save_to = NULL}) or an
#'   \code{HDF5Matrix} handle to the persisted result.
#'
#' @export
colVars <- function(x, ...) UseMethod("colVars")

#' @rdname colVars
#' @export
colVars.HDF5Matrix <- function(x,
                                paral   = NULL,
                                wsize   = NULL,
                                threads = NULL,
                                save_to = NULL,
                                overwrite = TRUE,
                                ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_colnames(x$colVars(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}

#' @rdname colVars
#' @export
rowVars <- function(x, ...) UseMethod("rowVars")

#' @rdname colVars
#' @export
rowVars.HDF5Matrix <- function(x,
                                paral   = NULL,
                                wsize   = NULL,
                                threads = NULL,
                                save_to = NULL,
                                overwrite = TRUE,
                                ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_rownames(x$rowVars(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}


# ---------------------------------------------------------------------------
# colSds / rowSds
# ---------------------------------------------------------------------------

#' Column and row standard deviations for HDF5Matrix
#'
#' Block-wise computation of column and row standard deviations
#' (Bessel's correction, n-1).
#'
#' @inheritParams colSums.HDF5Matrix
#' @return A numeric vector (when \code{save_to = NULL}) or an
#'   \code{HDF5Matrix} handle to the persisted result.
#'
#' @export
colSds <- function(x, ...) UseMethod("colSds")

#' @rdname colSds
#' @export
colSds.HDF5Matrix <- function(x,
                               paral   = NULL,
                               wsize   = NULL,
                               threads = NULL,
                               save_to = NULL,
                               overwrite = TRUE,
                               ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_colnames(x$colSds(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}

#' @rdname colSds
#' @export
rowSds <- function(x, ...) UseMethod("rowSds")

#' @rdname colSds
#' @export
rowSds.HDF5Matrix <- function(x,
                               paral   = NULL,
                               wsize   = NULL,
                               threads = NULL,
                               save_to = NULL,
                               overwrite = TRUE,
                               ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- .agg_apply_rownames(x$rowSds(paral, wsize, threads), x)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}


# ---------------------------------------------------------------------------
# Scalar: mean, sum, min, max, var, sd  (whole matrix as flat vector)
#
# R dispatches sum(), min(), max() via the Summary group generic and
# mean() via its own generic.  var() and sd() are not group generics;
# we register dedicated S3 methods.
# ---------------------------------------------------------------------------

#' Summary statistics for HDF5Matrix
#'
#' Scalar aggregations over all elements of an HDF5 matrix, computed
#' block-wise without loading the full data into RAM.
#'
#' @param x        An \code{HDF5Matrix} object.
#' @param ...      For \code{Summary} methods: additional objects (currently
#'                 not supported; raises an error if supplied).
#' @param na.rm    Ignored (included for generic compatibility).
#' @param paral    Logical or NULL. Enable OpenMP parallelisation.
#' @param wsize    Integer or NULL. Block size (NULL = auto).
#' @param threads  Integer or NULL. Number of threads (NULL = auto).
#' @param save_to  Optional persistence target (same format as
#'                 \code{\link{colSums.HDF5Matrix}}).
#' @param overwrite Logical.  Overwrite existing dataset when saving.
#'
#' @return A scalar numeric (when \code{save_to = NULL}) or an
#'   \code{HDF5Matrix} pointing to a 1×1 persisted dataset.
#'
#' @name HDF5Matrix-scalar-aggregations
NULL


# mean.HDF5Matrix -----------------------------------------------------------
#' @rdname HDF5Matrix-scalar-aggregations
#' @export
mean.HDF5Matrix <- function(x, na.rm = FALSE,
                             paral   = NULL,
                             wsize   = NULL,
                             threads = NULL,
                             save_to = NULL,
                             overwrite = TRUE,
                             ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- x$mean(paral, wsize, threads)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}


# Summary.HDF5Matrix --------------------------------------------------------
# Handles sum(), min(), max() via the Summary group generic.
#' @rdname HDF5Matrix-scalar-aggregations
#' @export
Summary.HDF5Matrix <- function(..., na.rm = FALSE) {
    # .Generic is set by R's dispatch mechanism
    args <- list(...)

    # Only single-object calls are supported; reject multi-object use
    # (e.g. sum(A, B) where both are HDF5Matrix) to avoid confusion.
    hdf5_args <- Filter(function(a) inherits(a, "HDF5Matrix"), args)
    if (length(hdf5_args) != 1L || length(args) != 1L)
        stop("Summary methods for HDF5Matrix support exactly one HDF5Matrix argument. ",
             "Use as.matrix() first for multi-object calls.")

    x <- hdf5_args[[1L]]
    if (!x$is_valid()) stop("Dataset is closed or invalid")

    switch(.Generic,
        "sum" = x$sum(),
        "min" = x$min(),
        "max" = x$max(),
        stop("Summary generic '", .Generic, "' is not implemented for HDF5Matrix")
    )
}


# var.HDF5Matrix ------------------------------------------------------------
#' Variance of all elements of an HDF5Matrix
#'
#' Equivalent to \code{var(as.vector(X))} — treats all matrix elements as a
#' single sample and uses Bessel's correction (N-1).
#' 
#' @param y   Ignored. Present for compatibility with \code{stats::var}.
#' @param use Ignored. Present for compatibility with \code{stats::var}.
#'
#' @inheritParams mean.HDF5Matrix
#' @return Scalar numeric or an \code{HDF5Matrix} when \code{save_to} is set.
#' @export
var.HDF5Matrix <- function(x,
                           y = NULL, 
                           na.rm = FALSE,
                           use,
                            paral   = NULL,
                            wsize   = NULL,
                            threads = NULL,
                            save_to = NULL,
                            overwrite = TRUE,
                            ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- x$var(paral, wsize, threads)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}


# sd.HDF5Matrix -------------------------------------------------------------
#' Standard deviation of all elements of an HDF5Matrix
#'
#' Equivalent to \code{sd(as.vector(X))} — uses Bessel's correction (N-1).
#'
#' @inheritParams mean.HDF5Matrix
#' @return Scalar numeric or an \code{HDF5Matrix} when \code{save_to} is set.
#' @export
sd.HDF5Matrix <- function(x,
                          na.rm = FALSE,
                          paral   = NULL,
                          wsize   = NULL,
                          threads = NULL,
                          save_to = NULL,
                          overwrite = TRUE,
                          ...) {
    if (!x$is_valid()) stop("Dataset is closed or invalid")
    result <- x$sd(paral, wsize, threads)
    if (is.null(save_to)) return(result)
    .agg_save_result(result, save_to, x)
}
