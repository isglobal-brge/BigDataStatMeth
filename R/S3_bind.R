# S3_bind.R
#
# S3 methods for cbind() and rbind() on HDF5Matrix objects (Phase 10).
#
# cbind and rbind are internal generics in R, so S3 dispatch works automatically
# when cbind.HDF5Matrix / rbind.HDF5Matrix are exported.
#
# Both methods accept one or more HDF5Matrix arguments plus standard R
# objects (plain matrices are written to a temp HDF5 dataset first).
# At least the first argument must be an HDF5Matrix.


# ── cbind() ─────────────────────────────────────────────────────────────────

#' Column-bind HDF5Matrix objects
#'
#' @description
#' Binds two or more \code{HDF5Matrix} objects by columns (appending columns
#' to the right). All matrices must have the same number of rows. The
#' operation is performed block-wise on disk.
#'
#' @param ...          One or more \code{HDF5Matrix} objects (all with the
#'   same number of rows). Plain R matrices are also accepted and will be
#'   written to a temporary HDF5 dataset automatically.
#' @param deparse.level Ignored (for S3 compatibility with base::cbind).
#' @param out_file     Output HDF5 file. \code{NULL} = same file as first argument.
#' @param out_group    Output group.   \code{NULL} = \code{"BIND"}.
#' @param out_dataset  Output dataset. \code{NULL} = auto-generated name.
#' @param block_rows   Integer. Rows per I/O block (default 1000).
#' @param overwrite    Logical. Overwrite existing output. Default \code{FALSE}.
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result datasets.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 6).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @return \code{HDF5Matrix} pointing to the combined dataset.
#'
#' @examples
#' \donttest{
#' A <- hdf5_matrix("data.h5", "grp/A")
#' B <- hdf5_matrix("data.h5", "grp/B")
#' C <- cbind(A, B)          # columns of A followed by columns of B
#' dim(C)                    # nrow(A) x (ncol(A) + ncol(B))
#' }
#'
#' @export
cbind.HDF5Matrix <- function(...,
                             deparse.level = 1,
                             out_file    = NULL,
                             out_group   = NULL,
                             out_dataset = NULL,
                             block_rows  = 1000L,
                             overwrite   = FALSE,
                             compression = NULL) {
    
    args <- list(...)
    if (length(args) == 0L) stop("cbind.HDF5Matrix: no arguments supplied")
    
    for (i in seq_along(args)) {
        if (!inherits(args[[i]], "HDF5Matrix"))
            stop("cbind.HDF5Matrix: argument ", i,
                 " is not an HDF5Matrix. Convert with hdf5_matrix() first.")
    }
    
    if (length(args) == 1L) return(args[[1L]])
    
    result <- args[[1L]]$cbind(
        args[[2L]],
        out_file    = out_file,
        out_group   = out_group,
        out_dataset = out_dataset,
        block_rows  = block_rows,
        overwrite   = overwrite,
        compression = compression
    )
    
    if (length(args) > 2L) {
        for (i in seq(3L, length(args))) {
            result <- result$cbind(args[[i]],
                                   out_file    = out_file,
                                   out_group   = out_group,
                                   block_rows  = block_rows,
                                   overwrite   = TRUE,
                                   compression = compression)
        }
    }
    
    result
}


# ── rbind() ─────────────────────────────────────────────────────────────────

#' Row-bind HDF5Matrix objects
#'
#' @description
#' Binds two or more \code{HDF5Matrix} objects by rows (appending rows
#' below). All matrices must have the same number of columns. The
#' operation is performed block-wise on disk.
#'
#' @param ...          One or more \code{HDF5Matrix} objects (all with the
#'   same number of columns).
#' @param deparse.level Ignored (for S3 compatibility with base::rbind).
#' @param out_file     Output HDF5 file. \code{NULL} = same file as first argument.
#' @param out_group    Output group.   \code{NULL} = \code{"BIND"}.
#' @param out_dataset  Output dataset. \code{NULL} = auto-generated name.
#' @param block_rows   Integer. Rows per I/O block (default 1000).
#' @param overwrite    Logical. Overwrite existing output. Default \code{FALSE}.
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result datasets.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 6).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @return \code{HDF5Matrix} pointing to the combined dataset.
#'
#' @examples
#' \donttest{
#' A <- hdf5_matrix("data.h5", "grp/A")
#' B <- hdf5_matrix("data.h5", "grp/B")
#' C <- rbind(A, B)          # rows of A followed by rows of B
#' dim(C)                    # (nrow(A) + nrow(B)) x ncol(A)
#' }
#'
#' @export
rbind.HDF5Matrix <- function(...,
                              deparse.level = 1,
                              out_file    = NULL,
                              out_group   = NULL,
                              out_dataset = NULL,
                              block_rows  = 1000L,
                              overwrite   = FALSE,
                              compression = NULL) {

    args <- list(...)
    if (length(args) == 0L) stop("rbind.HDF5Matrix: no arguments supplied")

    for (i in seq_along(args)) {
        if (!inherits(args[[i]], "HDF5Matrix"))
            stop("rbind.HDF5Matrix: argument ", i,
                 " is not an HDF5Matrix. Convert with hdf5_matrix() first.")
    }

    if (length(args) == 1L) return(args[[1L]])

    result <- args[[1L]]$rbind(
        args[[2L]],
        out_file    = out_file,
        out_group   = out_group,
        out_dataset = out_dataset,
        block_rows  = block_rows,
        overwrite   = overwrite,
        compression = compression
    )

    if (length(args) > 2L) {
        for (i in seq(3L, length(args))) {
            result <- result$rbind(args[[i]],
                                   out_file   = out_file,
                                   out_group  = out_group,
                                   block_rows = block_rows,
                                   overwrite  = TRUE,
                                   compression = compression)
        }
    }

    result
}
