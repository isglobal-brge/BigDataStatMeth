# S3_bind.R
#
# S3 methods for cbind() and rbind() on HDF5Matrix objects (Phase 10).
#
# cbind and rbind are internal generics in R, so S3 dispatch works automatically
# when cbind.HDF5Matrix / rbind.HDF5Matrix are exported.
#
# Both methods require every argument to be an HDF5Matrix. Plain R matrices
# are not auto-converted; convert them first with hdf5_matrix() /
# hdf5_create_matrix(). Dimensions must be conformable (same number of rows
# for cbind, same number of columns for rbind).


# Internal helper: 6-char hex hash of input dataset names.
# Deterministic: same inputs -> same hash; different inputs -> different hash.
.bind_uid <- function(n = 8L)
{
    paste(sample(c(0:9, letters[1:6]), n, replace = TRUE), collapse = "")
}

# ── cbind() ─────────────────────────────────────────────────────────────────

#' Column-bind HDF5Matrix objects
#'
#' @description
#' Binds two or more \code{HDF5Matrix} objects by columns (appending columns
#' to the right). All matrices must have the same number of rows. The
#' operation is performed block-wise on disk.
#'
#' @param ...          One or more \code{HDF5Matrix} objects (all with the
#'   same number of rows). Every argument must be an \code{HDF5Matrix};
#'   plain R matrices are not accepted -- convert them first with
#'   \code{\link{hdf5_matrix}} or \code{\link{hdf5_create_matrix}}.
#' @param deparse.level Ignored (for S3 compatibility with base::cbind).
#' @param out_file     Output HDF5 file. \code{NULL} = same file as first argument.
#' @param out_group    Output group.   \code{NULL} = \code{"BIND"}.
#' @param out_dataset  Output dataset name. \code{NULL} = auto-generated:
#'   for two inputs the name is \code{"A_cbind_B"}; for three or more inputs
#'   it is \code{"cbind_N"} where \code{N} is the number of inputs, to
#'   prevent unbounded name growth.
#' @param block_rows   Integer. Rows per I/O block (default 1000).
#' @param overwrite    Logical. Overwrite existing output. Default \code{FALSE}.
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result datasets.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 1).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @return \code{HDF5Matrix} pointing to the combined dataset.
#'
#' @examples
#' \donttest{
#' 
#' fn <- tempfile(fileext = ".h5")
#' 
#' A  <- hdf5_create_matrix(fn, "grp/A", data = matrix(rnorm(100), 10, 10))
#' B  <- hdf5_create_matrix(fn, "grp/B", data = matrix(rnorm(100), 10, 10))
#' 
#' A <- hdf5_matrix(fn, "grp/A")
#' B <- hdf5_matrix(fn, "grp/B")
#' C <- cbind(A, B)          # columns of A followed by columns of B
#' dim(C)                    # nrow(A) x (ncol(A) + ncol(B))
#' 
#' hdf5_close_all()
#' unlink(fn)
#' 
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
    
    # Same strategy as rbind: compact step names to prevent quadratic growth.
    n_args   <- length(args)
    
    # final_ds <- if (!is.null(out_dataset)) out_dataset
    # else if (n_args > 2L)      paste0("cbind_", n_args)
    # else                       NULL

    # final_ds <- if (!is.null(out_dataset)) out_dataset
    # else if (n_args > 2L)      paste0("cbind_", n_args, "_", .bind_hash(args))
    # else                       NULL
    #     
    # first_ds <- if (n_args > 2L) paste0("cbind_", n_args, "_s1") else final_ds
    
    uid      <- if (n_args > 2L) .bind_uid() else NULL
    final_ds <- if (!is.null(out_dataset)) out_dataset
    else if (n_args > 2L)      paste0("cbind_", n_args, "_", uid)
    else                       NULL
    first_ds <- if (n_args > 2L) paste0("cbind_", n_args, "_", uid, "_s1") else final_ds
    
    result <- args[[1L]]$cbind(
        args[[2L]],
        out_file    = out_file,
        out_group   = out_group,
        out_dataset = first_ds,
        block_rows  = block_rows,
        overwrite   = overwrite,
        compression = compression
    )
    
    if (n_args > 2L) {
        for (i in seq(3L, length(args))) {
            step_ds <- if (i == n_args) final_ds
            else             paste0("cbind_", n_args, "_s", i - 1L)
            result <- result$cbind(args[[i]],
                                   out_file    = out_file,
                                   out_group   = out_group,
                                   out_dataset = step_ds,
                                   block_rows  = block_rows,
                                   overwrite   = TRUE,
                                   compression = compression)
        }
    }
    
    result
    
    
    # result <- args[[1L]]$cbind(
    #     args[[2L]],
    #     out_file    = out_file,
    #     out_group   = out_group,
    #     out_dataset = out_dataset,
    #     block_rows  = block_rows,
    #     overwrite   = overwrite,
    #     compression = compression
    # )
    # 
    # if (length(args) > 2L) {
    #     for (i in seq(3L, length(args))) {
    #         result <- result$cbind(args[[i]],
    #                                out_file    = out_file,
    #                                out_group   = out_group,
    #                                block_rows  = block_rows,
    #                                overwrite   = TRUE,
    #                                compression = compression)
    #     }
    # }
    # 
    # result
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
#' @param out_dataset  Output dataset name. \code{NULL} = auto-generated:
#'   for two inputs the name is \code{"A_rbind_B"}; for three or more inputs
#'   it is \code{"rbind_N"} where \code{N} is the number of inputs, to
#'   prevent unbounded name growth.
#' @param block_rows   Integer. Rows per I/O block (default 1000).
#' @param overwrite    Logical. Overwrite existing output. Default \code{FALSE}.
#' @param compression Integer (0-9) or NULL. gzip compression level for the
#'   result datasets.  NULL uses the global option set by
#'   \code{\link{hdf5matrix_options}} (default 1).  Use \code{0} to disable
#'   compression (faster for benchmarks).
#' @return \code{HDF5Matrix} pointing to the combined dataset.
#'
#' @examples
#' \donttest{
#' 
#' fn <- tempfile(fileext = ".h5")
#' 
#' A  <- hdf5_create_matrix(fn, "grp/A", data = matrix(rnorm(100), 10, 10))
#' B  <- hdf5_create_matrix(fn, "grp/B", data = matrix(rnorm(100), 10, 10))
#' 
#' A <- hdf5_matrix(fn, "grp/A")
#' B <- hdf5_matrix(fn, "grp/B")
#' C <- rbind(A, B)          # rows of A followed by rows of B
#' dim(C)                    # (nrow(A) + nrow(B)) x ncol(A)
#' 
#' hdf5_close_all()
#' unlink(fn)
#' 
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
    
    # When binding N > 2 arguments, the auto-generated name would grow
    # quadratically because each step prepends the accumulated result name.
    # Instead: use compact step names "rbind_N_sK" for intermediates and
    # "rbind_N" for the final result.  The 2-arg case keeps the old behaviour
    # (user-supplied out_dataset or R6 auto "A_rbind_B").
    n_args   <- length(args)
    
    # final_ds <- if (!is.null(out_dataset)) out_dataset
    # else if (n_args > 2L)      paste0("rbind_", n_args)
    # else                       NULL
    
    # final_ds <- if (!is.null(out_dataset)) out_dataset
    # else if (n_args > 2L)      paste0("rbind_", n_args, "_", .bind_hash(args))
    # else                       NULL
    # 
    # first_ds <- if (n_args > 2L) paste0("rbind_", n_args, "_s1") else final_ds
    
    uid      <- if (n_args > 2L) .bind_uid() else NULL
    final_ds <- if (!is.null(out_dataset)) out_dataset
    else if (n_args > 2L)      paste0("rbind_", n_args, "_", uid)
    else                       NULL
    first_ds <- if (n_args > 2L) paste0("rbind_", n_args, "_", uid, "_s1") else final_ds
    
    result <- args[[1L]]$rbind(
        args[[2L]],
        out_file    = out_file,
        out_group   = out_group,
        out_dataset = first_ds,
        block_rows  = block_rows,
        overwrite   = overwrite,
        compression = compression
    )
    
    if (n_args > 2L) {
        for (i in seq(3L, n_args)) {
            step_ds <- if (i == n_args) final_ds
            else             paste0("rbind_", n_args, "_s", i - 1L)
            result <- result$rbind(args[[i]],
                                   out_file    = out_file,
                                   out_group   = out_group,
                                   out_dataset = step_ds,
                                   block_rows  = block_rows,
                                   overwrite   = TRUE,
                                   compression = compression)
        }
    }
    
    result

    # result <- args[[1L]]$rbind(
    #     args[[2L]],
    #     out_file    = out_file,
    #     out_group   = out_group,
    #     out_dataset = out_dataset,
    #     block_rows  = block_rows,
    #     overwrite   = overwrite,
    #     compression = compression
    # )
    # 
    # if (length(args) > 2L) {
    #     for (i in seq(3L, length(args))) {
    #         result <- result$rbind(args[[i]],
    #                                out_file   = out_file,
    #                                out_group  = out_group,
    #                                block_rows = block_rows,
    #                                overwrite  = TRUE,
    #                                compression = compression)
    #     }
    # }
    # 
    # result
}
