# HDF5Matrix_op_bind.R
#
# Adds $cbind() and $rbind() methods to the HDF5Matrix R6 class (Phase 10).
#
# Both methods follow the established pattern:
#   - Close own handle before calling C++ (string-based, same-file writes)
#   - C++ opens all handles internally (no conflicts)
#   - Reopen own handle after C++ finishes
#
# $cbind(other) -> HDF5Matrix  (same rows, appended columns)
# $rbind(other) -> HDF5Matrix  (same cols, appended rows)
# $bind_multi(others, func) -> HDF5Matrix  (internal: N-ary, single pass)


# ── $cbind() ────────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "cbind",
# @description
# Column-bind this matrix with another HDF5Matrix.
#'
# Appends the columns of \code{other} to the right of this matrix.
# Both matrices must have the same number of rows. The operation is
# performed block-wise on disk; no full matrix is loaded into RAM.
#'
# @param other      An \code{HDF5Matrix} with the same number of rows.
# @param out_file   Output HDF5 file path. \code{NULL} = same file as \code{self}.
# @param out_group  Output group. \code{NULL} = auto (\code{"BIND"}).
# @param out_dataset Output dataset name. \code{NULL} = auto.
# @param block_rows Integer. Rows per I/O block (default 1000).
# @param overwrite  Logical. Overwrite existing output. Default \code{FALSE}.
# @return \code{HDF5Matrix} pointing to the combined dataset.
function(other,
         out_file    = NULL,
         out_group   = NULL,
         out_dataset = NULL,
         block_rows  = 1000L,
         overwrite   = FALSE,
         compression = NULL) {

    if (!self$is_valid())  stop("HDF5Matrix is closed or invalid")
    if (!inherits(other, "HDF5Matrix") || !other$is_valid())
        stop("cbind: 'other' must be a valid HDF5Matrix")

    of  <- if (is.null(out_file))    private$filename else out_file
    og  <- if (is.null(out_group))   "BIND"           else out_group
    ods <- if (is.null(out_dataset))
               paste0(private$dataset, "_cbind_", other$.__enclos_env__$private$dataset)
           else out_dataset

    compression_eff <- .get_option("compression", default = NULL, override = compression)

    # ── Close own handle before C++ touches the same file; reopen on exit
    #    if C++ fails, so self stays valid ─────────────────────────────────
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    on.exit(if (is.null(private$ptr))
                private$ptr <- rcpp_hdf5dataset_open(
                    private$filename, private$group, private$dataset),
            add = TRUE)

    res <- rcpp_hdf5dataset_bind(
        private$filename, private$group, private$dataset,
        other$.__enclos_env__$private$filename,
        other$.__enclos_env__$private$group,
        other$.__enclos_env__$private$dataset,
        of, og, ods,
        func        = "cbind",
        overwrite   = isTRUE(overwrite),
        block_rows  = as.integer(block_rows),
        compression = compression_eff
    )

    # ── Reopen own handle ───────────────────────────────────────────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    hdf5_matrix(res$file, paste0(res$group, "/", res$dataset))
})


# ── $rbind() ────────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "rbind",
# @description
# Row-bind this matrix with another HDF5Matrix.
#'
# Appends the rows of \code{other} below this matrix. Both matrices must
# have the same number of columns. The operation is performed block-wise
# on disk; no full matrix is loaded into RAM.
#'
# @param other      An \code{HDF5Matrix} with the same number of columns.
# @param out_file   Output HDF5 file path. \code{NULL} = same file as \code{self}.
# @param out_group  Output group. \code{NULL} = auto (\code{"BIND"}).
# @param out_dataset Output dataset name. \code{NULL} = auto.
# @param block_rows Integer. Rows per I/O block (default 1000).
# @param overwrite  Logical. Overwrite existing output. Default \code{FALSE}.
# @return \code{HDF5Matrix} pointing to the combined dataset.
function(other,
         out_file    = NULL,
         out_group   = NULL,
         out_dataset = NULL,
         block_rows  = 1000L,
         overwrite   = FALSE,
         compression = NULL) {

    if (!self$is_valid())  stop("HDF5Matrix is closed or invalid")
    if (!inherits(other, "HDF5Matrix") || !other$is_valid())
        stop("rbind: 'other' must be a valid HDF5Matrix")

    of  <- if (is.null(out_file))    private$filename else out_file
    og  <- if (is.null(out_group))   "BIND"           else out_group
    ods <- if (is.null(out_dataset))
               paste0(private$dataset, "_rbind_", other$.__enclos_env__$private$dataset)
           else out_dataset

    compression_eff <- .get_option("compression", default = NULL, override = compression)

    # ── Close own handle before C++ touches the same file; reopen on exit
    #    if C++ fails, so self stays valid ─────────────────────────────────
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    on.exit(if (is.null(private$ptr))
                private$ptr <- rcpp_hdf5dataset_open(
                    private$filename, private$group, private$dataset),
            add = TRUE)

    res <- rcpp_hdf5dataset_bind(
        private$filename, private$group, private$dataset,
        other$.__enclos_env__$private$filename,
        other$.__enclos_env__$private$group,
        other$.__enclos_env__$private$dataset,
        of, og, ods,
        func        = "rbind",
        overwrite   = isTRUE(overwrite),
        block_rows  = as.integer(block_rows),
        compression = compression_eff
    )

    # ── Reopen own handle ───────────────────────────────────────────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    hdf5_matrix(res$file, paste0(res$group, "/", res$dataset))
})


# ── $bind_multi() ───────────────────────────────────────────────────────────

HDF5Matrix$set("public", "bind_multi",
# @description
# Internal: column- or row-bind this matrix with one or more HDF5Matrix
# objects in a single pass (used by cbind()/rbind() with 3+ arguments).
# The output is created once with its final dimensions; no intermediate
# datasets are written.
#'
# @param others     List of \code{HDF5Matrix} objects, bound after \code{self}
#   in order.
# @param func       \code{"cbind"} or \code{"rbind"}.
# @param out_file   Output HDF5 file path. \code{NULL} = same file as \code{self}.
# @param out_group  Output group. \code{NULL} = auto (\code{"BIND"}).
# @param out_dataset Output dataset name (required).
# @param block_rows Integer. Rows per I/O block (default 1000).
# @param overwrite  Logical. Overwrite existing output. Default \code{FALSE}.
# @return \code{HDF5Matrix} pointing to the combined dataset.
function(others,
         func        = c("cbind", "rbind"),
         out_file    = NULL,
         out_group   = NULL,
         out_dataset,
         block_rows  = 1000L,
         overwrite   = FALSE,
         compression = NULL) {

    func <- match.arg(func)
    if (!self$is_valid())  stop("HDF5Matrix is closed or invalid")
    for (i in seq_along(others)) {
        if (!inherits(others[[i]], "HDF5Matrix") || !others[[i]]$is_valid())
            stop(func, ": argument ", i + 1L, " must be a valid HDF5Matrix")
    }

    ops <- c(list(self), others)
    loc <- function(field)
        vapply(ops, function(m) m$.__enclos_env__$private[[field]], "")

    of  <- if (is.null(out_file))    private$filename else out_file
    og  <- if (is.null(out_group))   "BIND"           else out_group

    compression_eff <- .get_option("compression", default = NULL, override = compression)

    # ── Close own handle before C++ touches the same file; reopen on exit
    #    if C++ fails, so self stays valid ─────────────────────────────────
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    on.exit(if (is.null(private$ptr))
                private$ptr <- rcpp_hdf5dataset_open(
                    private$filename, private$group, private$dataset),
            add = TRUE)

    res <- rcpp_hdf5dataset_bind_multi(
        loc("filename"), loc("group"), loc("dataset"),
        of, og, out_dataset,
        func        = func,
        overwrite   = isTRUE(overwrite),
        block_rows  = as.integer(block_rows),
        compression = compression_eff
    )

    # ── Reopen own handle ───────────────────────────────────────────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    hdf5_matrix(res$file, paste0(res$group, "/", res$dataset))
})
