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

    # ── Close own handle before C++ touches the same file ──────────────
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL

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

    # ── Close own handle before C++ touches the same file ──────────────
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL

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
