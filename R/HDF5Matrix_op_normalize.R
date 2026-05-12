# HDF5Matrix_op_normalize.R
#
# Adds normalization / scaling methods to the HDF5Matrix R6 class.
#
# Public method added:
#   $normalize(center, scale, byrows, wsize, result_path)
#
# Follows the same pattern as hdf5_r6_arithmetic.cpp: passes filename /
# group / dataset strings to C++ instead of external pointers.  C++ opens
# all datasets internally with unique_ptr, so no HDF5 handle conflict can
# occur regardless of how many HDF5Matrix objects the user has open.


# ---------------------------------------------------------------------------
# R6 method: $normalize()
# ---------------------------------------------------------------------------

HDF5Matrix$set("public", "normalize",
# Block-wise normalization (centering and/or scaling) of the matrix.
# Equivalent to base R scale() but operates entirely on disk.
#
# @param center   Logical. Subtract column/row means. Default TRUE.
# @param scale    Logical. Divide by column/row SDs.  Default TRUE.
# @param byrows   Logical. If TRUE normalize row-wise; otherwise column-wise.
# @param wsize    Integer or NULL. Block size (NULL = auto).
# @param result_path Where to write the output dataset:
#   NULL (default): same file, path "NORMALIZED/<group>/<dataset>".
#   Character string "group/dataset": same file, given path.
#   Named list list(file="f.h5", path="group/ds"): different file.
# @return A new HDF5Matrix pointing to the normalized dataset, with
#   scaled:center and scaled:scale attributes attached.
function(center      = TRUE,
         scale       = TRUE,
         byrows      = FALSE,
         wsize       = NULL,
         result_path = NULL,
         compression = NULL) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # ── Resolve input location from private fields ─────────────────────────
    # Use private strings directly (same pattern as arithmetic) so C++ opens
    # its own handle — no ptr is passed, no handle conflict possible.
    in_file    <- private$filename
    in_group   <- private$group
    in_dataset <- private$dataset

    # ── Resolve output location ────────────────────────────────────────────
    if (is.null(result_path)) {
        # Default: mirror bdNormalize_hdf5 convention
        parts     <- strsplit(in_dataset, "/")[[1]]
        ds_name   <- parts[length(parts)]
        grp_name  <- if (length(parts) > 1)
                         paste(parts[-length(parts)], collapse = "/")
                     else in_group
        out_file  <- in_file
        out_path  <- paste0("NORMALIZED/", grp_name, "/", ds_name)
    } else if (is.character(result_path)) {
        out_file <- in_file
        out_path <- result_path
    } else if (is.list(result_path)) {
        if (is.null(result_path$file) || is.null(result_path$path))
            stop("result_path list must have 'file' and 'path' elements")
        out_file <- result_path$file
        out_path <- result_path$path
    } else {
        stop("result_path must be NULL, a character string, or list(file, path)")
    }

    # ── Split out_path into group and dataset name for C++ ─────────────────
    out_parts    <- strsplit(out_path, "/")[[1]]
    out_ds_name  <- out_parts[length(out_parts)]
    out_grp_name <- paste(out_parts[-length(out_parts)], collapse = "/")

    # ── Close own handle — C++ will open its own handles to the same file ──
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    
    # ── Call C++ block-wise normalizer ─────────────────────────────────────
    stats <- rcpp_hdf5dataset_normalize(
        in_file,
        in_group,
        in_dataset,
        out_file,
        out_grp_name,
        out_ds_name,
        center  = isTRUE(center),
        scale   = isTRUE(scale),
        byrows  = isTRUE(byrows),
        wsize   = wsize,
        compression = .get_option("compression", default = NULL, override = compression)
    )

    # ── Reopen after C++ has finished and released all its handles ──────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)
    
    # ── Open result and attach scaled:center / scaled:scale attributes ─────
    result <- hdf5_matrix(out_file, out_path)
    if (isTRUE(center) && length(stats$center) > 0)
        attr(result, "scaled:center") <- stats$center
    if (isTRUE(scale) && length(stats$scale) > 0)
        attr(result, "scaled:scale")  <- stats$scale

    result
})
