# HDF5Matrix_op_reduce_apply.R
#
# Adds $reduce() and $apply_function() to the HDF5Matrix R6 class.
# Both operate on HDF5 groups (lists of datasets), so they take filename/group
# strings and do NOT need an active XPtr.


# -- $reduce() ----------------------------------------------------------------

HDF5Matrix$set("public", "reduce",
# @description
# Reduces all datasets in the group of this HDF5Matrix by applying a binary
# accumulation operation ("+" or "-") across them.
#
# All datasets in the group must have the same dimensions.  The result is
# stored as a single dataset in out_group/out_dataset.
#
# This is the HDF5-group reduce (e.g. sum a collection of partial
# results that have been split or parallelised), not a row/column summary
# (use colSums(), rowSums(), etc. for those).
#
# @param out_group    Character. Output HDF5 group (default "REDUCED").
# @param out_dataset  Character. Output dataset name (default "reduced").
# @param func         Character. Accumulation function: "+" (default) or "-".
# @param overwrite    Logical. Overwrite existing output (default FALSE).
# @param remove_input Logical. Remove input datasets after reduction (default FALSE).
# @return A new HDF5Matrix pointing to the reduced dataset.
function(out_group    = "REDUCED",
         out_dataset  = "reduced",
         func         = "+",
         overwrite    = FALSE,
         remove_input = FALSE) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")
    func <- match.arg(func, c("+", "-"))

    # String-based call - close own handle to avoid HDF5 file conflicts
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    on.exit({
        if (is.null(private$ptr))
            private$ptr <- rcpp_hdf5dataset_open(
                private$filename, private$group, private$dataset)
    }, add = TRUE)

    res <- rcpp_hdf5dataset_reduce(
        filename     = private$filename,
        group        = private$group,
        out_group    = out_group,
        out_dataset  = out_dataset,
        func         = func,
        overwrite    = isTRUE(overwrite),
        remove_input = isTRUE(remove_input)
    )

    hdf5_matrix(res$filename, res$path)
})


# -- $apply_function() --------------------------------------------------------

HDF5Matrix$set("public", "apply_function",
# @description
# Applies a statistical or algebraic function to one or more named datasets
# in the same group as this HDF5Matrix and writes results to out_group.
#
# Valid func values:
#   "QR"               - QR decomposition
#   "CrossProd"        - Cross product A'A
#   "tCrossProd"       - Transposed cross product AA'
#   "invChol"          - Matrix inverse via Cholesky
#   "blockmult"        - Block matrix multiply (requires b_datasets)
#   "CrossProd_double" - Double-precision cross product
#   "tCrossProd_double"- Double-precision transposed cross product
#   "solve"            - Matrix solve
#   "normalize"        - Column normalisation
#   "sdmean"           - Standard deviation and mean
#   "descChol"         - Cholesky decomposition
#
# @param datasets     Character vector. Dataset names in group to process.
#   NULL (default) = use self only.
# @param func         Character. Function to apply (see above).
# @param out_group    Character. Output group (default "APPLIED").
# @param b_group      Character. Group for second-operand datasets (for "blockmult").
# @param b_datasets   Character vector or NULL. Second-operand dataset names.
# @param overwrite    Logical. Overwrite existing results (default FALSE).
# @param transp_a     Logical. Transpose input A (default FALSE).
# @param transp_b     Logical. Transpose input B (default FALSE).
# @param full_matrix  Logical. Return full symmetric matrix (default FALSE).
# @param byrows       Logical. Apply by rows (default FALSE).
# @param threads      Integer or NULL. OpenMP threads.
# @return Named list: filename, out_group, func, datasets.
function(datasets    = NULL,
         func        = "QR",
         out_group   = "APPLIED",
         b_group     = "",
         b_datasets  = NULL,
         overwrite   = FALSE,
         transp_a    = FALSE,
         transp_b    = FALSE,
         full_matrix = FALSE,
         byrows      = FALSE,
         threads     = NULL) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    valid_funcs <- c("QR", "CrossProd", "tCrossProd",
                     "invChol", "blockmult", "CrossProd_double",
                     "tCrossProd_double", "solve", "normalize",
                     "sdmean", "descChol")
    func <- match.arg(func, valid_funcs)

    ds <- if (is.null(datasets)) private$dataset else as.character(datasets)
    threads_eff <- .get_option("threads", default = NULL, override = threads)

    # # String-based call - close own handle
    # rcpp_hdf5dataset_close(private$ptr)
    # private$ptr <- NULL
    # on.exit({
    #     if (is.null(private$ptr))
    #         private$ptr <- rcpp_hdf5dataset_open(
    #             private$filename, private$group, private$dataset)
    # }, add = TRUE)

    res <- rcpp_hdf5dataset_apply_function(
        filename    = private$filename,
        group       = private$group,
        datasets    = ds,
        out_group   = out_group,
        func        = func,
        b_group     = b_group,
        b_datasets  = b_datasets,
        overwrite   = isTRUE(overwrite),
        transp_a    = isTRUE(transp_a),
        transp_b    = isTRUE(transp_b),
        full_matrix = isTRUE(full_matrix),
        byrows      = isTRUE(byrows),
        threads     = threads_eff
    )
    
    # # Rename outputs: header names them <ds>, test expects <func>_<ds>
    # for (d in basename(ds)) {
    #     src <- paste0(out_group, "/", d)
    #     dst <- paste0(out_group, "/", func, "_", d)
    #     tryCatch(
    #         bdmove_hdf5_dataset(private$filename, src, dst,
    #                             overwrite = isTRUE(overwrite)),
    #         error = function(e) invisible(NULL)
    #     )
    # }

    res
})
