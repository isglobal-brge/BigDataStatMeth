# HDF5Matrix_op_correlation.R
#
# Adds the $cor() method to the HDF5Matrix R6 class.
#
# Delegates to rcpp_hdf5dataset_cor() (hdf5_r6_correlation.cpp), which:
#   1. Runs bdCorr_hdf5() — always writes to in_file_x.
#   2. If out_file != in_file_x, copies result block-by-block via
#      hdf5Dataset::copyToFile() (defined in hdf5Datasets.hpp).
#
# Pattern: passes filename/group/dataset strings — never external pointers.


HDF5Matrix$set("public", "cor",
# @description
# Block-wise correlation of an HDF5Matrix.
# Computes \code{cor(self)} or \code{cor(self, y)} entirely on disk.
#'
# @param y         An \code{HDF5Matrix} for cross-correlation, or \code{NULL}
#                  (default) for auto-correlation.
# @param method    \code{"pearson"} (default) or \code{"spearman"}.
# @param use       \code{"everything"} (default) or \code{"complete.obs"}.
# @param trans_x   Logical. Correlate rows of self instead of columns.
# @param trans_y   Logical. Correlate rows of y instead of columns.
# @param compute_pvalues Logical. Compute and store p-values. Default TRUE.
# @param block_size Integer or NULL. Block size for computation (NULL = 1000).
# @param threads   Ignored (reserved for future OpenMP support).
# @param result_path Output location:
#   \code{NULL} — same file, auto group \code{"CORR/<dataset>"}.
#   Character — same file, that string as output group.
#   List \code{list(file=, group=)} — different file and/or group.
# @return A new \code{HDF5Matrix} pointing to the correlation matrix.
function(y               = NULL,
         method          = "pearson",
         use             = "everything",
         trans_x         = FALSE,
         trans_y         = FALSE,
         compute_pvalues = TRUE,
         block_size      = NULL,
         threads         = NULL,
         result_path     = NULL,
         compression     = NULL) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    method <- match.arg(method, c("pearson", "spearman"))
    use    <- match.arg(use, c("everything", "complete.obs",
                               "all.obs", "pairwise.complete.obs",
                               "na.or.complete"))
    use_complete <- use %in% c("complete.obs", "pairwise.complete.obs",
                               "na.or.complete")

    # ── Input (strings only — no ptr) ─────────────────────────────────────
    in_file_x  <- private$filename
    in_group_x <- private$group
    in_ds_x    <- private$dataset

    in_file_y  <- ""
    in_group_y <- ""
    in_ds_y    <- ""
    if (!is.null(y)) {
        if (!inherits(y, "HDF5Matrix")) stop("y must be an HDF5Matrix or NULL")
        if (!y$is_valid())             stop("y is closed or invalid")
        in_file_y  <- y$.__enclos_env__$private$filename
        in_group_y <- y$.__enclos_env__$private$group
        in_ds_y    <- y$.__enclos_env__$private$dataset
    }

    # ── Resolve output location ────────────────────────────────────────────
    out_file  <- ""   # "" -> rcpp_hdf5dataset_cor uses in_file_x
    out_group <- ""   # "" -> auto: "CORR/<dataset>"
    if (!is.null(result_path)) {
        if (is.character(result_path)) {
            out_group <- result_path
        } else if (is.list(result_path)) {
            if (!is.null(result_path$file))  out_file  <- result_path$file
            if (!is.null(result_path$group)) out_group <- result_path$group
        } else {
            stop("result_path must be NULL, a character string, or list(file=, group=)")
        }
    }

    bs <- if (is.null(block_size)) 1000L else as.integer(block_size)

    # ── Close own handle — C++ will open its own handles to the same file ──
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    
    compression_eff <- .get_option("compression", default = NULL, override = compression)
    
    # ── Call C++ wrapper ─────────────────────────────────────────────────
    res <- rcpp_hdf5dataset_cor(
        in_file_x, in_group_x, in_ds_x,
        in_file_y, in_group_y, in_ds_y,
        out_file, out_group,
        isTRUE(trans_x), isTRUE(trans_y),
        method, use_complete,
        isTRUE(compute_pvalues),
        bs,
        500L,
        threads,
        compression_eff
    )

    # ── Reopen after C++ has finished and released all its handles ──────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)
    
    # ── Open result as HDF5Matrix ─────────────────────────────────────────
    result_file <- res$filename
    result_full <- paste(res$group, res$correlation, sep = "/")

    result <- hdf5_matrix(result_file, result_full)

    attr(result, "cor.method")  <- res$method
    attr(result, "cor.type")    <- res$correlation_type
    attr(result, "cor.n.vars")  <- res$n_variables
    attr(result, "cor.n.obs")   <- res$n_observations
    if (isTRUE(compute_pvalues) && nchar(res$pvalues) > 0)
        attr(result, "cor.pvalues.path") <- paste(res$group, res$pvalues, sep = "/")

    result
})
