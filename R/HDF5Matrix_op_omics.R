# HDF5Matrix_op_omics.R
#
# Adds omics-specific methods to the HDF5Matrix R6 class (Phase 12):
#   $impute_snps()         — fill missing values using column/row means
#   $filter_low_coverage() — remove columns/rows with too many missing values
#   $filter_maf()          — remove columns/rows by Minor Allele Frequency
#
# Missing values follow the 0/1/2/3 genotype encoding of the omics layer:
# the value 3 marks a missing entry, R NAs are not recognised as such.
#
# All methods follow Rule 1.7: close own handle before calling C++,
# reopen after. Delegates to thin rcpp_hdf5dataset_* wrappers in
# hdf5_r6_omics.cpp, which call the header-level inline functions directly.
#
# All three methods propagate the dimension names of the input dataset to the
# result. The two filtering methods use the `kept` field returned by their
# wrapper (1-based surviving positions along the filtered axis) to subset the
# names of the axis that shrank.


# ── dimnames propagation helper ──────────────────────────────────────────────

# Copy dimension names onto the dataset wrapped by an HDF5Matrix.
#
# @param target   HDF5Matrix whose dataset receives the names.
# @param src      List with elements `rownames` / `colnames`, as returned by
#                 rcpp_hdf5dataset_read_dimnames(). Empty vectors mean the
#                 corresponding dimension has no names.
# @param keep_rows Integer vector of 1-based surviving row positions, or NULL
#                 when the row count did not change.
# @param keep_cols Integer vector of 1-based surviving column positions, or
#                 NULL when the column count did not change.
# @return NULL, invisibly. Nothing is written when the input carried no names.
# @keywords internal
# @noRd
.hdf5_propagate_dimnames <- function(target, src,
                                     keep_rows = NULL, keep_cols = NULL) {

    if (is.null(src)) return(invisible(NULL))

    rn <- src$rownames
    cn <- src$colnames
    if (is.null(rn)) rn <- character(0)
    if (is.null(cn)) cn <- character(0)

    if (length(rn) > 0 && !is.null(keep_rows)) rn <- rn[keep_rows]
    if (length(cn) > 0 && !is.null(keep_cols)) cn <- cn[keep_cols]

    # Nothing to inherit - leave the result without dimnames
    if (length(rn) == 0 && length(cn) == 0) return(invisible(NULL))

    rcpp_hdf5dataset_write_dimnames(
        target$.__enclos_env__$private$ptr,
        as.character(rn),
        as.character(cn))

    invisible(NULL)
}


# ── "filter removed everything" guard ────────────────────────────────────────

# Report whether a filter left no output dataset behind.
#
# The header-level filters only create the output dataset once a block with at
# least one surviving element is written: when the criterion removes every
# column/row nothing is created (Rcpp_Remove_Low_Data_hdf5() warns "All data
# removed ...", Rcpp_Remove_MAF_hdf5() returns silently), yet the wrapper still
# reports the intended output location. Opening it afterwards fails deep inside
# the C++ layer with the unhelpful "please create Dataset before proceed", so
# the R6 methods check here first and stop with an actionable message.
#
# Two independent signals, either one is conclusive:
#   * `kept` empty  — no element survived the filter, so nothing was written;
#   * the dataset is absent from its group in the file (also covers the case
#     where `overwrite = TRUE` wiped an earlier output that was not rewritten).
#
# @param filename Path of the HDF5 file holding the output.
# @param group    Output group as reported by the wrapper.
# @param dataset  Output dataset name as reported by the wrapper.
# @param kept     Integer vector of surviving 1-based positions.
# @return TRUE when the output dataset is missing.
# @keywords internal
# @noRd
.hdf5_filter_output_missing <- function(filename, group, dataset, kept) {

    if (length(kept) == 0L) return(TRUE)

    present <- tryCatch(
        bdgetDatasetsList_hdf5(filename, group = group),
        error = function(e) character(0))

    !(dataset %in% present)
}


# ── $impute_snps() ───────────────────────────────────────────────────────────

HDF5Matrix$set("public", "impute_snps",
# @description
# Impute missing values in SNP data stored in HDF5.
#
# Fills the entries coded as 3 by computing column or row means of the
# non-missing values (for 0/1/2/3-coded genotype data). By default the
# result overwrites the input dataset.
#
# @param out_group   Output group.   \code{NULL} = same as input.
# @param out_dataset Output dataset. \code{NULL} = same as input (in-place).
# @param by_cols     Logical. Impute by columns (\code{TRUE}, default) or rows.
# @param threads     Integer. Number of OpenMP threads (-1 = auto).
# @param overwrite   Logical. Overwrite existing output. Default \code{FALSE}.
# @return \code{HDF5Matrix} pointing to the imputed dataset. Dimension names
#   of the input are carried over (imputation does not change dimensions).
function(out_group   = NULL,
         out_dataset = NULL,
         by_cols     = TRUE,
         threads     = -1L,
         overwrite   = FALSE) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    og  <- if (is.null(out_group))   private$group   else out_group
    ods <- if (is.null(out_dataset)) private$dataset else out_dataset

    # Read the dimnames while our own handle is still open
    src_dn <- rcpp_hdf5dataset_read_dimnames(private$ptr)

    # ── Close own handle — C++ opens its own handles ────────────────────
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL

    res <- rcpp_hdf5dataset_impute_snps(
        in_file    = private$filename,
        in_group   = private$group,
        in_dataset = private$dataset,
        out_group  = og,
        out_dataset = ods,
        by_cols    = isTRUE(by_cols),
        threads    = as.integer(threads),
        overwrite  = isTRUE(overwrite)
    )

    # ── Reopen own handle ───────────────────────────────────────────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    out <- hdf5_matrix(res$file, paste0(res$group, "/", res$dataset))

    # In-place imputation keeps the names already stored for the dataset;
    # a new output dataset has to inherit them. Dimensions are unchanged,
    # so both axes are copied as-is.
    if (!identical(og, private$group) || !identical(ods, private$dataset))
        .hdf5_propagate_dimnames(out, src_dn)

    out
})


# ── $filter_low_coverage() ──────────────────────────────────────────────────

HDF5Matrix$set("public", "filter_low_coverage",
# @description
# Remove columns or rows with too many missing values.
#
# Missing values are the entries equal to 3, the missing-data code of the
# 0/1/2/3 encoding used by the omics functions; \code{NA} entries are not
# counted. A column or row is removed when its proportion of 3s is greater
# than or equal to \code{pcent}. Always writes to a new dataset.
#
# @param out_group   Output group (required).
# @param out_dataset Output dataset name (required).
# @param pcent       Numeric in \[0,1\]. Missing-data threshold.
#   Default \code{0.05} (5\%). Features reaching this proportion are removed.
# @param by_cols     Logical. Filter by columns (\code{TRUE}, default) or rows.
# @param overwrite   Logical. Overwrite existing output. Default \code{FALSE}.
# @return Named list:
#   \describe{
#     \item{result}{\code{HDF5Matrix} pointing to the filtered dataset, with
#       the dimension names of the surviving elements.}
#     \item{n_removed}{Integer. Number of columns/rows removed.}
#   }
function(out_group   = NULL,
         out_dataset = NULL,
         pcent       = 0.05,
         by_cols     = TRUE,
         overwrite   = FALSE) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # Auto-derive output location when not specified (same pattern as impute_snps)
    og  <- if (is.null(out_group))   private$group                        else out_group
    ods <- if (is.null(out_dataset)) paste0(private$dataset, "_filtered") else out_dataset

    # Read the dimnames while our own handle is still open
    src_dn <- rcpp_hdf5dataset_read_dimnames(private$ptr)

    # ── Close own handle ────────────────────────────────────────────────
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL

    res <- rcpp_hdf5dataset_filter_low_coverage(
        in_file     = private$filename,
        in_group    = private$group,
        in_dataset  = private$dataset,
        out_group   = og,
        out_dataset = ods,
        pcent       = as.double(pcent),
        by_cols     = isTRUE(by_cols),
        overwrite   = isTRUE(overwrite)
    )

    # ── Reopen own handle ───────────────────────────────────────────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    kept <- as.integer(res$kept)

    # Nothing survived the filter: the C++ layer never created the output
    # dataset, so fail here instead of letting hdf5_matrix() report the
    # internal "please create Dataset before proceed".
    if (.hdf5_filter_output_missing(res$file, res$group, res$dataset, kept))
        stop("filter_low_coverage() removed all ",
             if (isTRUE(by_cols)) "columns" else "rows",
             ": no output dataset was created. ",
             "Raise `pcent` or review the input data ",
             "(missing values must be coded as 3, not NA).",
             call. = FALSE)

    out <- hdf5_matrix(res$file, paste0(res$group, "/", res$dataset))

    # by_cols = TRUE filters HDF5 dimension 0, which is the R column axis;
    # by_cols = FALSE filters HDF5 dimension 1, the R row axis. Only the axis
    # that shrank gets its names subset.
    if (isTRUE(by_cols))
        .hdf5_propagate_dimnames(out, src_dn, keep_cols = kept)
    else
        .hdf5_propagate_dimnames(out, src_dn, keep_rows = kept)

    list(
        result    = out,
        n_removed = as.integer(res$n_removed)
    )
})


# ── $filter_maf() ────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "filter_maf",
# @description
# Remove SNPs whose Minor Allele Frequency is at or below \code{maf_threshold}.
#
# MAF formula for 0/1/2-coded diploid data:
#   maf = (n0/n) + 0.5*(n1/n); if maf > 0.5: maf = 1 - maf.
# SNPs with maf <= maf_threshold are removed, so only SNPs whose MAF is
# strictly above the threshold are kept. Always writes to a new dataset.
#
# @param out_group     Output group (required).
# @param out_dataset   Output dataset name (required).
# @param maf_threshold Numeric in \[0, 0.5\]. Default \code{0.05}.
#   SNPs with MAF at or **below** this threshold are removed.
# @param by_cols       Logical. Treat SNPs as columns (\code{TRUE}) or as rows
#   (\code{FALSE}, default).
# @param block_size    Integer. I/O block size. Default \code{100L}.
# @param overwrite     Logical. Overwrite existing output. Default \code{FALSE}.
# @return Named list:
#   \describe{
#     \item{result}{\code{HDF5Matrix} pointing to the filtered dataset, with
#       the dimension names of the surviving elements.}
#     \item{n_removed}{Integer. Number of SNPs removed.}
#   }
function(out_group     = NULL,
         out_dataset   = NULL,
         maf_threshold = 0.05,
         by_cols       = FALSE,
         block_size    = 100L,
         overwrite     = FALSE) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # Auto-derive output location when not specified
    og  <- if (is.null(out_group))   private$group                           else out_group
    ods <- if (is.null(out_dataset)) paste0(private$dataset, "_maf_filtered") else out_dataset

    # Read the dimnames while our own handle is still open
    src_dn <- rcpp_hdf5dataset_read_dimnames(private$ptr)

    # ── Close own handle ────────────────────────────────────────────────
    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL

    res <- rcpp_hdf5dataset_filter_maf(
        in_file       = private$filename,
        in_group      = private$group,
        in_dataset    = private$dataset,
        out_group     = og,
        out_dataset   = ods,
        maf_threshold = as.double(maf_threshold),
        by_cols       = isTRUE(by_cols),
        block_size    = as.integer(block_size),
        overwrite     = isTRUE(overwrite)
    )

    # ── Reopen own handle ───────────────────────────────────────────────
    private$ptr <- rcpp_hdf5dataset_open(
        private$filename, private$group, private$dataset)

    kept <- as.integer(res$kept)

    # Nothing survived the filter: the C++ layer never created the output
    # dataset, so fail here instead of letting hdf5_matrix() report the
    # internal "please create Dataset before proceed".
    if (.hdf5_filter_output_missing(res$file, res$group, res$dataset, kept))
        stop("filter_maf() removed all ",
             if (isTRUE(by_cols)) "columns" else "rows",
             ": no output dataset was created. ",
             "Lower `maf_threshold` or review the input data ",
             "(genotypes must be 0/1/2 coded).",
             call. = FALSE)

    out <- hdf5_matrix(res$file, paste0(res$group, "/", res$dataset))

    # by_cols = TRUE filters HDF5 dimension 0, which is the R column axis;
    # by_cols = FALSE filters HDF5 dimension 1, the R row axis. Only the axis
    # that shrank gets its names subset.
    if (isTRUE(by_cols))
        .hdf5_propagate_dimnames(out, src_dn, keep_cols = kept)
    else
        .hdf5_propagate_dimnames(out, src_dn, keep_rows = kept)

    list(
        result    = out,
        n_removed = as.integer(res$n_removed)
    )
})
