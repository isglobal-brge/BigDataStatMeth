# HDF5Matrix_op_omics.R
#
# Adds omics-specific methods to the HDF5Matrix R6 class (Phase 12):
#   $impute_snps()         — fill missing values (NAs) using column/row means
#   $filter_low_coverage() — remove columns/rows with too many NAs
#   $filter_maf()          — remove columns/rows by Minor Allele Frequency
#
# All methods follow Rule 1.7: close own handle before calling C++,
# reopen after. Delegates to thin rcpp_hdf5dataset_* wrappers in
# hdf5_r6_omics.cpp, which call the header-level inline functions directly.


# ── $impute_snps() ───────────────────────────────────────────────────────────

HDF5Matrix$set("public", "impute_snps",
# @description
# Impute missing values (NAs) in SNP data stored in HDF5.
#
# Fills NA entries by computing column or row means of non-missing values
# (for 0/1/2-coded genotype data). By default the result overwrites the
# input dataset.
#
# @param out_group   Output group.   \code{NULL} = same as input.
# @param out_dataset Output dataset. \code{NULL} = same as input (in-place).
# @param by_cols     Logical. Impute by columns (\code{TRUE}, default) or rows.
# @param threads     Integer. Number of OpenMP threads (-1 = auto).
# @param overwrite   Logical. Overwrite existing output. Default \code{FALSE}.
# @return \code{HDF5Matrix} pointing to the imputed dataset.
function(out_group   = NULL,
         out_dataset = NULL,
         by_cols     = TRUE,
         threads     = -1L,
         overwrite   = FALSE) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    og  <- if (is.null(out_group))   private$group   else out_group
    ods <- if (is.null(out_dataset)) private$dataset else out_dataset

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

    hdf5_matrix(res$file, paste0(res$group, "/", res$dataset))
})


# ── $filter_low_coverage() ──────────────────────────────────────────────────

HDF5Matrix$set("public", "filter_low_coverage",
# @description
# Remove columns or rows whose NA proportion exceeds \code{pcent}.
#
# For SNP data: removes SNPs (columns) or samples (rows) that exceed the
# missing-data threshold. Always writes to a new dataset.
#
# @param out_group   Output group (required).
# @param out_dataset Output dataset name (required).
# @param pcent       Numeric in \[0,1\]. Maximum allowed NA proportion.
#   Default \code{0.05} (5\%). Features above this threshold are removed.
# @param by_cols     Logical. Filter by columns (\code{TRUE}, default) or rows.
# @param overwrite   Logical. Overwrite existing output. Default \code{FALSE}.
# @return Named list:
#   \describe{
#     \item{result}{\code{HDF5Matrix} pointing to the filtered dataset.}
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

    list(
        result    = hdf5_matrix(res$file, paste0(res$group, "/", res$dataset)),
        n_removed = as.integer(res$n_removed)
    )
})


# ── $filter_maf() ────────────────────────────────────────────────────────────

HDF5Matrix$set("public", "filter_maf",
# @description
# Remove SNPs whose Minor Allele Frequency exceeds \code{maf_threshold}.
#
# MAF formula for 0/1/2-coded diploid data:
#   maf = (n0/n) + 0.5*(n1/n); if maf > 0.5: maf = 1 - maf.
# SNPs with maf > maf_threshold are removed. Always writes to a new dataset.
#
# @param out_group     Output group (required).
# @param out_dataset   Output dataset name (required).
# @param maf_threshold Numeric in \[0, 0.5\]. Default \code{0.05}.
#   SNPs with MAF **above** this threshold are removed.
# @param by_cols       Logical. Process by columns (\code{FALSE}, default =
#   SNPs are rows) or by columns.
# @param block_size    Integer. I/O block size. Default \code{100L}.
# @param overwrite     Logical. Overwrite existing output. Default \code{FALSE}.
# @return Named list:
#   \describe{
#     \item{result}{\code{HDF5Matrix} pointing to the filtered dataset.}
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

    list(
        result    = hdf5_matrix(res$file, paste0(res$group, "/", res$dataset)),
        n_removed = as.integer(res$n_removed)
    )
})
