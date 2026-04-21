# S3_omics.R
#
# S3 generics and HDF5Matrix methods for omics-specific operations (Phase 12):
#   impute_snps()         — fill NAs in SNP data
#   filter_low_coverage() — remove high-missingness features
#   filter_maf()          — remove low-MAF SNPs


# ── impute_snps() ────────────────────────────────────────────────────────────

#' Impute missing SNP values in an HDF5Matrix
#'
#' @description
#' Fills NA entries in SNP data by computing column or row means of non-missing
#' values. Intended for 0/1/2-coded diploid genotype matrices.
#'
#' @param x           An \code{HDF5Matrix} containing SNP data with NAs.
#' @param out_group   Output group. \code{NULL} = same as input (default).
#' @param out_dataset Output dataset name. \code{NULL} = same as input (default, in-place).
#' @param by_cols     Logical. Impute by columns (\code{TRUE}, default) or rows.
#' @param threads     Integer. Number of threads (-1 = auto).
#' @param overwrite   Logical. Overwrite existing output. Default \code{FALSE}.
#' @param ...         Ignored.
#' @return \code{HDF5Matrix} pointing to the imputed dataset.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' 
#' # SNP data: 0/1/2 coded, 3 = missing (not NA)
#' snps <- matrix(sample(c(0L, 1L, 2L, 3L), 100 * 20,
#'                        replace = TRUE,
#'                        prob    = c(0.3, 0.3, 0.3, 0.1)),
#'                nrow = 100, ncol = 20)
#' 
#' X   <- hdf5_create_matrix(tmp, "geno/raw", data = snps)
#' imp <- impute_snps(X, out_group = "geno", out_dataset = "imputed")
#' dim(imp)
#' 
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @export
impute_snps <- function(x, ...) UseMethod("impute_snps")

#' @rdname impute_snps
#' @export
impute_snps.HDF5Matrix <- function(x,
                                    out_group   = NULL,
                                    out_dataset = NULL,
                                    by_cols     = TRUE,
                                    threads     = -1L,
                                    overwrite   = FALSE,
                                    ...) {
    x$impute_snps(out_group   = out_group,
                  out_dataset = out_dataset,
                  by_cols     = by_cols,
                  threads     = threads,
                  overwrite   = overwrite)
}


# ── filter_low_coverage() ────────────────────────────────────────────────────

#' Remove high-missingness features from an HDF5Matrix
#'
#' @description
#' Removes columns (SNPs) or rows (samples) whose proportion of missing values
#' (NAs) exceeds \code{pcent}. Writes result to a new dataset.
#'
#' When \code{out_group}/\code{out_dataset} are \code{NULL} (default) the result
#' is written alongside the input dataset with the suffix \code{"_filtered"}.
#'
#' @param x           An \code{HDF5Matrix} containing SNP data.
#' @param out_group   Output group. \code{NULL} (default) = same group as input.
#' @param out_dataset Output dataset name. \code{NULL} (default) = input name + \code{"_filtered"}.
#' @param pcent       Numeric in \[0,1\]. Maximum allowed NA proportion
#'   (default \code{0.05}). Features above this are removed.
#' @param by_cols     Logical. Filter columns (\code{TRUE}, default) or rows.
#' @param overwrite   Logical. Overwrite existing output. Default \code{FALSE}.
#' @param ...         Ignored.
#' @return \code{HDF5Matrix} pointing to the filtered dataset.
#'
#' @examples
#' \donttest{
#' fn <- tempfile(fileext = ".h5")
#' snps <- matrix(sample(c(0, 1, 2, NA), 200, replace = TRUE,
#'                        prob = c(.25, .25, .25, .25)), 20, 10)
#' X   <- hdf5_create_matrix(fn, "geno/raw", data = snps)
#' 
#' # Filter with auto output path (adds "_filtered" suffix)
#' out <- filter_low_coverage(X, pcent = 0.1)
#' 
#' # Filter with explicit output
#' out2 <- filter_low_coverage(X, out_group = "geno",
#'                              out_dataset = "filtered", overwrite = TRUE)
#' hdf5_close_all()
#' unlink(fn)
#' }
#'
#' @export
filter_low_coverage <- function(x, ...) UseMethod("filter_low_coverage")

#' @rdname filter_low_coverage
#' @export
filter_low_coverage.HDF5Matrix <- function(x,
                                            out_group   = NULL,
                                            out_dataset = NULL,
                                            pcent       = 0.05,
                                            by_cols     = TRUE,
                                            overwrite   = FALSE,
                                            ...) {
    res <- x$filter_low_coverage(out_group   = out_group,
                                  out_dataset = out_dataset,
                                  pcent       = pcent,
                                  by_cols     = by_cols,
                                  overwrite   = overwrite)
    # R6 returns list(result=HDF5Matrix, n_removed=int)
    # S3 unwraps to HDF5Matrix directly (consistent with other S3 operations)
    if (is.list(res) && inherits(res$result, "HDF5Matrix")) res$result else res
}


# ── filter_maf() ─────────────────────────────────────────────────────────────

#' Remove SNPs by Minor Allele Frequency from an HDF5Matrix
#'
#' @description
#' Removes columns or rows whose Minor Allele Frequency (MAF) exceeds
#' \code{maf_threshold}. Designed for 0/1/2-coded diploid genotype matrices.
#'
#' When \code{out_group}/\code{out_dataset} are \code{NULL} (default) the
#' result is written alongside the input dataset with suffix \code{"_maf_filtered"}.
#'
#' @param x             An \code{HDF5Matrix} containing SNP data.
#' @param out_group     Output group. \code{NULL} (default) = same group as input.
#' @param out_dataset   Output dataset name. \code{NULL} (default) = input name + \code{"_maf_filtered"}.
#' @param maf_threshold Numeric in \[0, 0.5\]. MAF threshold (default \code{0.05}).
#'   SNPs with MAF **above** this value are removed.
#' @param by_cols       Logical. Process by columns (\code{FALSE}, default) or rows.
#' @param block_size    Integer. Block size for I/O. Default \code{100L}.
#' @param overwrite     Logical. Overwrite existing output. Default \code{FALSE}.
#' @param ...           Ignored.
#' @return \code{HDF5Matrix} pointing to the filtered dataset.
#'
#' @examples
#' \donttest{
#' fn <- tempfile(fileext = ".h5")
#' snps <- matrix(sample(c(0, 1, 2), 200, replace = TRUE,
#'                        prob = c(.6, .3, .1)), 20, 10)
#' X   <- hdf5_create_matrix(fn, "geno/raw", data = snps)
#' 
#' # Filter with auto output path (adds "_maf_filtered" suffix)
#' out <- filter_maf(X, maf_threshold = 0.05)
#' 
#' # Filter with explicit output
#' out2 <- filter_maf(X, out_group = "geno",
#'                    out_dataset = "maf_filtered", overwrite = TRUE)
#' hdf5_close_all()
#' unlink(fn)
#' }
#'
#' @export
filter_maf <- function(x, ...) UseMethod("filter_maf")

#' @rdname filter_maf
#' @export
filter_maf.HDF5Matrix <- function(x,
                                   out_group     = NULL,
                                   out_dataset   = NULL,
                                   maf_threshold = 0.05,
                                   by_cols       = FALSE,
                                   block_size    = 100L,
                                   overwrite     = FALSE,
                                   ...) {
    res <- x$filter_maf(out_group     = out_group,
                         out_dataset   = out_dataset,
                         maf_threshold = maf_threshold,
                         by_cols       = by_cols,
                         block_size    = block_size,
                         overwrite     = overwrite)
    # Unwrap list → HDF5Matrix
    if (is.list(res) && inherits(res$result, "HDF5Matrix")) res$result else res
}
