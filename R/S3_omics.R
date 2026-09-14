# S3_omics.R
#
# S3 generics and HDF5Matrix methods for omics-specific operations (Phase 12):
#   impute_snps()         — fill missing values in SNP data
#   filter_low_coverage() — remove high-missingness features
#   filter_maf()          — remove low-MAF SNPs
#
# Missing values follow the 0/1/2/3 genotype encoding of the omics layer:
# the value 3 marks a missing entry, R NAs are not recognised as such.


# ── impute_snps() ────────────────────────────────────────────────────────────

#' Impute missing SNP values in an HDF5Matrix
#'
#' @description
#' Fills the missing entries of SNP data by computing column or row means of the
#' non-missing values. Intended for 0/1/2-coded diploid genotype matrices in
#' which missing genotypes are stored as the value \code{3}; \code{NA} values
#' must be recoded to \code{3} before calling this function.
#'
#' Row and column names of \code{x}, when present, are carried over to the
#' result.
#'
#' @param x           An \code{HDF5Matrix} containing SNP data with missing
#'   values coded as \code{3}.
#' @param out_group   Output group. \code{NULL} = same as input (default).
#' @param out_dataset Output dataset name. \code{NULL} = same as input (default, in-place).
#' @param by_cols     Logical. Impute by columns (\code{TRUE}, default) or rows.
#' @param threads     Integer. Number of threads (-1 = auto).
#' @param overwrite   Logical. Overwrite existing output. Default \code{FALSE}.
#' @param ...         Ignored.
#' @return \code{HDF5Matrix} pointing to the imputed dataset, carrying the
#'   dimension names of \code{x}.
#'
#' @examples
#' \donttest{
#' set.seed(42)
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
#' Removes columns (SNPs) or rows (samples) that carry too many missing values.
#' Writes the result to a new dataset.
#'
#' Missing values are the entries equal to \code{3}, the missing-data code of
#' the 0/1/2/3 genotype encoding shared by the omics functions; \code{NA} or
#' \code{NaN} entries are \strong{not} counted, so recode them to \code{3}
#' beforehand. A column (or row) is removed when its proportion of \code{3}s is
#' greater than or equal to \code{pcent}.
#'
#' Row and column names of \code{x}, when present, are carried over to the
#' result: the names of the filtered axis are subset to the surviving elements.
#'
#' When \code{out_group}/\code{out_dataset} are \code{NULL} (default) the result
#' is written alongside the input dataset with the suffix \code{"_filtered"}.
#'
#' @param x           An \code{HDF5Matrix} containing SNP data with missing
#'   values coded as \code{3}.
#' @param out_group   Output group. \code{NULL} (default) = same group as input.
#' @param out_dataset Output dataset name. \code{NULL} (default) = input name + \code{"_filtered"}.
#' @param pcent       Numeric in \[0,1\]. Missing-data threshold
#'   (default \code{0.05}). Features whose proportion of \code{3}s reaches this
#'   value are removed.
#' @param by_cols     Logical. Filter columns (\code{TRUE}, default) or rows.
#' @param overwrite   Logical. Overwrite existing output. Default \code{FALSE}.
#' @param ...         Ignored.
#' @return \code{HDF5Matrix} pointing to the filtered dataset, carrying the
#'   dimension names of the surviving elements.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' fn <- tempfile(fileext = ".h5")
#' # 0/1/2 genotypes with about 10% missing values, coded as 3
#' snps <- matrix(sample(c(0, 1, 2, 3), 200, replace = TRUE,
#'                        prob = c(.3, .3, .3, .1)), 20, 10)
#' X   <- hdf5_create_matrix(fn, "geno/raw", data = snps)
#'
#' # Filter with auto output path (adds "_filtered" suffix)
#' out <- filter_low_coverage(X, pcent = 0.2)
#'
#' # Filter with explicit output
#' out2 <- filter_low_coverage(X, pcent = 0.2, out_group = "geno",
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
#' Removes the columns or rows whose Minor Allele Frequency (MAF) is at or below
#' \code{maf_threshold}, the usual way of dropping rare variants. Designed for
#' 0/1/2-coded diploid genotype matrices.
#'
#' MAF is computed as \code{maf = n0/n + 0.5 * n1/n}, taking \code{1 - maf} when
#' that value exceeds 0.5, where \code{n0} and \code{n1} count the 0s and 1s of
#' the feature. Only features with \code{maf > maf_threshold} are kept.
#'
#' Row and column names of \code{x}, when present, are carried over to the
#' result: the names of the filtered axis are subset to the surviving elements.
#'
#' When \code{out_group}/\code{out_dataset} are \code{NULL} (default) the
#' result is written alongside the input dataset with suffix \code{"_maf_filtered"}.
#'
#' @param x             An \code{HDF5Matrix} containing SNP data.
#' @param out_group     Output group. \code{NULL} (default) = same group as input.
#' @param out_dataset   Output dataset name. \code{NULL} (default) = input name + \code{"_maf_filtered"}.
#' @param maf_threshold Numeric in \[0, 0.5\]. MAF threshold (default \code{0.05}).
#'   SNPs with MAF at or **below** this value are removed.
#' @param by_cols       Logical. Treat SNPs as columns (\code{TRUE}) or as rows
#'   (\code{FALSE}, default).
#' @param block_size    Integer. Block size for I/O. Default \code{100L}.
#' @param overwrite     Logical. Overwrite existing output. Default \code{FALSE}.
#' @param ...           Ignored.
#' @return \code{HDF5Matrix} pointing to the filtered dataset, carrying the
#'   dimension names of the surviving elements.
#'
#' @examples
#' \donttest{
#' set.seed(42)
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
