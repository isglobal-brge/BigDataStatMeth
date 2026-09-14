# =============================================================================
# Regression test: dimname inheritance across the omics operations
# -----------------------------------------------------------------------------
# THE BUG THIS PINS DOWN. impute_snps(), filter_low_coverage() and filter_maf()
# wrote their result to a brand-new HDF5 dataset and stopped there: the hidden
# `.<dataset>_dimnames` group of the input was never carried over, so the result
# came back with dimnames() == NULL even when the input had row and column
# names. Every downstream step that matches features or samples by name (the
# usual reason those names exist) silently lost its keys.
#
# THE FIX. The two header-level filters now report, through an optional output
# vector, the 0-based indices along the filtered axis of the elements that
# SURVIVED; the Rcpp wrappers hand them to R as the 1-based `kept` field, and
# the R6 layer subsets the names of the axis that shrank before writing them on
# the result. impute_snps() does not change dimensions, so it copies both axes
# verbatim (and skips the write entirely when it imputed in place, where the
# names are already stored).
#
# THE TRANSPOSITION. HDF5 stores R matrices transposed: HDF5 dimension 0 is the
# R COLUMN axis. The filters iterate over HDF5 dimension 0 when by_cols = TRUE,
# so `by_cols = TRUE` shrinks the R column axis and its `kept` indices are R
# column positions; `by_cols = FALSE` shrinks the R row axis. Sections 2 and 3
# assert both axes against an in-memory computation of what should survive.
#
# THE MISSING-DATA CODE (see section 5). filter_low_coverage() counts the
# entries equal to the literal value 3 -- the missing code of the 0/1/2/3
# genotype encoding shared by the omics layer -- and NOT R NAs. A matrix whose
# gaps are NA is filtered as if it had no missing data at all. That is the real
# semantics, not a bug; the documentation used to say "NA proportion" and has
# been corrected. The same holds for the comparison being >= pcent, not > pcent.
#
# It lives in tests_regression/, git-tracked (reaches GitHub) but excluded from
# the built tarball via .Rbuildignore (^test\.*), so it does not run in
# R CMD check / ship to CRAN -- matching this repo's convention.
# Run with:  Rscript tests_regression/test_omics_dimnames.R
# =============================================================================

suppressMessages({
    library(BigDataStatMeth)
    library(testthat)
})

cat("BigDataStatMeth version:",
    as.character(packageVersion("BigDataStatMeth")), "\n\n")

# --- reference implementations of the two filter criteria --------------------

# Proportion of missing entries (the value 3) per column / per row.
# Unnamed on purpose: the expectations below compare bare index vectors.
miss_frac <- function(m, by_cols) {
    unname(if (by_cols) colMeans(m == 3) else rowMeans(m == 3))
}

# MAF exactly as calc_freq() computes it (hdf5Omics/hdf5OmicsUtils.hpp).
calc_freq_r <- function(v) {
    n   <- length(v)
    maf <- sum(v == 0) / n + 0.5 * (sum(v == 1) / n)
    if (maf > 0.5) maf <- 1 - maf
    maf
}
maf_per <- function(m, by_cols) {
    unname(apply(m, if (by_cols) 2L else 1L, calc_freq_r))
}

RN <- function(n) paste0("sample_", seq_len(n))
CN <- function(n) paste0("snp_",    seq_len(n))


# -----------------------------------------------------------------------------
# 1. impute_snps() keeps the dimnames identical (dimensions do not change)
# -----------------------------------------------------------------------------
test_that("impute_snps inherits the dimnames unchanged", {
    f <- tempfile(fileext = ".h5")
    set.seed(101)
    m <- matrix(sample(c(0, 1, 2), 20 * 10, TRUE), 20, 10)
    m[3, 4] <- 3; m[5, 6] <- 3; m[12, 4] <- 3   # missing entries to impute
    dimnames(m) <- list(RN(20), CN(10))

    X <- hdf5_create_matrix(f, "geno/raw", data = m, overwrite = TRUE)

    # by_cols = TRUE and FALSE both write to a NEW dataset here
    for (bc in c(TRUE, FALSE)) {
        out <- impute_snps(X, out_group = "geno",
                           out_dataset = paste0("imp_", bc),
                           by_cols = bc, overwrite = TRUE)
        expect_equal(dim(out), c(20L, 10L))
        expect_equal(rownames(out), RN(20))
        expect_equal(colnames(out), CN(10))
        out$close()
    }

    # In-place imputation must not lose the names it already had
    inp <- impute_snps(X, overwrite = TRUE)
    expect_equal(rownames(inp), RN(20))
    expect_equal(colnames(inp), CN(10))
    inp$close()

    X$close()
    unlink(f)
})


# -----------------------------------------------------------------------------
# 2. filter_low_coverage() keeps exactly the survivors' names, on both axes
# -----------------------------------------------------------------------------
test_that("filter_low_coverage inherits the surviving dimnames (columns)", {
    f <- tempfile(fileext = ".h5")
    set.seed(102)
    m <- matrix(sample(c(0, 1, 2), 20 * 10, TRUE), 20, 10)
    m[1:15, 3] <- 3          # 75% missing
    m[1:15, 7] <- 3          # 75% missing
    m[1,      9] <- 3        #  5% missing -> exactly at the >= threshold
    dimnames(m) <- list(RN(20), CN(10))

    X   <- hdf5_create_matrix(f, "geno/raw", data = m, overwrite = TRUE)
    out <- filter_low_coverage(X, out_group = "geno", out_dataset = "lowc",
                               pcent = 0.05, by_cols = TRUE, overwrite = TRUE)

    keep <- which(miss_frac(m, by_cols = TRUE) < 0.05)
    expect_equal(keep, c(1L, 2L, 4L, 5L, 6L, 8L, 10L))   # 3, 7 and 9 dropped

    expect_equal(dim(out), c(20L, length(keep)))
    expect_equal(unname(as.matrix(out)), unname(m[, keep]))

    # Rows untouched -> full rownames; columns filtered -> subset colnames
    expect_equal(rownames(out), RN(20))
    expect_equal(colnames(out), CN(10)[keep])

    out$close(); X$close()
    unlink(f)
})

test_that("filter_low_coverage inherits the surviving dimnames (rows)", {
    f <- tempfile(fileext = ".h5")
    set.seed(103)
    m <- matrix(sample(c(0, 1, 2), 20 * 10, TRUE), 20, 10)
    m[4,  1:8] <- 3          # 80% missing
    m[11, 1:8] <- 3          # 80% missing
    dimnames(m) <- list(RN(20), CN(10))

    X   <- hdf5_create_matrix(f, "geno/raw", data = m, overwrite = TRUE)
    out <- filter_low_coverage(X, out_group = "geno", out_dataset = "lowr",
                               pcent = 0.05, by_cols = FALSE, overwrite = TRUE)

    keep <- which(miss_frac(m, by_cols = FALSE) < 0.05)
    expect_equal(keep, setdiff(1:20, c(4L, 11L)))

    expect_equal(dim(out), c(length(keep), 10L))
    expect_equal(unname(as.matrix(out)), unname(m[keep, ]))

    expect_equal(rownames(out), RN(20)[keep])
    expect_equal(colnames(out), CN(10))

    out$close(); X$close()
    unlink(f)
})


# -----------------------------------------------------------------------------
# 3. filter_maf() keeps exactly the survivors' names, on both axes
# -----------------------------------------------------------------------------
# Semantics: a feature is REMOVED when maf <= maf_threshold, so only features
# with maf strictly above the threshold survive (rare-variant filtering).
test_that("filter_maf inherits the surviving dimnames (columns)", {
    f <- tempfile(fileext = ".h5")
    set.seed(104)
    m <- matrix(sample(c(0, 1, 2), 20 * 10, TRUE, prob = c(.4, .4, .2)), 20, 10)
    m[, 2] <- 0                       # monomorphic -> maf 0
    m[, 5] <- 0                       # monomorphic -> maf 0
    m[, 9] <- c(rep(0, 19), 1)        # maf 0.025 -> below 0.05
    dimnames(m) <- list(RN(20), CN(10))

    X   <- hdf5_create_matrix(f, "geno/raw", data = m, overwrite = TRUE)
    out <- filter_maf(X, out_group = "geno", out_dataset = "mafc",
                      maf_threshold = 0.05, by_cols = TRUE, overwrite = TRUE)

    keep <- which(maf_per(m, by_cols = TRUE) > 0.05)
    expect_equal(keep, c(1L, 3L, 4L, 6L, 7L, 8L, 10L))   # 2, 5 and 9 dropped

    expect_equal(dim(out), c(20L, length(keep)))
    expect_equal(unname(as.matrix(out)), unname(m[, keep]))

    expect_equal(rownames(out), RN(20))
    expect_equal(colnames(out), CN(10)[keep])

    out$close(); X$close()
    unlink(f)
})

test_that("filter_maf inherits the surviving dimnames (rows)", {
    f <- tempfile(fileext = ".h5")
    set.seed(105)
    m <- matrix(sample(c(0, 1, 2), 20 * 10, TRUE, prob = c(.4, .4, .2)), 20, 10)
    m[6,  ] <- 0                      # monomorphic row -> maf 0
    m[14, ] <- c(rep(0, 9), 1)        # maf 0.05 -> removed by the <= test
    dimnames(m) <- list(RN(20), CN(10))

    X   <- hdf5_create_matrix(f, "geno/raw", data = m, overwrite = TRUE)
    out <- filter_maf(X, out_group = "geno", out_dataset = "mafr",
                      maf_threshold = 0.05, by_cols = FALSE, overwrite = TRUE)

    keep <- which(maf_per(m, by_cols = FALSE) > 0.05)
    expect_false(6L  %in% keep)
    expect_false(14L %in% keep)

    expect_equal(dim(out), c(length(keep), 10L))
    expect_equal(unname(as.matrix(out)), unname(m[keep, ]))

    expect_equal(rownames(out), RN(20)[keep])
    expect_equal(colnames(out), CN(10))

    out$close(); X$close()
    unlink(f)
})


# -----------------------------------------------------------------------------
# 4. An input without dimnames yields a result without dimnames, no error
# -----------------------------------------------------------------------------
test_that("nothing is invented when the input has no dimnames", {
    f <- tempfile(fileext = ".h5")
    set.seed(106)
    m <- matrix(sample(c(0, 1, 2), 20 * 10, TRUE), 20, 10)
    m[1:15, 3] <- 3

    keep <- which(miss_frac(m, by_cols = TRUE) < 0.05)
    expect_equal(keep, setdiff(1:10, 3L))

    X <- hdf5_create_matrix(f, "geno/raw", data = m, overwrite = TRUE)
    expect_null(dimnames(X))

    imp <- impute_snps(X, out_group = "geno", out_dataset = "imp",
                       overwrite = TRUE)
    expect_null(dimnames(imp))

    lc <- filter_low_coverage(X, out_group = "geno", out_dataset = "lowc",
                              pcent = 0.05, by_cols = TRUE, overwrite = TRUE)
    expect_equal(dim(lc), c(20L, length(keep)))
    expect_null(dimnames(lc))

    mf <- filter_maf(X, out_group = "geno", out_dataset = "mafc",
                     maf_threshold = 0.05, by_cols = TRUE, overwrite = TRUE)
    expect_null(dimnames(mf))

    # A half-named input keeps the half it has
    colnames(X) <- CN(10)
    expect_null(rownames(X))
    lc2 <- filter_low_coverage(X, out_group = "geno", out_dataset = "lowc2",
                               pcent = 0.05, by_cols = TRUE, overwrite = TRUE)
    expect_null(rownames(lc2))
    expect_equal(colnames(lc2), CN(10)[keep])

    imp$close(); lc$close(); mf$close(); lc2$close(); X$close()
    unlink(f)
})


# -----------------------------------------------------------------------------
# 5. filter_low_coverage(): the documented semantics of pcent (diagnosis)
# -----------------------------------------------------------------------------
# This section pins down the behaviour that was first mistaken for a bug:
# filter_low_coverage(pcent = 0.05) on a matrix whose gaps are NA removed
# nothing, on either axis. It is not a bug -- missing values are the entries
# equal to 3, NAs are not counted -- and these assertions keep the semantics
# from drifting without the documentation following.
test_that("missing data means the value 3, not NA", {
    f <- tempfile(fileext = ".h5")
    set.seed(107)
    base <- matrix(sample(c(0, 1, 2), 20 * 10, TRUE), 20, 10)

    # (a) gaps stored as NA: nothing is removed, on either axis
    mna <- base
    mna[1:15, 3] <- NA
    mna[1:15, 7] <- NA
    dimnames(mna) <- list(RN(20), CN(10))
    Xna <- hdf5_create_matrix(f, "geno/na", data = mna, overwrite = TRUE)

    na_c <- filter_low_coverage(Xna, out_group = "geno", out_dataset = "na_c",
                                pcent = 0.05, by_cols = TRUE, overwrite = TRUE)
    na_r <- filter_low_coverage(Xna, out_group = "geno", out_dataset = "na_r",
                                pcent = 0.05, by_cols = FALSE, overwrite = TRUE)
    expect_equal(dim(na_c), c(20L, 10L))
    expect_equal(dim(na_r), c(20L, 10L))
    expect_equal(colnames(na_c), CN(10))
    expect_equal(rownames(na_r), RN(20))

    # (b) the same gaps recoded to 3: removed as expected
    m3 <- mna
    m3[is.na(m3)] <- 3
    X3 <- hdf5_create_matrix(f, "geno/coded", data = m3, overwrite = TRUE)
    c3 <- filter_low_coverage(X3, out_group = "geno", out_dataset = "c3",
                              pcent = 0.05, by_cols = TRUE, overwrite = TRUE)
    expect_equal(dim(c3), c(20L, 8L))
    expect_equal(colnames(c3), CN(10)[-c(3, 7)])

    # (c) the threshold is >= pcent, not > pcent: a column at exactly 10%
    #     missing is removed by pcent = 0.10
    meq <- base
    meq[1:2, 5] <- 3                         # 2/20 = 0.10
    dimnames(meq) <- list(RN(20), CN(10))
    Xeq <- hdf5_create_matrix(f, "geno/eq", data = meq, overwrite = TRUE)
    ceq <- filter_low_coverage(Xeq, out_group = "geno", out_dataset = "ceq",
                               pcent = 0.10, by_cols = TRUE, overwrite = TRUE)
    expect_equal(dim(ceq), c(20L, 9L))
    expect_equal(colnames(ceq), CN(10)[-5])

    na_c$close(); na_r$close(); c3$close(); ceq$close()
    Xna$close(); X3$close(); Xeq$close()
    unlink(f)
})

cat("\nAll omics dimname regression tests passed.\n")
