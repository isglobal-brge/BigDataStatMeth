# =============================================================================
# Regression test: omics filters that remove every feature
# -----------------------------------------------------------------------------
# THE BUG THIS PINS DOWN. filter_low_coverage() and filter_maf() delegate to
# header-level routines that create the output dataset lazily: the dataset is
# only created when the first block carrying at least one SURVIVING element is
# written. When the criterion removes every column (or every row) no block ever
# qualifies, so nothing is created -- Rcpp_Remove_Low_Data_hdf5() emits the
# warning "All data removed - please adjust pcent parameter or review data" and
# Rcpp_Remove_MAF_hdf5() returns silently. The R6 layer ignored that and went
# straight to hdf5_matrix() on the location the wrapper had reported, which
# does not exist, so the user got the internal, unactionable
#
#     "please create Dataset before proceed"
#
# instead of being told that their threshold had wiped out the whole matrix.
#
# THE FIX. .hdf5_filter_output_missing() (R/HDF5Matrix_op_omics.R) runs between
# the C++ call and the hdf5_matrix() call. It treats an empty `kept` vector (no
# element survived) or the absence of the dataset from its group as conclusive,
# and the two R6 methods then stop() with a message that names what was removed
# and which parameter to adjust. The C++ side is untouched; the warning may
# still be emitted first, but the visible behaviour is the error.
#
# WHAT IS CHECKED. Both filters, both axes (by_cols = TRUE / FALSE): the call
# must raise a condition whose message contains "removed all" and must NOT
# contain "please create Dataset". A non-degenerate call of each filter is kept
# alongside as a control, so a guard that fired unconditionally would fail too.
#
# It lives in tests_regression/, git-tracked (reaches GitHub) but excluded from
# the built tarball via .Rbuildignore (^test\.*), so it does not run in
# R CMD check / ship to CRAN -- matching this repo's convention.
# Run with:  Rscript tests_regression/test_omics_empty_result.R
# =============================================================================

suppressMessages({
    library(BigDataStatMeth)
    library(testthat)
})

cat("BigDataStatMeth version:",
    as.character(packageVersion("BigDataStatMeth")), "\n\n")

# Run `expr`, swallowing any warning, and return the error message or NA.
err_msg <- function(expr) {
    withCallingHandlers(
        tryCatch({ force(expr); NA_character_ },
                 error = function(e) conditionMessage(e)),
        warning = function(w) invokeRestart("muffleWarning"))
}

expect_removed_all <- function(msg, label) {
    expect_true(!is.na(msg),
                info = paste0(label, ": expected an error, none was raised"))
    expect_true(grepl("removed all", msg, fixed = TRUE),
                info = paste0(label, ": message lacks 'removed all': ", msg))
    expect_false(grepl("please create Dataset", msg, fixed = TRUE),
                 info = paste0(label, ": still the internal message: ", msg))
}


# -----------------------------------------------------------------------------
# 1. filter_low_coverage(): every column / every row is missing
# -----------------------------------------------------------------------------
test_that("filter_low_coverage errors when it removes every feature", {
    f <- tempfile(fileext = ".h5")
    set.seed(201)

    # Fully missing matrix: the missing code of the omics layer is the value 3.
    # Every column and every row is at 100% missing, so any pcent <= 1 removes
    # the lot along either axis.
    allmiss <- matrix(3, 20, 10)
    Xm <- hdf5_create_matrix(f, "geno/allmiss", data = allmiss, overwrite = TRUE)

    expect_removed_all(
        err_msg(filter_low_coverage(Xm, out_group = "geno",
                                    out_dataset = "lc_cols",
                                    pcent = 0.5, by_cols = TRUE,
                                    overwrite = TRUE)),
        "filter_low_coverage / by_cols = TRUE")

    expect_removed_all(
        err_msg(filter_low_coverage(Xm, out_group = "geno",
                                    out_dataset = "lc_rows",
                                    pcent = 0.5, by_cols = FALSE,
                                    overwrite = TRUE)),
        "filter_low_coverage / by_cols = FALSE")

    # The message points at the parameter to adjust and at the right axis
    expect_true(grepl("pcent",
                      err_msg(filter_low_coverage(Xm, out_group = "geno",
                                                  out_dataset = "lc_p",
                                                  pcent = 0.5, by_cols = TRUE,
                                                  overwrite = TRUE)),
                      fixed = TRUE))
    expect_true(grepl("removed all columns",
                      err_msg(filter_low_coverage(Xm, out_group = "geno",
                                                  out_dataset = "lc_c2",
                                                  pcent = 0.5, by_cols = TRUE,
                                                  overwrite = TRUE)),
                      fixed = TRUE))
    expect_true(grepl("removed all rows",
                      err_msg(filter_low_coverage(Xm, out_group = "geno",
                                                  out_dataset = "lc_r2",
                                                  pcent = 0.5, by_cols = FALSE,
                                                  overwrite = TRUE)),
                      fixed = TRUE))

    # The input object survives the failed call and stays usable
    expect_true(Xm$is_valid())
    expect_equal(dim(Xm), c(20L, 10L))

    Xm$close()
    unlink(f)
})


# -----------------------------------------------------------------------------
# 2. filter_maf(): every feature is monomorphic (MAF = 0)
# -----------------------------------------------------------------------------
test_that("filter_maf errors when it removes every feature", {
    f <- tempfile(fileext = ".h5")

    # All-zero genotypes: maf = n0/n = 1 -> 1 - 1 = 0 for every row and column,
    # so maf <= maf_threshold removes everything along either axis.
    mono <- matrix(0, 20, 10)
    Xz <- hdf5_create_matrix(f, "geno/mono", data = mono, overwrite = TRUE)

    expect_removed_all(
        err_msg(filter_maf(Xz, out_group = "geno", out_dataset = "maf_rows",
                           maf_threshold = 0.05, by_cols = FALSE,
                           overwrite = TRUE)),
        "filter_maf / by_cols = FALSE")

    expect_removed_all(
        err_msg(filter_maf(Xz, out_group = "geno", out_dataset = "maf_cols",
                           maf_threshold = 0.05, by_cols = TRUE,
                           overwrite = TRUE)),
        "filter_maf / by_cols = TRUE")

    expect_true(grepl("maf_threshold",
                      err_msg(filter_maf(Xz, out_group = "geno",
                                         out_dataset = "maf_p",
                                         maf_threshold = 0.05, by_cols = FALSE,
                                         overwrite = TRUE)),
                      fixed = TRUE))

    expect_true(Xz$is_valid())
    expect_equal(dim(Xz), c(20L, 10L))

    Xz$close()
    unlink(f)
})


# -----------------------------------------------------------------------------
# 3. Control: a filter that keeps something must NOT trip the new guard
# -----------------------------------------------------------------------------
test_that("the empty-result guard does not fire on a normal filter", {
    f <- tempfile(fileext = ".h5")
    set.seed(202)

    m <- matrix(sample(c(0, 1, 2), 20 * 10, TRUE), 20, 10)
    m[1:15, 3] <- 3                     # 75% missing -> this column goes
    X <- hdf5_create_matrix(f, "geno/raw", data = m, overwrite = TRUE)

    lc <- filter_low_coverage(X, out_group = "geno", out_dataset = "lc_ok",
                              pcent = 0.05, by_cols = TRUE, overwrite = TRUE)
    expect_s3_class(lc, "HDF5Matrix")
    expect_equal(dim(lc), c(20L, 9L))
    lc$close()

    mf <- filter_maf(X, out_group = "geno", out_dataset = "maf_ok",
                     maf_threshold = 0.05, by_cols = FALSE, overwrite = TRUE)
    expect_s3_class(mf, "HDF5Matrix")
    expect_equal(nrow(mf), 20L)
    mf$close()

    # impute_snps() changes no dimension and always writes its output, so it is
    # not guarded; assert that it still returns a usable matrix here.
    im <- impute_snps(X, out_group = "geno", out_dataset = "imp_ok",
                      overwrite = TRUE)
    expect_s3_class(im, "HDF5Matrix")
    expect_equal(dim(im), c(20L, 10L))
    im$close()

    X$close()
    unlink(f)
})

cat("\nAll omics empty-result regression tests passed.\n")
