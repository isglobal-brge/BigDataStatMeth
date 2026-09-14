# =============================================================================
# Regression test: zero-variance guard for svd()/prcomp() with scale = TRUE
# -----------------------------------------------------------------------------
# THE BUG THIS PINS DOWN. svd(X, center = TRUE, scale = TRUE) on a matrix with
# a constant column returned ALL singular values as 0 -- no error, no warning,
# no NaN in sight. Scaling divides by the column standard deviation, which is
# exactly 0 for a constant column; the resulting NaN spreads through the Gram
# matrix, LAPACK returns nothing usable, and the caller gets a plausible-looking
# vector of zeros. There was no sd == 0 check at any layer: the divisions in
# matrixNormalization.hpp were unconditional and the standard deviations from
# matrixSdMean.hpp were never inspected.
#
# THE FIX. A guard now runs after the standard deviations have been computed
# (they are already computed block-wise, so it costs nothing extra and stays
# out-of-core) and before any division, at three points in the SVD/PCA flow:
#
#   · RcppbdSVD_hdf5()                          exact / in-RAM LAPACK regime
#   · First_level_SvdBlock_decomposition_hdf5()  approximate / block regime
#   · RcppTypifyNormalizeHdf5()                 the PCA normalization step
#
# It throws "cannot rescale a constant/zero column to unit variance: ...",
# following the precedent of stats::prcomp(x, scale. = TRUE), and names how
# many columns are affected plus a few of their positions (1-based, in the
# orientation R sees, not the transposed HDF5 one).
#
# WHAT IS DELIBERATELY NOT GUARDED. scale.HDF5Matrix() keeps mimicking
# base::scale(), which produces NaN for a constant column without complaining.
# That is why the guard lives in the SVD/PCA flow and not inside the shared
# normalization routine. Section 5 pins that down.
#
# It lives in tests_regression/, git-tracked but excluded from the built tarball
# via .Rbuildignore (^test\.*), so it does not run in R CMD check / ship to
# CRAN -- matching this repo's convention.
# Run with:  Rscript tests_regression/test_zero_variance_scale.R
# =============================================================================

suppressMessages({
    library(BigDataStatMeth)
    library(testthat)
})

cat("BigDataStatMeth version:",
    as.character(packageVersion("BigDataStatMeth")), "\n\n")

# -----------------------------------------------------------------------------
# Fixtures: the same matrix with and without a constant column
# -----------------------------------------------------------------------------
set.seed(42)
NR <- 60L
NC <- 8L

x_ok  <- matrix(rnorm(NR * NC), NR, NC)   # every column has positive variance
x_bad <- x_ok
x_bad[, 3] <- 7                            # one constant column
x_bad2 <- x_ok
x_bad2[, 2] <- 0                           # zero column
x_bad2[, 6] <- -1.5                        # ... and a second constant one

# The message fragment every guard must produce, whatever the code path.
FRAGMENT <- "constant/zero column"

new_h5 <- function(data, path = "m/A") {
    f <- tempfile(fileext = ".h5")
    h <- hdf5_create_matrix(f, path, data = data, overwrite = TRUE)
    list(file = f, mat = h)
}

# -----------------------------------------------------------------------------
# 1. svd(scale = TRUE) on a constant column must error, not degenerate
# -----------------------------------------------------------------------------
test_that("svd() rejects a constant column when scale = TRUE", {
    h <- new_h5(x_bad)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    expect_error(svd(h$mat, center = TRUE, scale = TRUE), FRAGMENT)
})

test_that("svd(center = FALSE, scale = TRUE) rejects a zero column", {
    # With center = FALSE the divisor is the root mean square, not the sd --
    # base::scale() behaves the same way -- so a constant NON-zero column is
    # not a division by zero there, but an all-zero column always is.  The
    # guard mirrors that formula instead of assuming the centred one.
    h <- new_h5(x_bad2)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    expect_error(svd(h$mat, center = FALSE, scale = TRUE), FRAGMENT)
})

test_that("the message names the count and the R-side column positions", {
    h <- new_h5(x_bad2)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    msg <- tryCatch({
        svd(h$mat, center = TRUE, scale = TRUE)
        NA_character_
    }, error = conditionMessage)

    expect_true(grepl(FRAGMENT, msg, fixed = TRUE))
    expect_true(grepl("2 columns", msg, fixed = TRUE))
    # 1-based, R orientation: columns 2 and 6, NOT the HDF5 transposed indices.
    expect_true(grepl("columns 2, 6", msg, fixed = TRUE))
    cat("  reported: ", msg, "\n", sep = "")
})

# -----------------------------------------------------------------------------
# 2. The guard covers the approximate (block) regime too
# -----------------------------------------------------------------------------
# method = "blocks" takes a completely different route: the standard deviations
# come from get_HDF5_mean_sd_by_column() inside
# First_level_SvdBlock_decomposition_hdf5(), not from an in-RAM Eigen matrix.
test_that("the block SVD path rejects a constant column as well", {
    h <- new_h5(x_bad)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    expect_error(
        svd(h$mat, center = TRUE, scale = TRUE, method = "blocks"),
        FRAGMENT)
})

# -----------------------------------------------------------------------------
# 3. scale = FALSE is unaffected: a constant column is legitimate there
# -----------------------------------------------------------------------------
test_that("svd(scale = FALSE) still works on a matrix with a constant column", {
    h <- new_h5(x_bad)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    res <- suppressMessages(svd(h$mat, center = TRUE, scale = FALSE))

    ref <- base::svd(scale(x_bad, center = TRUE, scale = FALSE))$d
    expect_equal(res$d, ref, tolerance = 1e-8, ignore_attr = TRUE)

    # Centering alone turns the constant column into zeros, so the last
    # singular value is (numerically) zero -- and that is a valid answer.
    expect_true(all(is.finite(res$d)))
    expect_lt(res$d[NC], 1e-8 * res$d[1])
})

# -----------------------------------------------------------------------------
# 4. Well-conditioned input with scale = TRUE still gives the right answer
# -----------------------------------------------------------------------------
# The guard must not fire here, and must not perturb the result.
test_that("svd(scale = TRUE) matches base::svd(scale(x)) when no column is constant", {
    h <- new_h5(x_ok)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    res <- suppressMessages(svd(h$mat, center = TRUE, scale = TRUE))
    ref <- base::svd(scale(x_ok, center = TRUE, scale = TRUE))$d

    expect_equal(res$d, ref, tolerance = 1e-8, ignore_attr = TRUE)
    cat(sprintf("  max |d - d_ref| = %.3g\n", max(abs(res$d - ref))))
})

test_that("the block path also agrees with base::svd(scale(x)) at full rank", {
    h <- new_h5(x_ok)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    res <- suppressMessages(
        svd(h$mat, center = TRUE, scale = TRUE, method = "blocks"))
    ref <- base::svd(scale(x_ok, center = TRUE, scale = TRUE))$d

    # Looser than the exact path: the hierarchical merge is an approximation
    # (about 1e-15 relative at full rank, per ?svd.HDF5Matrix).
    expect_equal(res$d, ref, tolerance = 1e-6, ignore_attr = TRUE)
})

# -----------------------------------------------------------------------------
# 5. scale() keeps the base R semantics -- the guard must NOT reach it
# -----------------------------------------------------------------------------
# base::scale() returns NaN for a constant column without erroring, and
# scale.HDF5Matrix() must keep doing exactly that.  If this test starts failing
# with a "constant/zero column" error, the guard has leaked into the shared
# normalization routine.
test_that("scale() on a constant column still returns NaN, not an error", {
    h <- new_h5(x_bad)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    sc <- suppressMessages(scale(h$mat, center = TRUE, scale = TRUE))
    m  <- as.matrix(sc)

    expect_true(all(is.nan(m[, 3])))          # the constant column
    expect_true(all(is.finite(m[, -3])))      # everything else survived
})

# -----------------------------------------------------------------------------
# 6. prcomp() -- same three assertions through the PCA pipeline
# -----------------------------------------------------------------------------
# prcomp() normalizes first (RcppTypifyNormalizeHdf5) and only then runs the
# SVD with center/scale already applied, so it needs its own guard and its own
# test: the SVD-side guards above are never reached from here.
test_that("prcomp() rejects a constant column when scale. = TRUE", {
    h <- new_h5(x_bad)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    expect_error(prcomp(h$mat, center = TRUE, scale. = TRUE), FRAGMENT)
})

test_that("prcomp(scale. = FALSE) still works on a matrix with a constant column", {
    h <- new_h5(x_bad)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    pca <- suppressMessages(prcomp(h$mat, center = TRUE, scale. = FALSE))

    expect_true(all(is.finite(pca$sdev)))
    expect_gt(pca$sdev[1], 0)
})

test_that("prcomp(scale. = TRUE) is correct when no column is constant", {
    h <- new_h5(x_ok)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    pca <- suppressMessages(prcomp(h$mat, center = TRUE, scale. = TRUE))
    ref <- stats::prcomp(x_ok, center = TRUE, scale. = TRUE)$sdev

    expect_true(all(is.finite(pca$sdev)))
    expect_equal(pca$sdev, ref, tolerance = 1e-8, ignore_attr = TRUE)
    cat(sprintf("  max |sdev - sdev_ref| = %.3g\n",
                max(abs(pca$sdev - ref))))
})

hdf5_close_all()

cat("\nAll zero-variance guard regression tests passed.\n")
