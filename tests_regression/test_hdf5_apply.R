# =============================================================================
# Regression test: hdf5_apply(func = "CrossProd" / "tCrossProd") on inputs
# large enough to take the block-streaming path
# -----------------------------------------------------------------------------
# THE BUG. bdapply_Function_hdf5() -- the C++ worker behind hdf5_apply() --
# picks between two implementations of the (t)cross-product:
#
#     if ((dims[0] * dims[1]) > (MAXELEMSINBLOCK / 1024)) { ... }
#
# MAXELEMSINBLOCK is (2 << 29) - 1 = 1073741823, so the switch happens at
# 1073741823 / 1024 = 1048575 elements: an input with 1048576 elements or more
# (nrow_R * ncol_R, orientation irrelevant) left the in-RAM branch and entered
# the streaming branch, which calls BigDataStatMeth::crossprod() /
# ::tcrossprod() directly.
#
# That call passed a hardcoded `iblock_size = 0`. Both routines size their work
# as (N + hdf5_block - 1) / hdf5_block; an unsigned division by zero is 0 on
# AArch64 (UDIV does not trap), so the block loop ran zero iterations. The
# output dataset had already been created at the right dimensions and was left
# holding the HDF5 default fill value. Result: a correctly shaped, entirely
# zero-filled answer, no warning, no error, exit status 0, and fully
# deterministic -- including with threads = 1.
#
# THE FIX. The streaming branch now derives its block size the same way the
# crossprod()/tcrossprod() R6 bindings do, via getMaxBlockSize(..., factor 2).
# crossprod()/tcrossprod() additionally reject hdf5_block == 0 outright, so the
# silent-zeros failure mode cannot come back through another caller.
#
# WHAT THIS PINS DOWN. Sizes straddling the 1048575-element boundary, on both
# sides and immediately either side of it, for CrossProd and tCrossProd. The
# small sizes were already correct before the fix and must stay bit-comparable
# to base R; the large ones were the zero-filled ones.
#
# It lives in tests_regression/, git-tracked but excluded from the built tarball
# via .Rbuildignore (^test\.*), so it does not run in R CMD check / ship to
# CRAN -- matching this repo's convention.
# Run with:  Rscript tests_regression/test_hdf5_apply.R
# =============================================================================

suppressMessages({
    library(BigDataStatMeth)
    library(testthat)
})

cat("BigDataStatMeth version:",
    as.character(packageVersion("BigDataStatMeth")), "\n\n")

# The element count at which bdapply_Function_hdf5() switches implementations.
SWITCH_AT <- 1048575L   # (2 << 29 - 1) %/% 1024

TOL <- 1e-10

# Run hdf5_apply() on one freshly written matrix and hand back the result.
# Every call gets its own temp file so nothing carries over between checks.
apply_one <- function(m, func, threads = 1L) {
    f <- tempfile(fileext = ".h5")
    on.exit({ hdf5_close_all(); unlink(f) }, add = TRUE)

    hdf5_create_matrix(f, "in/M", data = m, overwrite = TRUE)
    hdf5_close_all()

    hdf5_apply(f, group = "in", datasets = "M", func = func,
               outgroup = "out", overwrite = TRUE, threads = threads)

    h   <- hdf5_matrix(f, "out/M")
    got <- as.matrix(h)
    close(h)
    got
}

# One size, one function: compare against base R and report where the size sits
# relative to the branch boundary.
check_size <- function(nr, nc, func, threads = 1L) {
    n_elems <- as.numeric(nr) * nc
    branch  <- if (n_elems > SWITCH_AT) "streaming" else "in-RAM"

    set.seed(11)
    m   <- matrix(rnorm(n_elems), nr, nc)
    ref <- if (func == "CrossProd") base::crossprod(m) else base::tcrossprod(m)
    got <- apply_one(m, func, threads)

    label <- sprintf("%s %dx%d (%.0f elems, %s)", func, nr, nc, n_elems, branch)

    expect_equal(dim(got), dim(ref), info = label)
    # The pre-fix failure was an all-zero matrix, so state that separately from
    # the numeric comparison: it makes a regression unmistakable in the output.
    expect_gt(sum(got != 0), 0L)
    expect_equal(got, ref, tolerance = TOL, ignore_attr = TRUE, info = label)

    cat(sprintf("  %-46s max |diff| = %.3g\n", label, max(abs(got - ref))))
    invisible(NULL)
}

# -----------------------------------------------------------------------------
# 1. CrossProd across the branch boundary
# -----------------------------------------------------------------------------
# 2000x524 = 1048000 elems -> in-RAM (last size below the switch)
# 2000x525 = 1050000 elems -> streaming (first size above it)
# 1024x1024 = 1048576 elems -> streaming by a single element
test_that("CrossProd is correct on both sides of the block-size boundary", {
    check_size(200L,  50L,  "CrossProd")    #      10000  in-RAM
    check_size(1500L, 500L, "CrossProd")    #     750000  in-RAM
    check_size(3000L, 300L, "CrossProd")    #     900000  in-RAM
    check_size(2000L, 524L, "CrossProd")    #    1048000  in-RAM, just below
    check_size(1024L, 1024L, "CrossProd")   #    1048576  streaming, just above
    check_size(2000L, 525L, "CrossProd")    #    1050000  streaming
    check_size(2000L, 800L, "CrossProd")    #    1600000  streaming (the report)
})

# -----------------------------------------------------------------------------
# 2. tCrossProd across the same boundary
# -----------------------------------------------------------------------------
# tCrossProd squares the ROW count, so the outputs stay a sane size only for
# short-and-wide inputs; the element counts still straddle the switch.
test_that("tCrossProd is correct on both sides of the block-size boundary", {
    check_size(50L,  200L,  "tCrossProd")   #      10000  in-RAM
    check_size(500L, 1500L, "tCrossProd")   #     750000  in-RAM
    check_size(524L, 2000L, "tCrossProd")   #    1048000  in-RAM, just below
    check_size(1024L, 1024L, "tCrossProd")  #    1048576  streaming, just above
    check_size(525L, 2000L, "tCrossProd")   #    1050000  streaming
    check_size(800L, 2000L, "tCrossProd")   #    1600000  streaming
})

# -----------------------------------------------------------------------------
# 3. The streaming branch is thread-count independent
# -----------------------------------------------------------------------------
# The fixed branch always runs with bparal = TRUE, so the block loop is an
# OpenMP parallel for. Neither the requested thread count nor the answer may
# depend on it.
test_that("the streaming branch gives the same answer at any thread count", {
    set.seed(11)
    m   <- matrix(rnorm(2000 * 800), 2000L, 800L)
    ref <- base::crossprod(m)

    for (th in c(1L, 4L)) {
        got <- apply_one(m, "CrossProd", threads = th)
        expect_equal(got, ref, tolerance = TOL, ignore_attr = TRUE,
                     info = sprintf("threads = %d", th))
        cat(sprintf("  CrossProd 2000x800 threads=%-3d              max |diff| = %.3g\n",
                    th, max(abs(got - ref))))
    }
})

# -----------------------------------------------------------------------------
# 4. blockmult over the boundary, and the S3 apply_function() route
# -----------------------------------------------------------------------------
# blockmult takes the same `> MAXELEMSINBLOCK / 1024` decision but hands the
# large case to multiplication(), which sizes its own blocks -- it was never
# affected. Pinned here so the two branches stay in step.
test_that("blockmult over the boundary is unaffected", {
    set.seed(11)
    A <- matrix(rnorm(2000 * 800), 2000L, 800L)   # 1600000 elems -> streaming
    B <- matrix(rnorm(800 * 20), 800L, 20L)

    f <- tempfile(fileext = ".h5")
    on.exit({ hdf5_close_all(); unlink(f) }, add = TRUE)

    hdf5_create_matrix(f, "in/A", data = A, overwrite = TRUE)
    hdf5_create_matrix(f, "in/B", data = B, overwrite = TRUE)
    hdf5_close_all()

    hdf5_apply(f, group = "in", datasets = "A", func = "blockmult",
               outgroup = "out", b_group = "in", b_datasets = "B",
               overwrite = TRUE, threads = 1L)

    h   <- hdf5_matrix(f, "out/A_B")
    got <- as.matrix(h)
    close(h)

    expect_equal(got, A %*% B, tolerance = TOL, ignore_attr = TRUE)
    cat(sprintf("  blockmult 2000x800 %%*%% 800x20                  max |diff| = %.3g\n",
                max(abs(got - A %*% B))))
})

# apply_function() (the S3 method on an open HDF5Matrix) reaches a *different*
# C++ worker -- rcpp_hdf5dataset_apply_function() -> RcppApplyFunctionHdf5() --
# which has no size switch at all and always computes in RAM. It never produced
# zeros; the check is here because the two entry points are documented as
# interchangeable and must agree.
test_that("apply_function() agrees with hdf5_apply() above the boundary", {
    set.seed(11)
    m <- matrix(rnorm(2000 * 800), 2000L, 800L)

    f <- tempfile(fileext = ".h5")
    on.exit({ hdf5_close_all(); unlink(f) }, add = TRUE)

    h <- hdf5_create_matrix(f, "in/M", data = m, overwrite = TRUE)
    suppressMessages(apply_function(h, func = "CrossProd", out_group = "s3out",
                                    overwrite = TRUE, threads = 1L))
    close(h)

    # apply_function() renames its outputs to <func>_<dataset>.
    res <- list_datasets(f, group = "s3out")
    expect_gt(length(res), 0L)

    hs3 <- hdf5_matrix(f, paste0("s3out/", basename(res[1])))
    got <- as.matrix(hs3)
    close(hs3)

    expect_equal(got, base::crossprod(m), tolerance = TOL, ignore_attr = TRUE)
    cat(sprintf("  apply_function() CrossProd 2000x800           max |diff| = %.3g\n",
                max(abs(got - base::crossprod(m)))))
})

hdf5_close_all()

cat("\nAll hdf5_apply regression tests passed.\n")
