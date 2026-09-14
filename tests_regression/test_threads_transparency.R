# =============================================================================
# Regression test: the `threads` argument of svd()/prcomp() is honest and effective
# -----------------------------------------------------------------------------
# TWO BUGS THIS PINS DOWN.
#
# (1) FLAT SCALING. The matrix products inside the block SVD were called with
#     bparal = R_NilValue hardcoded (matrixSvd.hpp). multiplication() only
#     enters its OpenMP streaming path when bparal is explicitly TRUE
#     (multiplication.hpp); with NULL it decides on its own whether to preload
#     both operands and let the BLAS backend thread, and never consults the
#     thread count. So `threads = N` travelled all the way down the call stack
#     and changed nothing.
#
#     THE FIX. The caller's intent is now propagated: bparal = TRUE is passed
#     to those products, but ONLY when a thread count was actually requested.
#     With threads = NULL / -1 the path selection is bit-for-bit what it was,
#     so nobody who does not ask for threads sees a different number or a
#     different runtime.
#
# (2) SILENT CLAMP. get_number_threads() (openme-utils.hpp) discards a
#     requested value above the system ceiling, and the default ceiling is only
#     50% of the detected CPUs (R_DATATABLE_NUM_PROCS_PERCENT). Asking for 64
#     threads on a 4-core laptop looked like it had worked.
#
#     THE FIX. The R6 layer compares the request against the ceiling -- read
#     from the same C++ function the algorithms use, exposed as the internal
#     rcpp_effective_threads() -- and warns once when fewer threads can be used.
#
# WHAT MUST NOT CHANGE. Threading may not alter the answer. Sections 3 and 4
# check that the results agree to 1e-10 across thread counts, on each of the
# three product call sites that were touched.
#
# It lives in tests_regression/, git-tracked but excluded from the built tarball
# via .Rbuildignore (^test\.*), so it does not run in R CMD check / ship to
# CRAN -- matching this repo's convention.
# Run with:  Rscript tests_regression/test_threads_transparency.R
# =============================================================================

suppressMessages({
    library(BigDataStatMeth)
    library(testthat)
})

eff_threads <- BigDataStatMeth:::rcpp_effective_threads

cat("BigDataStatMeth version:",
    as.character(packageVersion("BigDataStatMeth")), "\n")
cat("detected cores:", BigDataStatMeth::get_cpu_cores(),
    " usable threads:", eff_threads(), "\n\n")

ABSURD <- 9999L   # more threads than any machine this will run on

# Collect the warnings raised by `expr` without letting them abort anything,
# so a test can assert on exactly the warning it cares about and ignore the
# rest.  Written by hand rather than with expect_no_warning() so the file
# works on older testthat releases too.
warnings_of <- function(expr) {
    w <- character(0)
    withCallingHandlers(
        suppressMessages(force(expr)),
        warning = function(cond) {
            w <<- c(w, conditionMessage(cond))
            invokeRestart("muffleWarning")
        })
    w
}

# One fresh input dataset per decomposition: svd()/prcomp() write their results
# next to the input, so reusing a dataset would need overwrite = TRUE and would
# invalidate handles we still hold.
seq_file <- new.env(parent = emptyenv())
seq_file$n <- 0L
new_h5 <- function(data) {
    seq_file$n <- seq_file$n + 1L
    f <- tempfile(fileext = ".h5")
    h <- hdf5_create_matrix(f, paste0("m/A", seq_file$n),
                            data = data, overwrite = TRUE)
    list(file = f, mat = h)
}

set.seed(7)
x_tall <- matrix(rnorm(240L * 24L), 240L, 24L)   # nrow >= ncol
x_wide <- matrix(rnorm(24L * 240L), 24L, 240L)   # nrow <  ncol

# -----------------------------------------------------------------------------
# 1. The ceiling probe itself
# -----------------------------------------------------------------------------
test_that("rcpp_effective_threads() reports the ceiling the algorithms use", {
    auto <- eff_threads()
    expect_true(is.finite(auto))
    expect_gte(auto, 1L)

    expect_equal(eff_threads(1L), 1L)           # a request of 1 is always met
    expect_equal(eff_threads(ABSURD), auto)     # above the ceiling -> clamped
    expect_lt(eff_threads(ABSURD), ABSURD)
})

# -----------------------------------------------------------------------------
# 2. An impossible request warns; an unspecified one stays silent
# -----------------------------------------------------------------------------
test_that("svd() warns when it cannot honour the requested thread count", {
    skip_if(eff_threads(ABSURD) >= ABSURD,
            "this machine can actually provide 9999 threads")

    h <- new_h5(x_tall)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    w   <- character(0)
    res <- withCallingHandlers(
        suppressMessages(svd(h$mat, threads = ABSURD)),
        warning = function(cond) {
            w <<- c(w, conditionMessage(cond))
            invokeRestart("muffleWarning")
        })

    expect_true(any(grepl("requested threads", w)))
    expect_true(any(grepl(sprintf("requested threads = %d", ABSURD), w)))
    cat("  warning: ", grep("requested threads", w, value = TRUE)[1], "\n",
        sep = "")

    # ... and the impossible request must not damage the result.
    ref <- base::svd(scale(x_tall, center = TRUE, scale = TRUE))$d
    expect_equal(res$d, ref, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("threads = NULL and threads = -1 produce no thread warning", {
    h1 <- new_h5(x_tall)
    h2 <- new_h5(x_tall)
    on.exit({ hdf5_close_all(); unlink(c(h1$file, h2$file)) }, add = TRUE)

    w_null <- warnings_of(svd(h1$mat, threads = NULL))
    w_auto <- warnings_of(svd(h2$mat, threads = -1L))

    expect_length(grep("requested threads", w_null), 0L)
    expect_length(grep("requested threads", w_auto), 0L)
})

test_that("prcomp() reports the same capped thread count", {
    skip_if(eff_threads(ABSURD) >= ABSURD,
            "this machine can actually provide 9999 threads")

    h <- new_h5(x_tall)
    on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)

    w <- warnings_of(prcomp(h$mat, threads = ABSURD))
    expect_true(any(grepl("requested threads", w)))
})

# -----------------------------------------------------------------------------
# 3. Threading does not change the answer -- exact (in-RAM) regime
# -----------------------------------------------------------------------------
test_that("svd() gives identical results with 2 threads and with the default", {
    h_ser <- new_h5(x_tall)
    h_par <- new_h5(x_tall)
    on.exit({ hdf5_close_all(); unlink(c(h_ser$file, h_par$file)) }, add = TRUE)

    d_ser <- suppressMessages(svd(h_ser$mat, threads = -1L))$d
    d_par <- suppressMessages(svd(h_par$mat, threads = 2L))$d

    expect_equal(d_par, d_ser, tolerance = 1e-10, ignore_attr = TRUE)
})

# -----------------------------------------------------------------------------
# 4. Threading does not change the answer -- block regime, the products touched
# -----------------------------------------------------------------------------
# bparal = TRUE is now propagated to the multiplication() call sites of
# RcppbdSVD_hdf5_Block(), each reached by a different combination:
#
#   · raw dsA            method = "blocks", center = FALSE, scale = FALSE
#   · NORMALIZED_T       prcomp(method = "blocks")
#   · normalized matrix  method = "blocks" with nrow < ncol      (section 5)
#
# IMPORTANT: the singular values alone would NOT detect a difference here.
# They come from the LAPACK SVD of the joined block matrix, which the products
# never touch; the product produces the OTHER factor.  So each check compares
# the singular values together with the factor that the multiplication writes
# (|u| for the block SVD, |x| for the block PCA), taking absolute values
# because the sign of a singular vector is arbitrary.
compare_threaded <- function(data, label, run) {
    h_ser <- new_h5(data)
    h_par <- new_h5(data)
    on.exit({ hdf5_close_all(); unlink(c(h_ser$file, h_par$file)) }, add = TRUE)

    v_ser <- suppressMessages(run(h_ser$mat, -1L))
    v_par <- suppressMessages(run(h_par$mat,  2L))

    expect_equal(length(v_par), length(v_ser), info = label)
    expect_equal(v_par, v_ser, tolerance = 1e-10, ignore_attr = TRUE,
                 info = label)
    cat(sprintf("  %-22s %6d values   max |diff| = %.3g\n", label,
                length(v_par), max(abs(v_par - v_ser))))
    invisible(NULL)
}

test_that("block SVD on the raw matrix is thread-count independent", {
    compare_threaded(x_tall, "blocks, raw", function(m, th) {
        res <- svd(m, center = FALSE, scale = FALSE, method = "blocks",
                   threads = th)
        c(res$d, abs(as.vector(as.matrix(res$u))))
    })
})

test_that("block PCA is thread-count independent", {
    compare_threaded(x_tall, "blocks, PCA", function(m, th) {
        pca <- prcomp(m, center = TRUE, scale. = TRUE, method = "blocks",
                      threads = th)
        c(pca$sdev, abs(as.vector(as.matrix(pca$x))))
    })
})

# -----------------------------------------------------------------------------
# 5. Same check for the wide-matrix branch, skipped if that path is unavailable
# -----------------------------------------------------------------------------
# The third product call site is only reached with center/scale = TRUE and
# nrow < ncol.  That branch of the block algorithm carries an unresolved review
# note of its own in matrixSvdBlock.hpp (the dimensions passed to
# createDataset() for the normalized matrix), unrelated to threading, so a hard
# failure there is reported as a skip rather than as a threading regression.
test_that("block SVD on the normalized matrix is thread-count independent", {
    run <- function(m, th) {
        res <- svd(m, center = TRUE, scale = TRUE, method = "blocks",
                   threads = th)
        c(res$d, abs(as.vector(as.matrix(res$u))))
    }

    probe <- tryCatch({
        h <- new_h5(x_wide)
        on.exit({ hdf5_close_all(); unlink(h$file) }, add = TRUE)
        suppressMessages(run(h$mat, -1L))
        TRUE
    }, error = function(e) conditionMessage(e))

    if (!isTRUE(probe))
        skip(paste("wide block path unavailable:", probe))

    compare_threaded(x_wide, "blocks, normalized", run)
})

hdf5_close_all()

cat("\nAll thread-transparency regression tests passed.\n")
