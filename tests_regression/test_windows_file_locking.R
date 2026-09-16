# =============================================================================
# Regression test: removal of the pre-flight file-lock probe (Windows fix)
# -----------------------------------------------------------------------------
# THE BUG THIS PINS DOWN. Any downstream package performing a rapid sequence of
# close/reopen cycles on one HDF5 file failed on Windows with "HDF5 file is in
# use by another process" (deterministic, first run of a fresh session, no
# other process involved). openFile()/createFile() gated
# the real open behind lockedByOtherProcess(), a probe that tests the lock BY
# OPENING THE FILE. Windows file locks are mandatory and their release after
# H5Fclose is not immediate, so the probe could collide with the lock left by
# a handle closed an instant earlier -- including checkHDF5File()'s own RDONLY
# open, which runs immediately before the gate -- and report "locked" with no
# other process around. Two evaluations of the predicate microseconds apart
# disagreed (762-line trace, zero locked=1, yet the gate fired).
#
# THE FIX (hdf5Files.hpp, 2.0.6):
#   1. The lockedByOtherProcess() gate is gone from openFile() and createFile().
#      The real open is attempted directly; a genuine conflict surfaces as the
#      actual HDF5 error.
#   2. On _WIN32 the package now sets HDF5_USE_FILE_LOCKING=FALSE (it used to
#      force TRUE, i.e. mandatory locking, process-wide).
#   3. The probe kept for the advisory isLocked() uses a fapl with locking
#      disabled (H5Pset_file_locking(fapl, 0, 1)) so it can never take a lock.
#   4. Error paths in bdCreate_hdf5_matrix()/rcpp_hdf5_create_matrix() close
#      the borrowed dataset BEFORE deleting the hdf5File that owns the shared
#      file id (the order was reversed).
#
# Linux/macOS never showed the failure (advisory locks), so this file cannot
# reproduce the Windows symptom; what it CAN pin down on every platform is
# that the churn patterns that triggered it keep working, that no code path
# still emits the retired "in use by another process" gate error, and that the
# session stays healthy through error paths and live-handle scenarios.
# The Windows-side validation was done rebuilding the reporting downstream
# package against the patched headers and re-running its failing tests.
#
# It lives in tests_regression/, git-tracked but excluded from the built
# tarball via .Rbuildignore (^test\.*), so it does not run in R CMD check.
# Run with:  Rscript tests_regression/test_windows_file_locking.R
# =============================================================================

suppressMessages({
    library(BigDataStatMeth)
    library(testthat)
})

cat("BigDataStatMeth version:",
    as.character(packageVersion("BigDataStatMeth")), "\n\n")

GATE_MSG <- "in use by another process"

set.seed(1)
ids <- sprintf("id%02d", 1:40)
mk  <- function(p, tag) matrix(rnorm(40 * p), 40, p,
                               dimnames = list(ids, paste0(tag, "_v", seq_len(p))))

# -----------------------------------------------------------------------------
# 1. Downstream-style churn: many close-then-reopen transitions on ONE file
# -----------------------------------------------------------------------------
test_that("rapid close/reopen churn on one file never hits the retired gate", {
    f <- tempfile(fileext = ".h5")
    bdCreate_hdf5_matrix(filename = f, object = mk(12, "A"), "data", "A",
                         overwriteFile = TRUE, overwriteDataset = TRUE)

    for (i in 1:30) {
        # each cycle: open a handle, compute, close, overwrite, rewrite dimnames
        X <- hdf5_matrix(f, "data/A")
        cp <- crossprod(X, outgroup = "res", outdataset = paste0("cp", i))
        cp$close()
        X$close()
        expect_error(
            bdCreate_hdf5_matrix(filename = f, object = mk(12, "A"), "data", "A",
                                 overwriteFile = FALSE,
                                 overwriteDataset = TRUE),
            NA)
    }
    succeed()
})

# -----------------------------------------------------------------------------
# 2. Several live handles (several file ids) on the same physical file
# -----------------------------------------------------------------------------
test_that("operations run while other handles on the same file stay open", {
    f <- tempfile(fileext = ".h5")
    bdCreate_hdf5_matrix(filename = f, object = mk(9, "B"), "data", "B",
                         overwriteFile = TRUE, overwriteDataset = TRUE)
    bdCreate_hdf5_matrix(filename = f, object = mk(9, "C"), "data", "C",
                         overwriteFile = FALSE, overwriteDataset = TRUE)

    hB <- hdf5_matrix(f, "data/B")
    hC <- hdf5_matrix(f, "data/C")   # second live id on the same file

    res <- crossprod(hB, hC, outgroup = "res", outdataset = "BtC")
    expect_equal(dim(res), c(9L, 9L))

    res$close(); hC$close(); hB$close()
})

# -----------------------------------------------------------------------------
# 3. Overwriting the FILE while a live handle exists in this same process.
#    Either outcome (H5 refuses the truncate, or it succeeds) is acceptable;
#    what must hold is: no retired gate message, and the session recovers
#    fully once the handle is closed.
# -----------------------------------------------------------------------------
test_that("file overwrite with a live in-process handle fails clean or works", {
    f <- tempfile(fileext = ".h5")
    bdCreate_hdf5_matrix(filename = f, object = mk(5, "D"), "data", "D",
                         overwriteFile = TRUE, overwriteDataset = TRUE)
    h <- hdf5_matrix(f, "data/D")

    out <- tryCatch(
        { bdCreate_hdf5_matrix(filename = f, object = mk(5, "D"), "data", "D",
                               overwriteFile = TRUE,
                               overwriteDataset = TRUE); "ok" },
        error = function(e) conditionMessage(e))
    if (!identical(out, "ok")) {
        expect_false(grepl(GATE_MSG, out, fixed = TRUE),
                     info = paste("retired gate error resurfaced:", out))
    }

    h$close()
    expect_error(
        bdCreate_hdf5_matrix(filename = f, object = mk(5, "D"), "data", "D",
                             overwriteFile = TRUE, overwriteDataset = TRUE),
        NA)
    h2 <- hdf5_matrix(f, "data/D")
    expect_equal(dim(h2), c(40L, 5L))
    h2$close()
})

# -----------------------------------------------------------------------------
# 4. Error path of bdCreate_hdf5_matrix leaves the session healthy
#    (exercises the reordered checkClose_file / delete objFile cleanup)
# -----------------------------------------------------------------------------
test_that("failed create (existing dataset, no force) then normal use", {
    f <- tempfile(fileext = ".h5")
    bdCreate_hdf5_matrix(filename = f, object = mk(4, "E"), "data", "E",
                         overwriteFile = TRUE, overwriteDataset = TRUE)
    # existing dataset without overwriteDataset must fail...
    expect_error(
        bdCreate_hdf5_matrix(filename = f, object = mk(4, "E"), "data", "E",
                             overwriteFile = FALSE, overwriteDataset = FALSE))
    # ...and must not have wedged the file: same file usable right away
    X <- hdf5_matrix(f, "data/E")
    expect_equal(dim(X), c(40L, 4L))
    X$close()
})

# -----------------------------------------------------------------------------
# 5. hdf5_close_all() mid-session, then keep working on the same file
# -----------------------------------------------------------------------------
test_that("hdf5_close_all() then immediate reopen and compute", {
    f <- tempfile(fileext = ".h5")
    bdCreate_hdf5_matrix(filename = f, object = mk(7, "F"), "data", "F",
                         overwriteFile = TRUE, overwriteDataset = TRUE)
    h <- hdf5_matrix(f, "data/F")
    hdf5_close_all()
    expect_false(is_open(h))
    # a closed object must be inspectable without errors (IDE panes call these)
    expect_null(dim(h))
    expect_identical(length(h), 0L)
    expect_output(print(h), "CLOSED")

    X <- hdf5_matrix(f, "data/F")
    res <- tcrossprod(X, outgroup = "res", outdataset = "FFt")
    expect_equal(dim(res), c(40L, 40L))
    res$close(); X$close()
})

# -----------------------------------------------------------------------------
# 6. Cross-file operation with live handles on both files
# -----------------------------------------------------------------------------
test_that("cross-file multiply with live handles on both files", {
    f1 <- tempfile(fileext = ".h5")
    f2 <- tempfile(fileext = ".h5")
    bdCreate_hdf5_matrix(filename = f1, object = mk(6, "G"), "data", "G",
                         overwriteFile = TRUE, overwriteDataset = TRUE)
    bdCreate_hdf5_matrix(filename = f2, object = t(mk(6, "H")), "data", "H",
                         overwriteFile = TRUE, overwriteDataset = TRUE)

    A <- hdf5_matrix(f1, "data/G")   # 40 x 6
    B <- hdf5_matrix(f2, "data/H")   #  6 x 40
    res <- A %*% B
    expect_equal(dim(res), c(40L, 40L))
    res$close(); B$close(); A$close()
})

cat("\nAll file-locking regression checks passed.\n")
