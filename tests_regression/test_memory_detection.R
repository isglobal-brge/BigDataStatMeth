# =============================================================================
# Regression test: available-memory detection and the block-size budget
# -----------------------------------------------------------------------------
# getAvailableMemoryMB() (inst/include/Utilities/system-utils.hpp) is what every
# block-wise algorithm consults to decide (a) how big a processing block may be
# and (b) whether a matrix is small enough to preload entirely into RAM
# (PATH1, "fits within 20% of available RAM") instead of streaming it.
#
# THE BUG THIS PINS DOWN. The _WIN32 branch used to call R's memory.size().
# That function is defunct since R 4.2.0: it returns Inf with a warning. Inf is
# not NA and raises no exception, so it sailed past both guards, and
#
#     static_cast<size_t>(Inf * 0.6)
#
# is undefined behaviour -- in practice an astronomically large number. The
# adaptive threshold inflated accordingly and Windows effectively always chose
# the preload path, on any matrix, whatever the real memory pressure. The 4 GB
# fallback was unreachable because it only triggered on an exception.
#
# THE FIX. Windows now uses GlobalMemoryStatusEx()/ullAvailPhys (the same API
# SystemInfo.hpp already used), matching what macOS (Mach free+inactive pages)
# and Linux (/proc/meminfo MemAvailable) report; and *every* platform branch now
# funnels its result through sanitizeAvailableMemoryMB(), which rejects NaN,
# +/-Inf, non-positive and absurd values and substitutes the 4 GB fallback.
# Guarding the value rather than only the call is what makes the fallback
# reachable at all.
#
# The three assertions below are platform-independent. Section 4 is the one
# that would have failed on Windows before the fix.
#
# It lives in tests_regression/, git-tracked but excluded from the built tarball
# via .Rbuildignore (^test\.*), so it does not run in R CMD check / ship to
# CRAN -- matching this repo's convention.
# Run with:  Rscript tests_regression/test_memory_detection.R
# =============================================================================

suppressMessages({
    library(BigDataStatMeth)
    library(testthat)
})

cat("BigDataStatMeth version:",
    as.character(packageVersion("BigDataStatMeth")), "\n")
cat("Platform:", BigDataStatMeth::system_info()$os, "\n\n")

FALLBACK_MB <- 4000    # MEMORY_DETECTION_FALLBACK_MB
MAX_MB      <- 1024^3  # MEMORY_DETECTION_MAX_MB (1 PB)

sanitize <- BigDataStatMeth:::sanitize_memory_mb
budget   <- BigDataStatMeth:::get_block_memory_budget_mb
blockel  <- BigDataStatMeth:::get_optimal_block_elements

# -----------------------------------------------------------------------------
# 1. The guard rejects exactly the values that used to slip through
# -----------------------------------------------------------------------------
test_that("non-finite and implausible memory figures fall back to 4 GB", {
    # Inf is what R's defunct memory.size() returns -- the original failure.
    expect_equal(sanitize(Inf),       FALLBACK_MB)
    expect_equal(sanitize(-Inf),      FALLBACK_MB)
    expect_equal(sanitize(NaN),       FALLBACK_MB)
    expect_equal(sanitize(NA_real_),  FALLBACK_MB)

    # Non-positive: a detection that produced nothing.
    expect_equal(sanitize(0),         FALLBACK_MB)
    expect_equal(sanitize(-1),        FALLBACK_MB)

    # Absurdly large: an overflowed or garbage reading, not a machine.
    expect_equal(sanitize(MAX_MB),     FALLBACK_MB)
    expect_equal(sanitize(MAX_MB * 2), FALLBACK_MB)
    expect_equal(sanitize(1e300),      FALLBACK_MB)
})

test_that("plausible memory figures pass through unchanged", {
    for (mb in c(1, 512, 4000, 16000, 128000, 1e6))
        expect_equal(sanitize(mb), mb, info = paste("mb =", mb))
})

# -----------------------------------------------------------------------------
# 2. The live detection returns something usable on this platform
# -----------------------------------------------------------------------------
test_that("getAvailableMemoryMB() is finite, positive and plausible", {
    mb <- budget()
    expect_true(is.finite(mb))
    expect_gt(mb, 0)
    expect_lt(mb, MAX_MB)
    cat(sprintf("  detected available memory: %.0f MB (%.1f GB)\n",
                mb, mb / 1024))
})

# -----------------------------------------------------------------------------
# 3. The derived block size stays on one of the four documented tiers
# -----------------------------------------------------------------------------
test_that("getOptimalBlockElements() returns a documented tier", {
    tiers <- c(75e6, 125e6, 200e6, 500e6)   # <8GB, 8-16GB, 16-64GB, >64GB
    el <- blockel()
    expect_true(el %in% tiers)
    cat(sprintf("  block budget: %.0f elements (%.1f GB of doubles)\n",
                el, el * 8 / 1024^3))
})

# -----------------------------------------------------------------------------
# 4. The two independent detectors must now agree
# -----------------------------------------------------------------------------
# SystemInfo::getAvailableRAM_MB() (exposed as get_available_ram(), in GB) and
# getAvailableMemoryMB() measure the same quantity by the same API on each
# platform, so they should land within a factor of two of each other -- the
# slack covers only the seconds of drift between the two calls.
#
# THIS IS THE WINDOWS REGRESSION. Before the fix the two disagreed by orders of
# magnitude on Windows, because SystemInfo already used GlobalMemoryStatusEx
# while system-utils used the defunct memory.size().
test_that("the two memory detectors agree within a factor of two", {
    a <- budget()                      # system-utils.hpp  (MB)
    b <- BigDataStatMeth::get_available_ram() * 1024   # SystemInfo.hpp (GB -> MB)

    skip_if(b <= 0, "SystemInfo could not detect available RAM on this platform")

    ratio <- max(a, b) / min(a, b)
    cat(sprintf("  system-utils: %.0f MB   SystemInfo: %.0f MB   ratio: %.2f\n",
                a, b, ratio))
    expect_lt(ratio, 2)
})

# -----------------------------------------------------------------------------
# 5. The budget actually reaches the algorithms (smoke test)
# -----------------------------------------------------------------------------
# A multiplication exercises the preload-vs-stream decision that consumes the
# budget. It must produce the same answer either way; this only checks that a
# sane budget does not break the path selection.
test_that("block-wise multiplication is unaffected by the detection change", {
    f <- tempfile(fileext = ".h5")
    set.seed(11)
    A <- matrix(rnorm(200 * 60), 200, 60)
    B <- matrix(rnorm(60 * 40),   60, 40)

    hA <- hdf5_create_matrix(f, "m/A", data = A, overwrite = TRUE)
    hB <- hdf5_create_matrix(f, "m/B", data = B, overwrite = TRUE)
    hC <- hA %*% hB

    expect_equal(as.matrix(hC), A %*% B, tolerance = 1e-10,
                 ignore_attr = TRUE)

    hA$close(); hB$close(); hC$close()
    unlink(f)
})

cat("\nAll memory-detection regression tests passed.\n")
