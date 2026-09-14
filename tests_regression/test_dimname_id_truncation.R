# =============================================================================
# Regression test: HDF5 dimname / string-vector id integrity (MAXSTRING)
# -----------------------------------------------------------------------------
# Row/column names and string vectors are stored in a fixed-length compound
# `char chr[MAXSTRING]`. Historically MAXSTRING = 20 truncated every id to 19
# characters, silently. The fix raised MAXSTRING to 64 (so ids up to 63 BYTES
# round-trip exactly) and made the full read block-wise.
#
# HARDENED BEHAVIOUR (this test's focus): identifiers (dimnames / rownames /
# colnames) are IDENTIFIERS -- a silently truncated id can collide with another
# and corrupt cross-layer sample matching. So an id whose length in BYTES
# exceeds MAXSTRING-1 (= 63) now raises a HARD ERROR (fail-fast) instead of
# being truncated. Length is counted in bytes, so multibyte (UTF-8) ids are
# handled correctly.
#
# Boundary cases asserted below:
#   - 63 bytes            -> preserved exactly (rows AND cols)
#   - 64 bytes            -> explicit error (no truncation)
#   - two ids sharing the first 63 bytes (len 64) -> error (would have collided)
#   - multibyte (UTF-8)   -> 64 bytes errors; 62 bytes preserved exactly
#   - write -> read -> close -> reopen -> read  -> still exact
#
# It lives in tests_regression/, git-tracked (reaches GitHub) but excluded from
# the built tarball via .Rbuildignore (^test\.*), so it does not run in
# R CMD check / ship to CRAN -- matching this repo's convention.
# Run with:  Rscript tests_regression/test_dimname_id_truncation.R
# =============================================================================

suppressMessages({
    library(BigDataStatMeth)
    library(testthat)
    library(rhdf5)
})

cat("BigDataStatMeth version:",
    as.character(packageVersion("BigDataStatMeth")), "\n")

MAX <- 63L   # MAXSTRING - 1, the maximum storable identifier length in bytes

new_file <- function() {
    f <- tempfile(fileext = ".h5")
    f
}

# -----------------------------------------------------------------------------
# 1. Exact preservation at and below the 63-byte boundary (rows AND cols),
#    unsorted, with a wide column vector (> MAXSTRBLOCK = 32) to exercise the
#    block-wise read, plus a write -> close -> reopen round-trip.
# -----------------------------------------------------------------------------
test_that("ids up to 63 bytes round-trip exactly, incl. reopen", {
    id_63a <- paste0("ROWA-", strrep("A", MAX - 5))   # exactly 63 bytes
    id_63b <- paste0("ROWB-", strrep("B", MAX - 5))   # exactly 63 bytes
    tcga   <- "TCGA-OR-A5J1-01A-11R-A29S-07"          # 28 bytes
    # 19-char-prefix pair that used to collide under the old cap (both <= 63):
    coll_1 <- paste0(strrep("A", 19), "-ONE")
    coll_2 <- paste0(strrep("A", 19), "-TWO")

    n_rows <- 120
    filler_rn <- vapply(seq_len(n_rows), function(i)
        sprintf("ROW-%05d-%s", i, paste0(sample(LETTERS, 20, TRUE),
                                         collapse = "")), character(1))
    rn <- c(id_63a, id_63b, tcga, coll_1, coll_2, "SHORT-1", filler_rn)
    rn <- rn[!duplicated(rn)]
    rn <- sample(rn)                                  # unsorted
    stopifnot(!any(duplicated(rn)), max(nchar(rn, "bytes")) == 63L)

    n_cols <- 150                                     # > MAXSTRBLOCK
    cn <- vapply(seq_len(n_cols), function(i)
        sprintf("FEATURE-%06d-%s", i, paste0(sample(LETTERS, 12, TRUE),
                                             collapse = "")), character(1))
    cn <- cn[!duplicated(cn)]

    X <- matrix(rnorm(length(rn) * length(cn)),
                nrow = length(rn), ncol = length(cn), dimnames = list(rn, cn))
    f <- new_file()
    on.exit({ suppressMessages(h5closeAll())
              if (file.exists(f)) file.remove(f) }, add = TRUE)

    Xh <- hdf5_create_matrix(f, "IN/mat", data = X, overwrite = TRUE)
    dn <- dimnames(Xh)
    expect_identical(dn[[1]], rn)
    expect_identical(dn[[2]], cn)
    expect_false(any(duplicated(dn[[1]])))
    rm(Xh); gc(); suppressMessages(h5closeAll())

    # write -> close -> REOPEN -> read : must still be exact
    Xr <- hdf5_matrix(f, "IN/mat")
    dn2 <- dimnames(Xr)
    expect_identical(dn2[[1]], rn)
    expect_identical(dn2[[2]], cn)
    rm(Xr); gc(); suppressMessages(h5closeAll())

    # independent raw read (rhdf5) confirms on-disk bytes are intact
    raw_rn <- as.character(rhdf5::h5read(f, "IN/.mat_dimnames/1")$chr)
    raw_cn <- as.character(rhdf5::h5read(f, "IN/.mat_dimnames/2")$chr)
    suppressMessages(h5closeAll())
    expect_identical(raw_rn, rn)
    expect_identical(raw_cn, cn)
})

# -----------------------------------------------------------------------------
# 2. 64-byte identifier -> hard error, on rows and on columns.
# -----------------------------------------------------------------------------
test_that("a 64-byte rowname raises a hard error (no silent truncation)", {
    too_long <- strrep("A", 64)                       # 64 bytes
    X <- matrix(rnorm(3), 1, 3, dimnames = list(too_long, c("a","b","c")))
    f <- new_file(); on.exit(if (file.exists(f)) file.remove(f), add = TRUE)
    expect_error(hdf5_create_matrix(f, "IN/w", data = X, overwrite = TRUE),
                 regexp = "identifier|63 bytes")
    suppressMessages(h5closeAll())
})

test_that("a 64-byte colname raises a hard error", {
    too_long <- strrep("C", 64)
    X <- matrix(rnorm(3), 3, 1, dimnames = list(c("r1","r2","r3"), too_long))
    f <- new_file(); on.exit(if (file.exists(f)) file.remove(f), add = TRUE)
    expect_error(hdf5_create_matrix(f, "IN/w", data = X, overwrite = TRUE),
                 regexp = "identifier|63 bytes")
    suppressMessages(h5closeAll())
})

# -----------------------------------------------------------------------------
# 3. Two ids sharing the first 63 bytes (length 64) -> error. Under silent
#    truncation these would have collapsed to the same 63-byte id.
# -----------------------------------------------------------------------------
test_that("ids sharing the first 63 bytes fail fast instead of colliding", {
    base <- strrep("Z", 63)
    a <- paste0(base, "1")                            # 64 bytes
    b <- paste0(base, "2")                            # 64 bytes, shares 63
    expect_true(substr(a, 1, 63) == substr(b, 1, 63))
    X <- matrix(rnorm(6), 2, 3, dimnames = list(c(a, b), c("a","b","c")))
    f <- new_file(); on.exit(if (file.exists(f)) file.remove(f), add = TRUE)
    expect_error(hdf5_create_matrix(f, "IN/w", data = X, overwrite = TRUE),
                 regexp = "identifier|63 bytes")
    suppressMessages(h5closeAll())
})

# -----------------------------------------------------------------------------
# 4. Multibyte (UTF-8): length is counted in BYTES, not visible characters.
# -----------------------------------------------------------------------------
test_that("multibyte id of 64 bytes errors (byte count, not char count)", {
    mb <- enc2utf8(strrep("á", 32))              # 32 x 'a-acute' = 64 bytes
    expect_equal(nchar(mb, type = "bytes"), 64L)
    expect_equal(nchar(mb, type = "chars"), 32L)      # only 32 visible chars
    X <- matrix(rnorm(3), 1, 3, dimnames = list(mb, c("a","b","c")))
    f <- new_file(); on.exit(if (file.exists(f)) file.remove(f), add = TRUE)
    expect_error(hdf5_create_matrix(f, "IN/w", data = X, overwrite = TRUE),
                 regexp = "identifier|63 bytes")
    suppressMessages(h5closeAll())
})

test_that("multibyte id of 62 bytes round-trips exactly", {
    mb <- enc2utf8(strrep("á", 31))              # 31 x 'a-acute' = 62 bytes
    expect_equal(nchar(mb, type = "bytes"), 62L)
    X <- matrix(rnorm(3), 1, 3, dimnames = list(mb, c("a","b","c")))
    f <- new_file()
    on.exit({ suppressMessages(h5closeAll())
              if (file.exists(f)) file.remove(f) }, add = TRUE)
    Xh <- hdf5_create_matrix(f, "IN/mb", data = X, overwrite = TRUE)
    got <- dimnames(Xh)[[1]]
    rm(Xh); gc(); suppressMessages(h5closeAll())
    expect_equal(nchar(got, type = "bytes"), 62L)
    expect_identical(enc2utf8(got), mb)
})

cat("\nAll dimname id-integrity regression checks passed.\n")
