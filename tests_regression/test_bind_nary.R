## test_bind_nary.R
## -----------------------------------------------------------------------------
## Regression test: cbind() / rbind() on three or more HDF5Matrix objects.
##
## Before the fix, cbind()/rbind() with N > 2 arguments chained binary binds
## and wrote intermediate datasets ("cbind_N_<uid>_s1", "cbind_N_s2", ...) to
## the output group. The "_s2"+ names had no unique suffix and were written
## with overwrite = TRUE, so concurrent or repeated executions could overwrite
## each other's intermediates, and all intermediates were left in the file.
## The result is now built in a single pass, with no intermediate datasets.
##
## Checks: exact values against base cbind()/rbind() (3 and 4 operands,
## several I/O blocks, operands in different files), explicit out_dataset,
## no intermediates in /BIND after two executions, the 2-argument path, and a
## clean error (no crash, no partial output) for non-conformable operands.
## Also checks that a failing 2-argument bind no longer leaves the first
## matrix closed and unusable.
##
## Self-contained: writes its HDF5 files to tempdir(). Each case runs in a
## separate R process (callr) so a crash is reported as a FAIL.
##
## Run with:  Rscript tests_regression/test_bind_nary.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(BigDataStatMeth)
    library(rhdf5)
})

cat("BigDataStatMeth", as.character(packageVersion("BigDataStatMeth")),
    "from", find.package("BigDataStatMeth"), "\n\n")

workdir <- file.path(tempdir(), "bdsm_bind_nary")
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

## Run 'fun(workdir)' in a fresh R process. 'fun' returns a character vector
## of failure messages (empty = PASS).
run_case <- function(name, fun) {
    wd <- file.path(workdir, name)
    unlink(wd, recursive = TRUE)
    dir.create(wd, recursive = TRUE)
    body <- function(fun, wd) {
        suppressPackageStartupMessages({
            library(BigDataStatMeth)
            library(rhdf5)
        })
        set.seed(1)
        msg <- tryCatch(fun(wd), error = function(e)
            paste("unexpected error:", conditionMessage(e)))
        try(hdf5_close_all(), silent = TRUE)
        msg
    }
    msg <- tryCatch(
        callr::r(body, args = list(fun = fun, wd = wd), libpath = .libPaths()),
        error = function(e) paste("R process failed:", conditionMessage(e)))
    ok <- length(msg) == 0
    cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", name))
    if (!ok) cat(paste0("       ", msg, collapse = "\n"), "\n")
    ok
}

## ── Case bodies (run in the child process) ──────────────────────────────────

case_cbind_3 <- function(wd) {
    f <- file.path(wd, "t.h5")
    m <- list(matrix(rnorm(7 * 2), 7), matrix(rnorm(7 * 3), 7),
              matrix(rnorm(7 * 1), 7))
    h <- lapply(seq_along(m), function(i)
        hdf5_create_matrix(f, paste0("data/M", i), data = m[[i]]))
    r <- cbind(h[[1]], h[[2]], h[[3]], block_rows = 3L)
    msg <- character(0)
    if (!identical(as.matrix(r), do.call(cbind, m)))
        msg <- c(msg, "values differ from base::cbind")
    if (!grepl("^cbind_3_[0-9a-f]{8}$", r$get_dataset()))
        msg <- c(msg, paste("unexpected name", r$get_dataset()))
    if (r$get_group() != "BIND") msg <- c(msg, paste("group", r$get_group()))
    if (!is.null(dimnames(r)))   msg <- c(msg, "dimnames not NULL")
    msg
}

case_cbind_4 <- function(wd) {
    f <- file.path(wd, "t.h5")
    m <- list(matrix(rnorm(5 * 2), 5), matrix(rnorm(5 * 4), 5),
              matrix(rnorm(5 * 1), 5), matrix(rnorm(5 * 3), 5))
    h <- lapply(seq_along(m), function(i)
        hdf5_create_matrix(f, paste0("data/M", i), data = m[[i]]))
    r <- cbind(h[[1]], h[[2]], h[[3]], h[[4]])
    msg <- character(0)
    if (!identical(as.matrix(r), do.call(cbind, m)))
        msg <- c(msg, "values differ from base::cbind")
    if (!grepl("^cbind_4_[0-9a-f]{8}$", r$get_dataset()))
        msg <- c(msg, paste("unexpected name", r$get_dataset()))
    msg
}

case_rbind_3 <- function(wd) {
    f <- file.path(wd, "t.h5")
    m <- list(matrix(rnorm(3 * 4), 3), matrix(rnorm(5 * 4), 5),
              matrix(rnorm(1 * 4), 1))
    h <- lapply(seq_along(m), function(i)
        hdf5_create_matrix(f, paste0("data/M", i), data = m[[i]]))
    r <- rbind(h[[1]], h[[2]], h[[3]], block_rows = 2L)
    msg <- character(0)
    if (!identical(as.matrix(r), do.call(rbind, m)))
        msg <- c(msg, "values differ from base::rbind")
    if (!grepl("^rbind_3_[0-9a-f]{8}$", r$get_dataset()))
        msg <- c(msg, paste("unexpected name", r$get_dataset()))
    if (r$get_group() != "BIND") msg <- c(msg, paste("group", r$get_group()))
    msg
}

case_rbind_4 <- function(wd) {
    f <- file.path(wd, "t.h5")
    m <- list(matrix(rnorm(2 * 3), 2), matrix(rnorm(4 * 3), 4),
              matrix(rnorm(3 * 3), 3), matrix(rnorm(6 * 3), 6))
    h <- lapply(seq_along(m), function(i)
        hdf5_create_matrix(f, paste0("data/M", i), data = m[[i]]))
    r <- rbind(h[[1]], h[[2]], h[[3]], h[[4]])
    msg <- character(0)
    if (!identical(as.matrix(r), do.call(rbind, m)))
        msg <- c(msg, "values differ from base::rbind")
    if (!grepl("^rbind_4_[0-9a-f]{8}$", r$get_dataset()))
        msg <- c(msg, paste("unexpected name", r$get_dataset()))
    msg
}

## Same matrix passed more than once; integer input stored as double output.
case_repeated_operand <- function(wd) {
    f <- file.path(wd, "t.h5")
    a <- matrix(1:12, 4); b <- matrix(rnorm(4 * 2), 4)
    A <- hdf5_create_matrix(f, "data/A", data = a)
    B <- hdf5_create_matrix(f, "data/B", data = b)
    rc <- cbind(A, B, A)
    rr <- rbind(A, A, A, A)
    msg <- character(0)
    if (!identical(as.matrix(rc), cbind(a + 0, b, a + 0)))
        msg <- c(msg, "cbind(A, B, A) values differ")
    if (!identical(as.matrix(rr), rbind(a, a, a, a) + 0))
        msg <- c(msg, "rbind(A, A, A, A) values differ")
    if (!A$is_valid()) msg <- c(msg, "A no longer valid")
    msg
}

## Operands in three different files: result goes to the first file.
case_cross_file <- function(wd) {
    fs <- file.path(wd, paste0("f", 1:3, ".h5"))
    mc <- list(matrix(rnorm(6 * 2), 6), matrix(rnorm(6 * 3), 6),
               matrix(rnorm(6 * 2), 6))
    mr <- list(matrix(rnorm(2 * 5), 2), matrix(rnorm(3 * 5), 3),
               matrix(rnorm(4 * 5), 4))
    hc <- lapply(1:3, function(i)
        hdf5_create_matrix(fs[i], "data/C", data = mc[[i]]))
    hr <- lapply(1:3, function(i)
        hdf5_create_matrix(fs[i], "data/R", data = mr[[i]]))
    rc <- cbind(hc[[1]], hc[[2]], hc[[3]])
    rr <- rbind(hr[[1]], hr[[2]], hr[[3]])
    msg <- character(0)
    if (!identical(as.matrix(rc), do.call(cbind, mc)))
        msg <- c(msg, "cross-file cbind values differ")
    if (!identical(as.matrix(rr), do.call(rbind, mr)))
        msg <- c(msg, "cross-file rbind values differ")
    if (normalizePath(rc$get_filename()) != normalizePath(fs[1]))
        msg <- c(msg, paste("cbind result in", rc$get_filename()))
    if (normalizePath(rr$get_filename()) != normalizePath(fs[1]))
        msg <- c(msg, paste("rbind result in", rr$get_filename()))
    msg
}

## Explicit out_dataset / out_group / out_file / overwrite are respected.
case_explicit_output <- function(wd) {
    f  <- file.path(wd, "t.h5")
    fo <- file.path(wd, "out.h5")
    m <- list(matrix(rnorm(4 * 2), 4), matrix(rnorm(4 * 2), 4),
              matrix(rnorm(4 * 2), 4))
    h <- lapply(seq_along(m), function(i)
        hdf5_create_matrix(f, paste0("data/M", i), data = m[[i]]))
    r1 <- cbind(h[[1]], h[[2]], h[[3]], out_dataset = "my_cbind")
    r2 <- rbind(h[[1]], h[[2]], h[[3]], out_dataset = "my_rbind",
                out_group = "RES", out_file = fo)
    msg <- character(0)
    if (r1$get_path() != "BIND/my_cbind" && r1$get_path() != "/BIND/my_cbind")
        msg <- c(msg, paste("cbind path", r1$get_path()))
    if (!identical(as.matrix(r1), do.call(cbind, m)))
        msg <- c(msg, "cbind values differ")
    if (!(r2$get_dataset() == "my_rbind" && r2$get_group() == "RES" &&
          normalizePath(r2$get_filename()) == normalizePath(fo)))
        msg <- c(msg, paste("rbind location", r2$get_filename(), r2$get_path()))
    if (!identical(as.matrix(r2), do.call(rbind, m)))
        msg <- c(msg, "rbind values differ")
    ## overwrite = FALSE (default) refuses an existing output and keeps it
    e <- tryCatch({
        cbind(h[[3]], h[[2]], h[[1]], out_dataset = "my_cbind")
        NULL
    }, error = function(e) conditionMessage(e))
    if (is.null(e))
        msg <- c(msg, "existing output overwritten with overwrite = FALSE")
    if (!h[[1]]$is_valid()) msg <- c(msg, "first operand invalid after error")
    kept <- as.matrix(hdf5_matrix(f, "BIND/my_cbind"))
    if (!identical(kept, do.call(cbind, m)))
        msg <- c(msg, "existing output changed by a refused bind")
    ## overwrite = TRUE replaces the explicit output
    r3 <- cbind(h[[3]], h[[2]], h[[1]], out_dataset = "my_cbind",
                overwrite = TRUE)
    if (!identical(as.matrix(r3), cbind(m[[3]], m[[2]], m[[1]])))
        msg <- c(msg, "overwrite = TRUE did not replace the output")
    msg
}

## Two executions in the same file: /BIND holds only the final results.
case_no_intermediates <- function(wd) {
    f <- file.path(wd, "t.h5")
    m <- lapply(1:4, function(i) matrix(rnorm(4 * 4), 4))
    h <- lapply(seq_along(m), function(i)
        hdf5_create_matrix(f, paste0("data/M", i), data = m[[i]]))
    results <- character(0)
    for (run in 1:2) {
        rs <- list(cbind(h[[1]], h[[2]], h[[3]]),
                   cbind(h[[1]], h[[2]], h[[3]], h[[4]]),
                   rbind(h[[1]], h[[2]], h[[3]]),
                   rbind(h[[1]], h[[2]], h[[3]], h[[4]]))
        results <- c(results, vapply(rs, function(r) r$get_dataset(), ""))
    }
    hdf5_close_all()
    ls  <- rhdf5::h5ls(f, recursive = TRUE)
    rhdf5::h5closeAll()
    bind <- ls$name[ls$group == "/BIND" & ls$otype == "H5I_DATASET"]
    msg <- character(0)
    if (length(unique(results)) != 8L)
        msg <- c(msg, "result names are not unique across executions")
    extra <- setdiff(bind, results)
    if (length(extra))
        msg <- c(msg, paste("unexpected datasets in /BIND:",
                            paste(extra, collapse = ", ")))
    if (any(grepl("_s[0-9]+$", bind)))
        msg <- c(msg, "intermediate '_sK' datasets found")
    if (length(setdiff(results, bind)))
        msg <- c(msg, "some results missing from /BIND")
    msg
}

## Two arguments: unchanged behaviour (name "A_cbind_B", values).
case_two_args <- function(wd) {
    f <- file.path(wd, "t.h5")
    a <- matrix(rnorm(5 * 3), 5); b <- matrix(rnorm(5 * 3), 5)
    A <- hdf5_create_matrix(f, "data/A", data = a)
    B <- hdf5_create_matrix(f, "data/B", data = b)
    rc <- cbind(A, B); rr <- rbind(A, B)
    rx <- cbind(A, B, out_dataset = "AB")
    msg <- character(0)
    if (rc$get_dataset() != "A_cbind_B") msg <- c(msg, rc$get_dataset())
    if (rr$get_dataset() != "A_rbind_B") msg <- c(msg, rr$get_dataset())
    if (rx$get_dataset() != "AB")        msg <- c(msg, rx$get_dataset())
    if (!identical(as.matrix(rc), cbind(a, b))) msg <- c(msg, "cbind values")
    if (!identical(as.matrix(rr), rbind(a, b))) msg <- c(msg, "rbind values")
    if (!identical(as.matrix(rx), cbind(a, b))) msg <- c(msg, "AB values")
    if (!identical(cbind(A, A)$get_dataset(), "A_cbind_A"))
        msg <- c(msg, "cbind(A, A) name")
    msg
}

## Two arguments, failing bind (existing output, overwrite = FALSE): clean
## error, and the first matrix stays valid and usable (it used to be left
## closed).
case_two_args_error_keeps_self <- function(wd) {
    f <- file.path(wd, "t.h5")
    a <- matrix(rnorm(4 * 2), 4); b <- matrix(rnorm(3 * 2), 3)
    A <- hdf5_create_matrix(f, "data/A", data = a)
    B <- hdf5_create_matrix(f, "data/B", data = b)
    msg <- character(0)
    for (fn in c("cbind", "rbind")) {
        bind <- get(fn)
        Y <- if (fn == "cbind") A else B   # cbind(A, A), rbind(A, B)
        y <- if (fn == "cbind") a else b
        ds <- paste0("fixed_", fn)
        bind(A, Y, out_dataset = ds)
        e <- tryCatch({ bind(A, Y, out_dataset = ds); NULL },
                      error = function(e) conditionMessage(e))
        if (is.null(e))
            msg <- c(msg, paste(fn, "rerun without overwrite did not fail"))
        if (any(grepl("c\\+\\+ exception", e, ignore.case = TRUE)))
            msg <- c(msg, "error message mentions 'c++ exception'")
        if (!A$is_valid()) {
            msg <- c(msg, paste("A not valid after failed", fn))
            next
        }
        if (!identical(as.matrix(A), a))
            msg <- c(msg, paste("A not readable after failed", fn))
        r <- tryCatch(bind(A, Y, out_dataset = ds, overwrite = TRUE),
                      error = function(e) NULL)
        if (is.null(r) || !identical(as.matrix(r), base::get(fn)(a, y)))
            msg <- c(msg, paste("A not usable in a new", fn, "after failure"))
    }
    msg
}

## Non-conformable 3rd operand: clean error, inputs still valid, no output.
case_bad_dims <- function(wd) {
    f <- file.path(wd, "t.h5")
    A <- hdf5_create_matrix(f, "data/A", data = matrix(rnorm(6), 3))
    B <- hdf5_create_matrix(f, "data/B", data = matrix(rnorm(6), 3))
    D <- hdf5_create_matrix(f, "data/D", data = matrix(rnorm(12), 4))
    msg <- character(0)
    ec <- tryCatch({ cbind(A, B, D); NULL },
                   error = function(e) conditionMessage(e))
    er <- tryCatch({ rbind(A, B, D); NULL },
                   error = function(e) conditionMessage(e))
    if (is.null(ec) || !grepl("nrow", ec))
        msg <- c(msg, paste("cbind error:", if (is.null(ec)) "none" else ec))
    if (is.null(er) || !grepl("ncol", er))
        msg <- c(msg, paste("rbind error:", if (is.null(er)) "none" else er))
    if (any(grepl("c\\+\\+ exception", c(ec, er), ignore.case = TRUE)))
        msg <- c(msg, "error message mentions 'c++ exception'")
    if (!all(vapply(list(A, B, D), function(x) x$is_valid(), TRUE)))
        msg <- c(msg, "an input is no longer valid after the error")
    ## The handles still work after the error
    if (!identical(dim(cbind(A, B, A)), c(3L, 6L)))
        msg <- c(msg, "bind after error failed")
    hdf5_close_all()
    ls <- rhdf5::h5ls(f, recursive = TRUE)
    rhdf5::h5closeAll()
    bind <- ls$name[ls$group == "/BIND" & ls$otype == "H5I_DATASET"]
    if (length(bind) != 1L)
        msg <- c(msg, paste("datasets left in /BIND:",
                            paste(bind, collapse = ", ")))
    msg
}

cases <- list(
    cbind_3_operands          = case_cbind_3,
    cbind_4_operands          = case_cbind_4,
    rbind_3_operands          = case_rbind_3,
    rbind_4_operands          = case_rbind_4,
    repeated_operand          = case_repeated_operand,
    cross_file                = case_cross_file,
    explicit_output           = case_explicit_output,
    no_intermediates_2_runs   = case_no_intermediates,
    two_arguments_unchanged   = case_two_args,
    two_args_error_keeps_self = case_two_args_error_keeps_self,
    nonconformable_3rd_clean  = case_bad_dims
)

results <- vapply(names(cases), function(nm) run_case(nm, cases[[nm]]),
                  logical(1))

unlink(workdir, recursive = TRUE)

cat(sprintf("\nN-ary bind regression: %d/%d PASS\n",
            sum(results), length(results)))
if (!all(results)) {
    stop("N-ary bind regression test FAILED: ",
         paste(names(cases)[!results], collapse = ", "))
}
cat("ALL PASS\n")
