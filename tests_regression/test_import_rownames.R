## test_import_rownames.R
## -----------------------------------------------------------------------------
## Regression test: bdImportTextFile_hdf5() with rownames = TRUE.
##
## Before the fix, bdImportTextFile_hdf5(..., rownames = TRUE) always failed
## with "attempt to set index N/N in SET_STRING_ELT". The rownames-column check
## in Rcpp_Import_File_to_hdf5() reads one extra data row before the newline
## count, so the row count (and the rownames buffer) was one short and writing
## the last row name overflowed it. CRLF cases check that no trailing "\r" is
## left in the dimnames or values; 'sci_*' and 'special_chars_*' cases check
## that data rows are split on the separator only (values such as 1.5e+02 and
## row names containing '+', '[', ']' or '^' used to be split as well).
##
## Self-contained: writes its CSV files to tempdir(), imports them, reads the
## result back with rhdf5 and checks dims, dimnames and values.
##
## Run with:  Rscript tests_regression/test_import_rownames.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(BigDataStatMeth)
    library(rhdf5)
})

cat("BigDataStatMeth", as.character(packageVersion("BigDataStatMeth")),
    "from", find.package("BigDataStatMeth"), "\n\n")

workdir <- file.path(tempdir(), "bdsm_import_rownames")
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

expected_values <- matrix(c(0.1, 0.2,
                            0.3, 0.4,
                            0.5, 0.6), nrow = 3, byrow = TRUE)
expected_rows <- c("cg1", "cg2", "cg3")
expected_cols <- c("s1", "s2")

## Write 'lines' to a file, with or without a trailing newline. 'eol' sets
## the exact line ending ("\n" for LF, "\r\n" for CRLF).
write_csv_lines <- function(path, lines, trailing_newline = TRUE, eol = "\n") {
    txt <- paste(lines, collapse = eol)
    if (trailing_newline) txt <- paste0(txt, eol)
    writeChar(txt, path, eos = NULL, useBytes = TRUE)
    invisible(path)
}

## Import 'csv' into 'h5f' and read it back. Runs in a separate R process
## (callr) so that a crash in the importer is reported as a FAIL instead of
## aborting the whole script.
import_and_read <- function(csv, h5f, header, rownames) {
    suppressPackageStartupMessages({
        library(BigDataStatMeth)
        library(rhdf5)
    })
    ## One dimnames component ("1" = rows, "2" = cols); NULL if absent.
    read_dimnames <- function(which) {
        path <- paste0("/IMPORT/.data_dimnames/", which)
        contents <- rhdf5::h5ls(h5f, recursive = TRUE)
        full <- gsub("/+", "/", paste0(contents$group, "/", contents$name))
        if (!(path %in% full)) return(NULL)
        x <- rhdf5::h5read(h5f, path)
        if (is.list(x) || is.data.frame(x)) x <- x[[1]]
        as.character(x)
    }
    err <- tryCatch({
        bdImportTextFile_hdf5(filename = csv, outputfile = h5f,
                              outGroup = "IMPORT", outDataset = "data",
                              sep = ",", header = header, rownames = rownames,
                              overwrite = TRUE)
        NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(err)) return(list(err = err))
    rhdf5::h5closeAll()
    out <- list(err = NULL, vals = rhdf5::h5read(h5f, "IMPORT/data"),
                rn = read_dimnames("1"), cn = read_dimnames("2"))
    rhdf5::h5closeAll()
    out
}

run_case <- function(name, lines, header, rownames, trailing_newline,
                     check_rows, check_cols, eol = "\n",
                     exp_vals = expected_values, exp_rows = expected_rows,
                     exp_cols = expected_cols, exact = FALSE) {
    csv <- file.path(workdir, paste0(name, ".csv"))
    h5f <- file.path(workdir, paste0(name, ".h5"))
    if (file.exists(h5f)) unlink(h5f)
    write_csv_lines(csv, lines, trailing_newline, eol)

    msg <- character(0)
    out <- tryCatch(
        callr::r(import_and_read,
                 args = list(csv = csv, h5f = h5f,
                             header = header, rownames = rownames),
                 libpath = .libPaths()),
        error = function(e) list(err = paste("R process failed:",
                                             conditionMessage(e))))

    if (!is.null(out$err)) {
        msg <- paste("import error:", out$err)
    } else {
        vals <- out$vals
        if (!identical(dim(vals), dim(exp_vals))) {
            msg <- c(msg, sprintf("dims %s, expected %s",
                                  paste(dim(vals), collapse = "x"),
                                  paste(dim(exp_vals), collapse = "x")))
        } else if (if (exact) !identical(c(vals), c(exp_vals))
                   else max(abs(vals - exp_vals)) > 1e-12) {
            msg <- c(msg, sprintf("values differ: %s",
                                  paste(format(c(vals)), collapse = " ")))
        }
        if (check_rows && !identical(out$rn, exp_rows))
            msg <- c(msg, sprintf("rownames = [%s]",
                                  paste(encodeString(out$rn), collapse = ",")))
        if (check_cols && !identical(out$cn, exp_cols))
            msg <- c(msg, sprintf("colnames = [%s]",
                                  paste(encodeString(out$cn), collapse = ",")))
    }

    ok <- length(msg) == 0
    cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", name))
    if (!ok) cat(paste0("       ", msg, collapse = "\n"), "\n")
    ok
}

## Values written in scientific notation in the 'sci_*' cases
sci_values <- matrix(c(150, 0.015,
                       2000, -0.45,
                       1, 6.25e10), nrow = 3, byrow = TRUE)

data_rn   <- c("cg1,0.1,0.2", "cg2,0.3,0.4", "cg3,0.5,0.6")
data_norn <- c("0.1,0.2", "0.3,0.4", "0.5,0.6")

cases <- list(
    ## Minimal reproducer
    list(name = "header_T_rownames_T", lines = c("id,s1,s2", data_rn),
         header = TRUE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = TRUE),
    ## Header without a field for the rownames column (write.table style)
    list(name = "header_T_rownames_T_short_header", lines = c("s1,s2", data_rn),
         header = TRUE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = TRUE),
    list(name = "header_F_rownames_T", lines = data_rn,
         header = FALSE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = FALSE),
    ## No trailing newline
    list(name = "header_T_rownames_T_no_final_newline",
         lines = c("id,s1,s2", data_rn),
         header = TRUE, rownames = TRUE, nl = FALSE, rows = TRUE, cols = TRUE),
    list(name = "header_F_rownames_T_no_final_newline", lines = data_rn,
         header = FALSE, rownames = TRUE, nl = FALSE, rows = TRUE, cols = FALSE),
    ## Non-regression: no rownames column
    list(name = "header_T_rownames_F", lines = c("s1,s2", data_norn),
         header = TRUE, rownames = FALSE, nl = TRUE, rows = FALSE, cols = TRUE),
    list(name = "header_T_rownames_F_no_final_newline",
         lines = c("s1,s2", data_norn),
         header = TRUE, rownames = FALSE, nl = FALSE, rows = FALSE, cols = TRUE),
    list(name = "header_F_rownames_F", lines = data_norn,
         header = FALSE, rownames = FALSE, nl = TRUE, rows = FALSE, cols = FALSE),
    ## CRLF line endings: no trailing '\r' in names or values
    list(name = "crlf_header_T_rownames_T", lines = c("id,s1,s2", data_rn),
         header = TRUE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = TRUE,
         eol = "\r\n"),
    list(name = "crlf_header_F_rownames_T", lines = data_rn,
         header = FALSE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = FALSE,
         eol = "\r\n"),
    list(name = "crlf_header_T_rownames_F", lines = c("s1,s2", data_norn),
         header = TRUE, rownames = FALSE, nl = TRUE, rows = FALSE, cols = TRUE,
         eol = "\r\n"),
    list(name = "crlf_header_T_rownames_T_no_final_newline",
         lines = c("id,s1,s2", data_rn),
         header = TRUE, rownames = TRUE, nl = FALSE, rows = TRUE, cols = TRUE,
         eol = "\r\n"),
    ## Scientific notation: values must not be split at '+' / '-'
    list(name = "sci_header_T_rownames_T",
         lines = c("id,s1,s2", "cg1,1.5e+02,1.5e-02", "cg2,2E+3,-4.5e-1",
                   "cg3,1e+00,6.25e+10"),
         header = TRUE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = TRUE,
         vals = sci_values, exact = TRUE),
    list(name = "sci_header_T_rownames_F",
         lines = c("s1,s2", "1.5e+02,1.5e-02", "2E+3,-4.5e-1", "1e+00,6.25e+10"),
         header = TRUE, rownames = FALSE, nl = TRUE, rows = FALSE, cols = TRUE,
         vals = sci_values, exact = TRUE),
    list(name = "crlf_sci_header_T_rownames_T",
         lines = c("id,s1,s2", "cg1,1.5e+02,1.5e-02", "cg2,2E+3,-4.5e-1",
                   "cg3,1e+00,6.25e+10"),
         header = TRUE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = TRUE,
         eol = "\r\n", vals = sci_values, exact = TRUE),
    ## Names containing '+', '[', ']' or '^' (characters of the old split set)
    list(name = "special_chars_in_names",
         lines = c("id,s+1,s[2]", "cg+1,0.1,0.2", "cg[2],0.3,0.4", "cg^3,0.5,0.6"),
         header = TRUE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = TRUE,
         rnames = c("cg+1", "cg[2]", "cg^3"), cnames = c("s+1", "s[2]")),
    list(name = "special_chars_in_names_header_F",
         lines = c("cg+1,0.1,0.2", "cg[2],0.3,0.4", "cg^3,0.5,0.6"),
         header = FALSE, rownames = TRUE, nl = TRUE, rows = TRUE, cols = FALSE,
         rnames = c("cg+1", "cg[2]", "cg^3"))
)

results <- vapply(cases, function(cs)
    run_case(cs$name, cs$lines, cs$header, cs$rownames, cs$nl,
             cs$rows, cs$cols,
             eol = if (is.null(cs$eol)) "\n" else cs$eol,
             exp_vals = if (is.null(cs$vals)) expected_values else cs$vals,
             exp_rows = if (is.null(cs$rnames)) expected_rows else cs$rnames,
             exp_cols = if (is.null(cs$cnames)) expected_cols else cs$cnames,
             exact = isTRUE(cs$exact)),
    logical(1))

unlink(workdir, recursive = TRUE)

cat(sprintf("\nImport rownames regression: %d/%d PASS\n",
            sum(results), length(results)))
if (!all(results)) {
    stop("Import rownames regression test FAILED: ",
         paste(vapply(cases[!results], `[[`, "", "name"), collapse = ", "))
}
cat("ALL PASS\n")
