# S3_standalone.R
#
# Standalone functions for group-level HDF5 operations.
# These operate on an entire HDF5 group, not on a single HDF5Matrix,
# so they are plain functions (no R6 method, no S3 dispatch).
#
#   hdf5_reduce()   — reduce all datasets in a group via "+" or "-"
#   hdf5_apply()    — apply one operation to multiple datasets in a group
#   list_datasets() — lists all datasets within a group of an HDF5Matrix file.

# ── hdf5_reduce() ─────────────────────────────────────────────────────────────

#' Reduce all datasets in an HDF5 group by a binary operation
#'
#' @description
#' Standalone function that applies a binary reduction (\code{"+"} or
#' \code{"-"}) across all datasets stored in a given HDF5 group and writes
#' the accumulated result as a new dataset. No open \code{HDF5Matrix} object
#' is required.
#'
#' @details
#' Use \code{hdf5_reduce()} when you only have the file path and group name.
#' If you already have an open \code{HDF5Matrix}, use
#' \code{\link{reduce}(x, ...)} instead — both produce the same result.
#'
#' All datasets in \code{group} must have the same dimensions.
#'
#' @param filename   Path to the HDF5 file.
#' @param group      Group path containing the datasets to reduce.
#' @param func       Character. Reduction operator: \code{"+"} (default)
#'   or \code{"-"}.
#' @param outgroup   Character or \code{NULL}. Output group
#'   (default: same as \code{group}).
#' @param outdataset Character or \code{NULL}. Output dataset name
#'   (default: same as \code{group}).
#' @param overwrite  Logical. Overwrite existing output dataset.
#' @param remove     Logical. Remove input datasets after reduction.
#' @return An \code{HDF5Matrix} pointing to the result dataset.
#'
#' @examples
#' \donttest{
#' fn <- tempfile(fileext = ".h5")
#' hdf5_create_matrix(fn, "blocks/A", data = matrix(1:6, 2, 3))
#' hdf5_create_matrix(fn, "blocks/B", data = matrix(1:6, 2, 3))
#' hdf5_create_matrix(fn, "blocks/C", data = matrix(1:6, 2, 3))
#' result <- hdf5_reduce(fn, group = "blocks", func = "+")
#' as.matrix(result)
#' hdf5_close_all()
#' unlink(fn)
#' }
#'
#' @seealso \code{\link{hdf5_apply}}, \code{\link{cbind.HDF5Matrix}}
#' @export
hdf5_reduce <- function(filename,
                        group,
                        func        = "+",
                        outgroup    = NULL,
                        outdataset  = NULL,
                        overwrite   = FALSE,
                        remove      = FALSE) {

    func <- match.arg(func, c("+", "-"))

    res <- bdReduce_hdf5_dataset(
        filename       = filename,
        group          = group,
        reducefunction = func,
        outgroup       = outgroup,
        outdataset     = outdataset,
        overwrite      = overwrite,
        remove         = remove
    )

    hdf5_matrix(res$fn, res$ds)
}


# ── hdf5_apply() ──────────────────────────────────────────────────────────────

#' Apply a mathematical operation to multiple HDF5 datasets
#'
#' @description
#' Standalone function that applies a predefined algebraic or statistical
#' operation to a list of datasets stored in an HDF5 group, writing results
#' to \code{outgroup}. No open \code{HDF5Matrix} object is required.
#'
#' \strong{This is not equivalent to \code{base::apply()}.} It dispatches
#' built-in C++ operations (QR, cross-product, Cholesky, etc.) to a batch of
#' named datasets in the file.
#' 
#' 
#'
#' @details
#' Use \code{hdf5_apply()} when you only have the file path and group name.
#' If you already have an open \code{HDF5Matrix}, use
#' \code{\link{apply_function}(x, ...)} instead — both produce the same result.
#'
#' \describe{
#'   \item{\code{"CrossProd"}}{Compute \eqn{A^T A} for each dataset.}
#'   \item{\code{"tCrossProd"}}{Compute \eqn{A A^T} for each dataset.}
#'   \item{\code{"CrossProd_double"}}{Double-precision cross-product.}
#'   \item{\code{"tCrossProd_double"}}{Double-precision transposed cross-product.}
#'   \item{\code{"blockmult"}}{Block-wise \eqn{A \times B}
#'     (requires \code{b_datasets}).}
#'   \item{\code{"QR"}}{QR decomposition for each dataset.}
#'   \item{\code{"invChol"}}{Inverse via Cholesky for each dataset.}
#'   \item{\code{"solve"}}{Solve \eqn{AX = B}
#'     (requires \code{b_datasets}).}
#'   \item{\code{"normalize"}}{Column-wise normalization.}
#'   \item{\code{"sdmean"}}{Compute SD and mean.}
#'   \item{\code{"descChol"}}{Cholesky decomposition.}
#' }
#'
#' @param filename         Path to the HDF5 file.
#' @param group            Group path containing \code{datasets}.
#' @param datasets         Character vector of dataset names to process.
#' @param func             Character. Operation to apply (see Details).
#' @param outgroup         Character. Output group path for results.
#' @param b_group          Character or \code{NULL}. Group of B datasets.
#' @param b_datasets       Character vector or \code{NULL}. Names of B datasets.
#' @param overwrite        Logical. Overwrite existing output datasets.
#' @param transp_dataset   Logical. Transpose A datasets before operation.
#' @param transp_bdataset  Logical. Transpose B datasets before operation.
#' @param fullMatrix       Logical. Return full matrix (not triangular).
#' @param byrows           Logical. Apply by rows (for normalize/sdmean).
#' @param threads          Integer or \code{NULL}. OpenMP threads.
#' @return Invisibly \code{NULL}. Results written to \code{outgroup}.
#'   Open them with \code{\link{hdf5_matrix}()}.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' A <- hdf5_create_matrix(tmp, "inp/A", data = matrix(rnorm(25), 5, 5))
#' B <- hdf5_create_matrix(tmp, "inp/B", data = matrix(rnorm(25), 5, 5))
#' hdf5_apply(tmp, group = "inp", datasets = c("A", "B"), 
#' func = "CrossProd", outgroup = "out")
#' res <- list_datasets(tmp)
#' res
#' res_A <- hdf5_matrix(tmp, res[3])
#' dim(res_A)   # 5 x 5
#' close(res_A)
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{apply_function}} for the S3 method on an open
#'   \code{HDF5Matrix}; \code{\link{hdf5_reduce}} for group-level reduction.
#'   
#' @export
hdf5_apply <- function(filename,
                       group,
                       datasets,
                       func,
                       outgroup,
                       b_group         = NULL,
                       b_datasets      = NULL,
                       overwrite       = FALSE,
                       transp_dataset  = FALSE,
                       transp_bdataset = FALSE,
                       fullMatrix      = FALSE,
                       byrows          = FALSE,
                       threads         = NULL) {

    valid_funcs <- c("QR", "CrossProd", "tCrossProd",
                     "invChol", "blockmult",
                     "CrossProd_double", "tCrossProd_double",
                     "solve", "normalize", "sdmean", "descChol")
    func <- match.arg(func, valid_funcs)

    bdapply_Function_hdf5(
        filename        = filename,
        group           = group,
        datasets        = datasets,
        outgroup        = outgroup,
        func            = func,
        b_group         = b_group,
        b_datasets      = if (is.null(b_datasets)) NULL else b_datasets,
        overwrite       = overwrite,
        transp_dataset  = transp_dataset,
        transp_bdataset = transp_bdataset,
        fullMatrix      = fullMatrix,
        byrows          = byrows,
        threads         = if (is.null(threads)) 2L else as.integer(threads)
    )

    invisible(NULL)
}


# ──────────────────────────────────────────────────────────────────────────────
#   hdf5_remove() — delete a dataset from an HDF5 file by path.
# ──────────────────────────────────────────────────────────────────────────────

#' Remove (delete) a dataset from an HDF5 file
#'
#' @description
#' Deletes a dataset from an HDF5 file, given the file path and the full
#' internal dataset path (\code{"group/name"}). No open \code{HDF5Matrix}
#' object is required. If you already have an open \code{HDF5Matrix}, you can
#' equivalently call its \code{$remove()} method.
#'
#' @details
#' \strong{Space is not reclaimed.} HDF5 only \emph{unlinks} the dataset: the
#' link is removed so the dataset can no longer be opened, but the disk space
#' it occupied is \strong{not} returned to the operating system until the file
#' is rewritten. To physically reclaim the space, repack the file with the
#' HDF5 command-line tool, e.g. \code{h5repack old.h5 new.h5}, or copy the
#' datasets you want to keep into a fresh file. Consequently, repeatedly
#' creating and removing datasets in the same file can make it grow on disk
#' even though the logical contents shrink.
#'
#' Removing a dataset that does not exist raises an error.
#'
#' @param filename Path to the HDF5 file.
#' @param dataset  Full internal path of the dataset to remove, in the form
#'   \code{"group/name"} (e.g. \code{"data/matrix"} or
#'   \code{"grp/sub/matrix"}).
#'
#' @return \code{TRUE} invisibly on success.
#'
#' @examples
#' \donttest{
#' fn <- tempfile(fileext = ".h5")
#' hdf5_create_matrix(fn, "grp/keep", data = matrix(1:6, 2, 3))
#' hdf5_create_matrix(fn, "grp/tmp",  data = matrix(1:6, 2, 3))
#' hdf5_remove(fn, "grp/tmp")
#' list_datasets(fn, group = "grp")   # only "keep" remains
#' hdf5_close_all()
#' unlink(fn)
#' }
#'
#' @seealso \code{\link{hdf5_matrix}}, \code{\link{hdf5_create_matrix}},
#'   \code{\link{list_datasets}}
#' @export
hdf5_remove <- function(filename, dataset) {
    if (!is.character(filename) || length(filename) != 1) {
        stop("filename must be a single string")
    }
    if (!is.character(dataset) || length(dataset) != 1) {
        stop("dataset must be a single string")
    }
    if (!file.exists(filename)) {
        stop("File does not exist: ", filename)
    }

    parts <- strsplit(dataset, "/")[[1]]
    parts <- parts[nzchar(parts)]
    if (length(parts) < 2) {
        stop("dataset must be in the format 'group/name' or ",
             "'group/subgroup/name'")
    }
    name  <- parts[length(parts)]
    group <- paste(parts[-length(parts)], collapse = "/")

    # Best-effort: invalidate any live HDF5Matrix handles pointing at this
    # dataset so their is_valid() flips to FALSE after the unlink.
    tryCatch(
        rcpp_hdf5_close_at_paths(filename, paste(group, name, sep = "/")),
        error = function(e) NULL
    )

    invisible(rcpp_hdf5_remove_dataset(filename, group, name))
}
