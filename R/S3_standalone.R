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
#' Applies a binary reduction (\code{"+"} or \code{"-"}) across all
#' datasets stored in a given HDF5 group and writes the result as a new
#' dataset. Delegates to \code{bdReduce_hdf5_dataset()}.
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
#' @return An \code{\link{HDF5Matrix}} pointing to the result dataset.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' bdCreate_hdf5_matrix(tmp, matrix(1:6, 2, 3), "blocks", "A",
#'                      overwriteFile = FALSE)
#' bdCreate_hdf5_matrix(tmp, matrix(1:6, 2, 3), "blocks", "B",
#'                      overwriteFile = FALSE)
#' bdCreate_hdf5_matrix(tmp, matrix(1:6, 2, 3), "blocks", "C",
#'                      overwriteFile = FALSE)
#' result <- hdf5_reduce(tmp, group = "blocks", func = "+")
#' as.matrix(result)
#' close(result)
#' unlink(tmp)
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
#' Applies one of several supported operations to a list of datasets
#' stored in an HDF5 group. Delegates to \code{bdapply_Function_hdf5()}.
#'
#' @details
#' Supported values for \code{func}:
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
#' bdCreate_hdf5_matrix(tmp, matrix(rnorm(20), 4, 5), "inp", "A",
#'                      overwriteFile = FALSE)
#' bdCreate_hdf5_matrix(tmp, matrix(rnorm(20), 4, 5), "inp", "B",
#'                      overwriteFile = FALSE)
#' hdf5_apply(tmp, group = "inp", datasets = c("A", "B"),
#'            func = "CrossProd", outgroup = "out")
#' res_A <- hdf5_matrix(tmp, "out/A")
#' dim(res_A)   # 5 x 5
#' close(res_A)
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{hdf5_reduce}},
#'   \code{\link{crossprod.HDF5Matrix}}, \code{\link{qr.HDF5Matrix}}
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

# ── list_datasets() ───────────────────────────────────────────────────────────

#' List datasets in an HDF5 group
#'
#' @description
#' Lists all datasets within a group of an \code{HDF5Matrix} file.
#'
#' @param x   An \code{HDF5Matrix} object, or a character string (file path).
#' @param group Character. Group path inside the HDF5 file. If \code{x} is
#'   an \code{HDF5Matrix}, defaults to the object's own group.
#' @param prefix Optional character. Only return datasets starting with this prefix.
#'
#' @return Character vector of dataset names.
#' @export
list_datasets <- function(x, group = NULL, prefix = NULL) {
    if (inherits(x, "HDF5Matrix")) {
        filename <- x$get_filename()
        if (is.null(group)) group <- x$get_group()
    } else if (is.character(x)) {
        filename <- x
        if (is.null(group)) stop("group must be specified when x is a filename")
    } else {
        stop("x must be an HDF5Matrix or a file path")
    }
    bdgetDatasetsList_hdf5(filename, group, prefix = prefix)
}