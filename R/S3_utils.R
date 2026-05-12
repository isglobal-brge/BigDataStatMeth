#' List datasets in an HDF5 file or group
#'
#' @description
#' Lists datasets within an HDF5 file. If no group is specified, the entire
#' file is traversed recursively and full relative paths are returned
#' (e.g. \code{"INPUT/A"}, \code{"RESULTS/SVD/d"}). If a group is given,
#' only the datasets in that group are listed unless \code{recursive = TRUE}.
#'
#' @param x         An \code{HDF5Matrix} object, or a character string
#'   (file path).
#' @param group     Character or \code{NULL}. Group path inside the HDF5
#'   file. If \code{NULL} and \code{x} is a file path, the whole file is
#'   listed recursively. If \code{NULL} and \code{x} is an
#'   \code{HDF5Matrix}, the object's own group is used (non-recursive).
#' @param prefix    Optional character. Only return datasets whose name
#'   starts with this prefix.
#' @param recursive Logical. If \code{TRUE}, recurse into subgroups and
#'   return full relative paths. Default \code{FALSE}.
#'
#' @return Character vector of dataset names or relative paths.
#'
#' @examples
#' \donttest{
#' fn <- tempfile(fileext = ".h5")
#' hdf5_create_matrix(fn, "INPUT/A",    data = matrix(rnorm(100), 10, 10))
#' hdf5_create_matrix(fn, "INPUT/B",    data = matrix(rnorm(100), 10, 10))
#' hdf5_create_matrix(fn, "RESULTS/C",  data = matrix(rnorm(100), 10, 10))
#'
#' # All datasets in the file (recursive from root)
#' list_datasets(fn)
#'
#' # Only INPUT group
#' list_datasets(fn, group = "INPUT")
#'
#' # From an HDF5Matrix object (uses object's own group)
#' X <- hdf5_matrix(fn, "INPUT/A")
#' list_datasets(X)
#'
#' hdf5_close_all()
#' unlink(fn)
#' }
#'
#' @seealso \code{\link{hdf5_matrix}}, \code{\link{hdf5_create_matrix}}
#' @export
list_datasets <- function(x, group = NULL, prefix = NULL, recursive = FALSE) {
    if (inherits(x, "HDF5Matrix")) {
        filename  <- x$get_filename()
        if (is.null(group)) {
            group     <- x$get_group()
            recursive <- FALSE
        }
    } else if (is.character(x)) {
        filename <- x
        if (is.null(group)) {
            group     <- "/"
            recursive <- TRUE
        }
    } else {
        stop("x must be an HDF5Matrix or a file path")
    }
    bdgetDatasetsList_hdf5(filename, group = group,
                           prefix = prefix, recursive = recursive)
}