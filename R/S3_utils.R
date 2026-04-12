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