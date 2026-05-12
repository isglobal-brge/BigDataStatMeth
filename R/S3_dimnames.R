# S3 methods - Dimension names (rownames / colnames)
# dimnames, rownames, colnames and their assignment forms for HDF5Matrix


#' Get dimension names of an HDF5Matrix
#'
#' @description
#' Returns the row and column names stored alongside the HDF5 dataset,
#' following the BigDataStatMeth convention. Returns \code{NULL} when no
#' names have been stored for a given dimension.
#'
#' @param x An \code{HDF5Matrix} object
#'
#' @return A list of length 2 with elements \code{[[1]]} (rownames) and
#'   \code{[[2]]} (colnames), or \code{NULL} for dimensions without names.
#'   Returns \code{NULL} when neither dimension has names.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' 
#' m <- matrix(1:6, 2, 3,
#'             dimnames = list(c("r1","r2"), c("c1","c2","c3")))
#' 
#' X <- hdf5_create_matrix(tmp, "data/mat", data = m)
#' 
#' dimnames(X)
#' rownames(X)
#' colnames(X)
#' rownames(X) <- c("row1", "row2")
#' rownames(X)
#' hdf5_close_all()
#' unlink(tmp)
#' 
#' }
#'
#' @export
dimnames.HDF5Matrix <- function(x) {
    if (!x$is_valid()) stop("HDF5Matrix object is closed or invalid")
    dn <- rcpp_hdf5dataset_read_dimnames(x$.__enclos_env__$private$ptr)
    rn <- if (length(dn$rownames) > 0) dn$rownames else NULL
    cn <- if (length(dn$colnames) > 0) dn$colnames else NULL
    if (is.null(rn) && is.null(cn)) return(NULL)
    list(rn, cn)
}


# #' @rdname dimnames.HDF5Matrix
# #' @export
# rownames.HDF5Matrix <- function(x, do.NULL = TRUE, prefix = "row") {
#     dn <- dimnames(x)
#     if (is.null(dn)) NULL else dn[[1]]
# }

#' @rdname dimnames.HDF5Matrix
#' @param do.NULL Logical. If \code{FALSE} and names are \code{NULL},
#'   generate names using \code{prefix}. Passed for compatibility with
#'   \code{base::rownames}; currently ignored.
#' @param prefix Character. Prefix for generated names when
#'   \code{do.NULL = FALSE}. Currently ignored.
#' @export
rownames.HDF5Matrix <- function(x, do.NULL = TRUE, prefix = "row") {
    dn <- dimnames(x)
    if (is.null(dn)) NULL else dn[[1]]
}


# #' @rdname dimnames.HDF5Matrix
# #' @export
# colnames.HDF5Matrix <- function(x, do.NULL = TRUE, prefix = "col") {
#     dn <- dimnames(x)
#     if (is.null(dn)) NULL else dn[[2]]
# }
#' @rdname dimnames.HDF5Matrix
#' @param do.NULL Logical. Ignored; present for base compatibility.
#' @param prefix Character. Ignored; present for base compatibility.
#' @export
colnames.HDF5Matrix <- function(x, do.NULL = TRUE, prefix = "col") {
    dn <- dimnames(x)
    if (is.null(dn)) NULL else dn[[2]]
}



#' Set dimension names on an HDF5Matrix
#'
#' @description
#' Writes row and/or column names to the HDF5 file alongside the dataset.
#' Setting a dimension name to \code{NULL} removes names for that dimension.
#'
#' @param x An \code{HDF5Matrix} object
#' @param value A list of length 2: \code{list(rownames, colnames)}.
#'   Either element may be \code{NULL} to skip that dimension.
#'
#' @return \code{x} invisibly
#'
#' @export
`dimnames<-.HDF5Matrix` <- function(x, value) {
    if (!x$is_valid()) stop("HDF5Matrix object is closed or invalid")
    if (!is.list(value) || length(value) != 2)
        stop("value must be a list of length 2")
    # info <- x$info()
    rn <- if (is.null(value[[1]])) character(0) else as.character(value[[1]])
    cn <- if (is.null(value[[2]])) character(0) else as.character(value[[2]])
    # bdWrite_hdf5_dimnames(info$filename, info$group, info$dataset, rn, cn)
    rcpp_hdf5dataset_write_dimnames(x$.__enclos_env__$private$ptr, rn, cn)
    invisible(x)
}


#' @rdname dimnames.HDF5Matrix
#' @param value Character vector of row names, or \code{NULL} to remove.
#' @usage \method{rownames}{HDF5Matrix}(x) <- value
#' @export
`rownames<-.HDF5Matrix` <- function(x, value) {
    if (!x$is_valid()) stop("HDF5Matrix object is closed or invalid")
    # info <- x$info()
    rn <- if (is.null(value)) character(0) else as.character(value)
    cn <- if (!is.null(colnames(x))) colnames(x) else character(0)
    # bdWrite_hdf5_dimnames(info$filename, info$group, info$dataset, rn, cn)
    rcpp_hdf5dataset_write_dimnames(x$.__enclos_env__$private$ptr, rn, cn)
    invisible(x)
}


#' @rdname dimnames.HDF5Matrix
#' @param value Character vector of column names, or \code{NULL} to remove.
#' @usage \method{colnames}{HDF5Matrix}(x) <- value 
#' @export
`colnames<-.HDF5Matrix` <- function(x, value) {
    if (!x$is_valid()) stop("HDF5Matrix object is closed or invalid")
    # info <- x$info()
    rn <- if (!is.null(rownames(x))) rownames(x) else character(0)
    cn <- if (is.null(value)) character(0) else as.character(value)
    # bdWrite_hdf5_dimnames(info$filename, info$group, info$dataset, rn, cn)
    rcpp_hdf5dataset_write_dimnames(x$.__enclos_env__$private$ptr, rn, cn)
    invisible(x)
}
