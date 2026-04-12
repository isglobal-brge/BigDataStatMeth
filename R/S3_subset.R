# S3 methods - Subsetting
# [.HDF5Matrix


#' Subset an HDF5Matrix
#'
#' @param x An \code{HDF5Matrix} object
#' @param i Row indices: numeric, integer, logical, or missing
#' @param j Column indices: numeric, integer, logical, or missing
#' @param drop Logical, whether to drop dimensions for single row/column
#'   results (default \code{TRUE})
#' @param ... Ignored
#' @return Numeric matrix, or vector when \code{drop = TRUE} and one dimension
#'   has length 1
#'
#' @details
#' All standard R indexing modes are supported:
#' \itemize{
#'   \item Contiguous ranges: \code{X[1:100, 1:50]}
#'   \item Non-contiguous: \code{X[c(1,3,5), c(2,4)]}
#'   \item Negative: \code{X[-c(1,2), ]} (all except rows 1 and 2)
#'   \item Logical: \code{X[row_mask, col_mask]}
#'   \item Missing: \code{X[, ]} (entire dataset)
#' }
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' rhdf5::h5createFile(tmp)
#' rhdf5::h5createGroup(tmp, "data")
#' rhdf5::h5write(matrix(1:100, 10, 10), tmp, "data/matrix")
#' X <- hdf5_matrix(tmp, "data/matrix")
#'
#' X[1:5, 1:3]               # submatrix
#' X[1, ]                    # single row as vector
#' X[1, , drop = FALSE]      # single row as matrix
#' X[, 2]                    # single column as vector
#' X[-c(1, 10), ]            # all except first and last row
#' X[c(TRUE, FALSE), ]       # logical row index
#' X[, ]                     # entire dataset
#'
#' X$close()
#' unlink(tmp)
#' }
#'
#' @export
`[.HDF5Matrix` <- function(x, i, j, drop = TRUE, ...) {
    
    dims  <- dim(x)
    nrows <- dims[1]
    ncols <- dims[2]
    
    # Handle missing indices
    if (missing(i)) i <- seq_len(nrows)
    if (missing(j)) j <- seq_len(ncols)
    
    # Logical -> numeric
    if (is.logical(i)) {
        if (length(i) != nrows) stop("Logical row index must have length ", nrows)
        i <- which(i)
    }
    if (is.logical(j)) {
        if (length(j) != ncols) stop("Logical column index must have length ", ncols)
        j <- which(j)
    }
    
    # Negative indices
    if (is.numeric(i) && any(i < 0)) {
        if (any(i > 0)) stop("Cannot mix positive and negative indices")
        i <- seq_len(nrows)[i]
    }
    if (is.numeric(j) && any(j < 0)) {
        if (any(j > 0)) stop("Cannot mix positive and negative indices")
        j <- seq_len(ncols)[j]
    }
    
    i <- as.integer(i)
    j <- as.integer(j)
    
    # Bounds check
    if (any(i < 1L) || any(i > nrows)) stop("Row indices out of bounds")
    if (any(j < 1L) || any(j > ncols)) stop("Column indices out of bounds")
    
    # Read dimnames BEFORE subset to avoid any HDF5 handle state side-effects
    dn <- dimnames(x)
    
    result <- x$subset(i, j)
    
    # Attach selected dimnames via structure() - no S3 dispatch, always works
    if (!is.null(dn)) {
        result <- structure(result, dimnames = list(
            if (!is.null(dn[[1]])) dn[[1]][i] else NULL,
            if (!is.null(dn[[2]])) dn[[2]][j] else NULL
        ))
    }
    
    # drop single dimensions
    if (drop && (length(i) == 1L || length(j) == 1L)) {
        return(as.vector(result))
    }
    
    result
}
