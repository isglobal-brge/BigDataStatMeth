# S3_split.R
#
# S3 generic split_dataset() and HDF5Matrix method.
# For the base::split() S3 method see S3_operations.R.

#' Split an HDF5Matrix into multiple block datasets
#'
#' @description
#' Splits an \code{\link{HDF5Matrix}} into equal-sized sub-matrices
#' stored as separate datasets in the same HDF5 file.
#'
#' Output datasets are named \code{<out_group>/<out_dataset>.0},
#' \code{<out_group>/<out_dataset>.1}, \ldots (0-based index).
#'
#' Exactly one of \code{n_blocks} or \code{block_size} must be provided.
#'
#' @param x         An \code{HDF5Matrix}.
#' @param n_blocks  Integer or \code{NULL}. Number of blocks.
#' @param block_size Integer or \code{NULL}. Rows or columns per block.
#' @param bycols    Logical. Split by columns (\code{TRUE}) or rows (default \code{FALSE}).
#' @param out_group Character. Output HDF5 group (default \code{"SPLIT"}).
#' @param out_dataset Character or NULL. Base dataset name.
#' @param overwrite Logical. Overwrite existing blocks (default \code{FALSE}).
#' @param \dots     Ignored.
#' @return A named list of \code{HDF5Matrix} objects.
#'
#' @examples
#' \donttest{
#' tmp  <- tempfile(fileext = ".h5")
#' M    <- hdf5_create_matrix(tmp, "data/M", data = matrix(1:60, 6, 10))
#' blks <- split_dataset(M, n_blocks = 3L)
#' length(blks)
#' lapply(blks, close)
#' close(M)
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link[BigDataStatMeth]{cbind.HDF5Matrix}}
#' @export
split_dataset <- function(x, n_blocks = NULL, block_size = NULL,
                          bycols = FALSE, ...) UseMethod("split_dataset")

#' @rdname split_dataset
#' @export
split_dataset.HDF5Matrix <- function(x, n_blocks = NULL, block_size = NULL,
                                     bycols = FALSE,
                                     out_group   = "SPLIT",
                                     out_dataset = NULL,
                                     overwrite   = FALSE,
                                     ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")
    n  <- if (is.null(n_blocks))   -1L else as.integer(n_blocks)
    bs <- if (is.null(block_size)) -1L else as.integer(block_size)
    x$split(n_blocks    = n,
            block_size  = bs,
            bycols      = bycols,
            out_group   = out_group,
            out_dataset = out_dataset,
            overwrite   = overwrite)
}
