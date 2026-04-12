# S3_operations.R
#
# S3 generics and methods for:
#   multiply_sparse()   — sparse-aware matrix multiplication
#   split.HDF5Matrix()  — split dataset into blocks
#   reduce()            — reduce a group of datasets
#   apply_function()    — apply algebraic/statistical function to datasets
#
# Design notes:
#   - multiply_sparse():  new generic (no base R equivalent).
#   - split():            base::split(x, f, drop, ...) has "..."
#                         → no override needed.
#   - reduce():           new generic (no base R equivalent).
#   - apply_function():   new generic; name chosen to avoid collision with
#                         base::apply().


# ── multiply_sparse() — new generic ─────────────────────────────────────────

#' Sparse-aware matrix multiplication (generic)
#'
#' @description
#' Generic function for block-wise sparse matrix multiplication.
#' The method for \code{HDF5Matrix} computes \code{x \%*\% y} using the
#' BigDataStatMeth sparse multiplication algorithm, which skips all-zero
#' blocks and is more efficient when one or both matrices are highly sparse.
#'
#' @param x   An \code{HDF5Matrix}.
#' @param y   An \code{HDF5Matrix}. Must be in the same HDF5 file as \code{x}.
#' @param ... Additional arguments forwarded to the method.
#' @return A new \code{HDF5Matrix} containing the product.
#'
#' @examples
#' \donttest{
#' A <- hdf5_matrix("data.h5", "data/A")
#' B <- hdf5_matrix("data.h5", "data/B")
#' C <- multiply_sparse(A, B)
#' }
#'
#' @seealso \code{\link{multiply_sparse.HDF5Matrix}}
#'
#' @export
multiply_sparse <- function(x, y, ...) UseMethod("multiply_sparse")


#' Sparse-aware matrix multiplication for HDF5Matrix
#'
#' @description
#' Computes \code{x \%*\% y} block-wise using BigDataStatMeth's sparse
#' algorithm.
#'
#' @param x           An \code{HDF5Matrix}.
#' @param y           An \code{HDF5Matrix}. Same HDF5 file as \code{x}.
#' @param block_size  Integer. Block size hint; -1 = auto (default).
#' @param mix_block   Integer. Memory block size for parallel path; -1 = auto.
#' @param paral       Logical or NULL.
#' @param threads     Integer or NULL.
#' @param compression Integer (0-9) or NULL.
#' @param ...         Ignored.
#' @return A new \code{HDF5Matrix}.
#'
#' @exportS3Method
multiply_sparse.HDF5Matrix <- function(x, y,
                                        block_size  = -1L,
                                        mix_block   = -1L,
                                        paral       = NULL,
                                        threads     = NULL,
                                        compression = NULL,
                                        ...) {
    x$multiply_sparse(y,
                      block_size  = block_size,
                      mix_block   = mix_block,
                      paral       = paral,
                      threads     = threads,
                      compression = compression)
}


# ── split.HDF5Matrix ─────────────────────────────────────────────────────────

#' Split an HDF5Matrix into a list of blocks
#'
#' @description
#' S3 method of \code{base::split()} for \code{HDF5Matrix} objects.
#' Divides the matrix into blocks along rows (default) or columns.
#'
#' Provide exactly ONE of \code{n_blocks} or \code{block_size}.
#'
#' @param x           An \code{HDF5Matrix}.
#' @param f           Ignored (kept for S3 signature compatibility).
#' @param n_blocks    Integer. Number of (roughly equal) blocks; -1 = unused.
#' @param block_size  Integer. Max rows (or cols) per block; -1 = unused.
#' @param bycols      Logical. If \code{TRUE}, split by columns (default = by rows).
#' @param out_group   Character. HDF5 group for output blocks (default \code{"SPLIT"}).
#' @param out_dataset Character or NULL. Base dataset name.
#' @param overwrite   Logical. Overwrite existing blocks (default \code{FALSE}).
#' @param drop        Ignored (S3 compatibility).
#' @param ...         Ignored.
#' @return Named list of \code{HDF5Matrix} objects:
#'   \code{block_0}, \code{block_1}, …
#'
#' @examples
#' \donttest{
#' X      <- hdf5_matrix("data.h5", "data/X")   # 20000 × 1000
#' blocks <- split(X, n_blocks = 4)             # 4 row-blocks of ~5000 rows each
#' }
#'
#' @exportS3Method
split.HDF5Matrix <- function(x,
                              f            = NULL,
                              n_blocks     = -1L,
                              block_size   = -1L,
                              bycols       = FALSE,
                              out_group    = "SPLIT",
                              out_dataset  = NULL,
                              overwrite    = FALSE,
                              drop         = FALSE,
                              ...) {
    x$split(n_blocks   = n_blocks,
            block_size = block_size,
            bycols     = bycols,
            out_group  = out_group,
            out_dataset = out_dataset,
            overwrite  = overwrite)
}


# ── reduce() — new generic ───────────────────────────────────────────────────

#' Reduce a group of HDF5 datasets by accumulation (generic)
#'
#' @description
#' Generic function for reducing (accumulating) all datasets in the same HDF5
#' group as \code{x} into a single dataset using a binary operation.
#'
#' @param x   An \code{HDF5Matrix}.
#' @param ... Additional arguments forwarded to the method.
#' @return A new \code{HDF5Matrix} containing the accumulated result.
#'
#' @examples
#' \donttest{
#' partial <- hdf5_matrix("data.h5", "partials/chunk_0")
#' total   <- reduce(partial, func = "+")   # sums all datasets in "partials/"
#' }
#'
#' @seealso \code{\link{reduce.HDF5Matrix}}
#'
#' @export
reduce <- function(x, ...) UseMethod("reduce")


#' @exportS3Method
reduce.HDF5Matrix <- function(x,
                               out_group    = "REDUCED",
                               out_dataset  = "reduced",
                               func         = "+",
                               overwrite    = FALSE,
                               remove_input = FALSE,
                               ...) {
    x$reduce(out_group    = out_group,
             out_dataset  = out_dataset,
             func         = func,
             overwrite    = overwrite,
             remove_input = remove_input)
}


# ── apply_function() — new generic ──────────────────────────────────────────

#' Apply a statistical or algebraic function to HDF5 datasets (generic)
#'
#' @description
#' Generic function that applies one of BigDataStatMeth's algebraic or
#' statistical functions to a list of datasets in the same HDF5 group as
#' \code{x}.
#'
#' Valid \code{func} values: \code{"QR"}, \code{"CrossProd"},
#' \code{"tCrossProd"}, \code{"invChol"}, \code{"blockmult"},
#' \code{"CrossProd_double"}, \code{"tCrossProd_double"}, \code{"solve"},
#' \code{"normalize"}, \code{"sdmean"}, \code{"descChol"}.
#'
#' @param x   An \code{HDF5Matrix}.
#' @param ... Additional arguments forwarded to the method.
#' @return Named list with elements \code{filename}, \code{out_group},
#'   \code{func}, \code{datasets}.
#'
#' @examples
#' \donttest{
#' X   <- hdf5_matrix("data.h5", "data/A")
#' res <- apply_function(X, func = "CrossProd", out_group = "RESULTS")
#' }
#'
#' @seealso \code{\link{apply_function.HDF5Matrix}}
#'
#' @export
apply_function <- function(x, ...) UseMethod("apply_function")


#' @exportS3Method
apply_function.HDF5Matrix <- function(x,
                                       datasets    = NULL,
                                       func        = "QR",
                                       out_group   = "APPLIED",
                                       b_group     = "",
                                       b_datasets  = NULL,
                                       overwrite   = FALSE,
                                       transp_a    = FALSE,
                                       transp_b    = FALSE,
                                       full_matrix = FALSE,
                                       byrows      = FALSE,
                                       threads     = NULL,
                                       ...) {
    x$apply_function(
        datasets    = datasets,
        func        = func,
        out_group   = out_group,
        b_group     = b_group,
        b_datasets  = b_datasets,
        overwrite   = overwrite,
        transp_a    = transp_a,
        transp_b    = transp_b,
        full_matrix = full_matrix,
        byrows      = byrows,
        threads     = threads
    )
}
