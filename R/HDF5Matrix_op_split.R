# HDF5Matrix_op_split.R
#
# Adds $split() to the HDF5Matrix R6 class.
# Delegates to rcpp_hdf5dataset_split() -> RcppSplit_matrix_hdf5() (header).


HDF5Matrix$set("public", "split",
# @description
# Splits this HDF5Matrix into a list of HDF5Matrix blocks along rows
# (bycols = FALSE) or columns (bycols = TRUE).
#
# Provide exactly ONE of n_blocks or block_size:
#   n_blocks:   split into N approximately equal pieces.
#   block_size: each piece has at most this many rows/columns.
#
# Output blocks are stored in the same HDF5 file at
#   <out_group>/<out_dataset>.0, <out_group>/<out_dataset>.1, ...
#
# @param n_blocks    Integer. Number of blocks (-1 = use block_size).
# @param block_size  Integer. Max rows (or cols) per block (-1 = use n_blocks).
# @param bycols      Logical. Split by columns (default FALSE = by rows).
# @param out_group   Character. HDF5 group for output (default "SPLIT").
# @param out_dataset Character or NULL. Base dataset name (default = input name).
# @param overwrite   Logical. Overwrite existing blocks (default FALSE).
# @return A named list of HDF5Matrix objects: block_0, block_1, ...
function(n_blocks    = -1L,
         block_size  = -1L,
         bycols      = FALSE,
         out_group   = "SPLIT",
         out_dataset = NULL,
         overwrite   = FALSE) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    if (n_blocks == -1L && block_size == -1L)
        stop("split: provide one of n_blocks or block_size")
    if (n_blocks > 0L && block_size > 0L)
        stop("split: provide exactly one of n_blocks or block_size, not both")

    ods <- if (is.null(out_dataset)) private$dataset else out_dataset

    # C++ opens independent handles; no need to close own ptr
    res <- rcpp_hdf5dataset_split(
        ptr         = private$ptr,
        bycols      = isTRUE(bycols),
        n_blocks    = as.integer(n_blocks),
        block_size  = as.integer(block_size),
        out_group   = out_group,
        out_dataset = ods,
        overwrite   = isTRUE(overwrite)
    )

    # Build list of HDF5Matrix objects, one per block
    nb   <- res$n_blocks
    base <- res$out_dataset
    grp  <- res$out_group
    fn   <- res$filename

    out <- vector("list", nb)
    for (i in seq_len(nb)) {
        path     <- paste0(grp, "/", base, ".", i - 1L)
        out[[i]] <- hdf5_matrix(fn, path)
    }
    names(out) <- paste0("block_", seq_len(nb) - 1L)
    out
})
