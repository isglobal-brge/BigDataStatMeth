# HDF5Matrix_op_sparse.R
#
# Adds $multiply_sparse() to the HDF5Matrix R6 class.
# Delegates to rcpp_hdf5dataset_multiply_sparse() -> multiplicationSparse() (header).


HDF5Matrix$set("public", "multiply_sparse",
# @description
# Sparse-aware block-wise matrix multiplication A %*% B, where one
# or both matrices may have many zero blocks.  Uses BigDataStatMeth's
# sparse multiplication algorithm which skips all-zero blocks.
# Both matrices must be in the same HDF5 file.
#
# @param other      An HDF5Matrix. Must be in the same HDF5 file.
# @param outgroup  Character or \code{NULL}. Output group in the HDF5 file.
#   Default \code{"OUTPUT"}.
# @param outdataset Character or \code{NULL}. Output dataset name.
#   Default is constructed from the operation and input nam
# @param block_size Integer. Block size hint; -1 = auto (recommended).
# @param mix_block  Integer. Memory block size for parallel path; -1 = auto.
# @param paral      Logical or NULL. Enable OpenMP parallelisation.
# @param threads    Integer or NULL. Thread count.
# @param compression Integer (0-9) or NULL. gzip level; NULL inherits from self.
# @return A new HDF5Matrix with the product A x B.
function(other,
         outgroup    = NULL,
         outdataset  = NULL,
         block_size  = -1L,
         mix_block   = -1L,
         paral       = NULL,
         threads     = NULL,
         compression = NULL) {

    if (!self$is_valid())         stop("HDF5Matrix is closed or invalid")
    if (!inherits(other, "HDF5Matrix")) stop("other must be an HDF5Matrix")
    if (!other$is_valid())        stop("other is closed or invalid")
    if (private$filename != other$get_filename())
        stop("multiply_sparse: both objects must be in the same HDF5 file")

    paral_eff       <- .get_option("paral",       default = NULL, override = paral)
    threads_eff     <- .get_option("threads",     default = NULL, override = threads)
    compression_eff <- .get_option("compression", default = NULL, override = compression)

    other_ptr <- other$.__enclos_env__$private$ptr

    res <- rcpp_hdf5dataset_multiply_sparse(
        ptr_a       = private$ptr,
        ptr_b       = other_ptr,
        block_size  = as.integer(block_size),
        mix_block   = as.integer(mix_block),
        paral       = paral_eff,
        threads     = threads_eff,
        compression = compression_eff,
        outgroup   = outgroup,
        outdataset = outdataset
    )

    hdf5_matrix(res$filename, res$path)
})
