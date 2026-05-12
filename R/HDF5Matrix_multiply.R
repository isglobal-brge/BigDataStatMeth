# HDF5Matrix - Algebra methods
# multiply (used by %*%, crossprod, tcrossprod S3 methods)
# Added via R6Class$set() to allow file partitioning


# @title Multiply two HDF5Matrix objects
#
# @description
# Computes matrix multiplication \code{A \%*\% B} (or transposed variants)
# using BigDataStatMeth's block-wise algorithm with optional parallelization.
#
# @param other HDF5Matrix object (right-hand side)
# @param transpose_self Logical. Transpose this matrix before multiplying?
# @param transpose_other Logical. Transpose \code{other} before multiplying?
# @param paral Logical or NULL. Enable OpenMP parallelization?
# @param block_size Integer or NULL. Number of elements per processing block.
# @param threads Integer or NULL. Number of OpenMP threads.
# @param compression Integer (0-9) or NULL. gzip compression level for the
#   result dataset. NULL uses the global option (default 6). Use 0 to disable.
#
# @return HDF5Matrix pointing to the result in group "OUTPUT",
#   dataset named "A_x_B" where A and B are the input dataset names.
HDF5Matrix$set("public", "multiply", function(other,
                                               transpose_self  = FALSE,
                                               transpose_other = FALSE,
                                               paral       = NULL,
                                               block_size  = NULL,
                                               threads     = NULL,
                                               compression = NULL) {
  # Validation
  if (!self$is_valid()) {
    stop("This dataset is closed or invalid")
  }
  if (!inherits(other, "HDF5Matrix")) {
    stop("other must be an HDF5Matrix object")
  }
  if (!other$is_valid()) {
    stop("other dataset is closed or invalid")
  }

  # Get effective parameter values (method args > global options > defaults)
  paral_eff       <- .get_option("paral",       default = NULL, override = paral)
  block_size_eff  <- .get_option("block_size",  default = NULL, override = block_size)
  threads_eff     <- .get_option("threads",     default = NULL, override = threads)
  compression_eff <- .get_option("compression", default = NULL,
                                             override = compression)

  # Get pointers
  other_ptr <- other$.__enclos_env__$private$ptr

  # Call C++ wrapper with effective parameters
  result_info <- rcpp_hdf5dataset_multiply(
    private$ptr,
    other_ptr,
    transpose_self,
    transpose_other,
    paral_eff,
    block_size_eff,
    threads_eff,
    compression_eff
  )

  # Return new HDF5Matrix wrapping the result
  hdf5_matrix(result_info$filename, result_info$path)
})
