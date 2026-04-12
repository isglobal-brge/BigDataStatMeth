# HDF5Matrix - Elementwise arithmetic methods
# add, subtract, multiply_ew, divide_ew
# Added via R6Class$set() following the same pattern as HDF5Matrix_multiply.R


# @description Element-wise addition: C = A + B
#
# @param other HDF5Matrix (right-hand side)
# @param paral Logical or NULL. Enable OpenMP parallelisation.
# @param block_size Integer or NULL. Block size for processing.
# @param threads Integer or NULL. Number of OpenMP threads.
# @return New HDF5Matrix containing the result.
HDF5Matrix$set("public", "add", function(other,
                                          paral       = NULL,
                                          block_size  = NULL,
                                          threads     = NULL,
                                          compression = NULL) {
    if (!self$is_valid())          stop("This dataset is closed or invalid")
    if (!inherits(other, "HDF5Matrix")) stop("other must be an HDF5Matrix object")
    if (!other$is_valid())         stop("other dataset is closed or invalid")

    paral_eff       <- .get_option("paral",       default = NULL, override = paral)
    block_size_eff  <- .get_option("block_size",  default = NULL, override = block_size)
    threads_eff     <- .get_option("threads",     default = NULL, override = threads)
    compression_eff <- .get_option("compression", default = NULL, override = compression)

    result_info <- rcpp_hdf5dataset_add(
        private$ptr,
        other$.__enclos_env__$private$ptr,
        paral_eff, block_size_eff, threads_eff, compression_eff
    )
    hdf5_matrix(result_info$filename, result_info$path)
})


# @description Element-wise subtraction: C = A - B
#
# @param other HDF5Matrix (right-hand side)
# @param paral Logical or NULL. Enable OpenMP parallelisation.
# @param block_size Integer or NULL. Block size for processing.
# @param threads Integer or NULL. Number of OpenMP threads.
# @return New HDF5Matrix containing the result.
HDF5Matrix$set("public", "subtract", function(other,
                                               paral       = NULL,
                                               block_size  = NULL,
                                               threads     = NULL,
                                               compression = NULL) {
    if (!self$is_valid())          stop("This dataset is closed or invalid")
    if (!inherits(other, "HDF5Matrix")) stop("other must be an HDF5Matrix object")
    if (!other$is_valid())         stop("other dataset is closed or invalid")

    paral_eff       <- .get_option("paral",       default = NULL, override = paral)
    block_size_eff  <- .get_option("block_size",  default = NULL, override = block_size)
    threads_eff     <- .get_option("threads",     default = NULL, override = threads)
    compression_eff <- .get_option("compression", default = NULL, override = compression)

    result_info <- rcpp_hdf5dataset_subtract(
        private$ptr,
        other$.__enclos_env__$private$ptr,
        paral_eff, block_size_eff, threads_eff, compression_eff
    )
    hdf5_matrix(result_info$filename, result_info$path)
})


# @description Element-wise multiplication (Hadamard product): C = A * B
#
# @param other HDF5Matrix (right-hand side)
# @param paral Logical or NULL. Enable OpenMP parallelisation.
# @param block_size Integer or NULL. Block size for processing.
# @param threads Integer or NULL. Number of OpenMP threads.
# @return New HDF5Matrix containing the result.
HDF5Matrix$set("public", "multiply_ew", function(other,
                                                   paral       = NULL,
                                                   block_size  = NULL,
                                                   threads     = NULL,
                                                   compression = NULL) {
    if (!self$is_valid())          stop("This dataset is closed or invalid")
    if (!inherits(other, "HDF5Matrix")) stop("other must be an HDF5Matrix object")
    if (!other$is_valid())         stop("other dataset is closed or invalid")

    paral_eff       <- .get_option("paral",       default = NULL, override = paral)
    block_size_eff  <- .get_option("block_size",  default = NULL, override = block_size)
    threads_eff     <- .get_option("threads",     default = NULL, override = threads)
    compression_eff <- .get_option("compression", default = NULL, override = compression)

    result_info <- rcpp_hdf5dataset_mul_ew(
        private$ptr,
        other$.__enclos_env__$private$ptr,
        paral_eff, block_size_eff, threads_eff, compression_eff
    )
    hdf5_matrix(result_info$filename, result_info$path)
})


# @description Element-wise division: C = A / B
#
# Division by zero produces NaN or Inf, matching base R behaviour.
#
# @param other HDF5Matrix (right-hand side)
# @param paral Logical or NULL. Enable OpenMP parallelisation.
# @param block_size Integer or NULL. Block size for processing.
# @param threads Integer or NULL. Number of OpenMP threads.
# @return New HDF5Matrix containing the result.
HDF5Matrix$set("public", "divide_ew", function(other,
                                                paral       = NULL,
                                                block_size  = NULL,
                                                threads     = NULL,
                                                compression = NULL) {
    if (!self$is_valid())          stop("This dataset is closed or invalid")
    if (!inherits(other, "HDF5Matrix")) stop("other must be an HDF5Matrix object")
    if (!other$is_valid())         stop("other dataset is closed or invalid")

    paral_eff       <- .get_option("paral",       default = NULL, override = paral)
    block_size_eff  <- .get_option("block_size",  default = NULL, override = block_size)
    threads_eff     <- .get_option("threads",     default = NULL, override = threads)
    compression_eff <- .get_option("compression", default = NULL, override = compression)

    result_info <- rcpp_hdf5dataset_div_ew(
        private$ptr,
        other$.__enclos_env__$private$ptr,
        paral_eff, block_size_eff, threads_eff, compression_eff
    )
    hdf5_matrix(result_info$filename, result_info$path)
})
