# HDF5Matrix_op_pseudoinv.R
#
# Adds $pseudoinverse() to the HDF5Matrix R6 class.
# Delegates to rcpp_hdf5dataset_pseudoinv() -> RcppPseudoinvHdf5() (header).


HDF5Matrix$set("public", "pseudoinverse",
# @description
# Moore-Penrose pseudoinverse of an HDF5Matrix via SVD.
# For A (n x m), returns A+ (m x n).
#
# @param out_group   Character. Output HDF5 group (default "PseudoInverse").
# @param out_dataset Character or NULL. Output dataset name (default = input name).
# @param overwrite   Logical. Overwrite existing result (default FALSE).
# @param threads     Integer. OpenMP threads (-1 = auto).
# @param compression Integer (0-9) or NULL. gzip level; NULL inherits from input.
# @return A new HDF5Matrix containing the pseudoinverse (m x n).
function(out_group   = "PseudoInverse",
         out_dataset = NULL,
         overwrite   = FALSE,
         threads     = -1L,
         compression = NULL) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    compression_eff <- .get_option("compression", default = NULL,
                                   override = compression)

    rcpp_hdf5dataset_close(private$ptr)
    private$ptr <- NULL
    on.exit({
        if (is.null(private$ptr))
            private$ptr <- rcpp_hdf5dataset_open(
                private$filename, private$group, private$dataset)
    }, add = TRUE)

    res <- rcpp_hdf5dataset_pseudoinv(
        filename    = private$filename,
        group       = private$group,
        dataset     = private$dataset,
        out_group   = out_group,
        out_dataset = if (is.null(out_dataset)) "" else out_dataset,
        overwrite   = isTRUE(overwrite),
        threads     = as.integer(threads),
        compression = if (is.null(compression_eff)) -1L
                      else as.integer(compression_eff)
    )

    hdf5_matrix(res$file, res$path)
})
