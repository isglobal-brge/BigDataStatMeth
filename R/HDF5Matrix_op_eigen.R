# HDF5Matrix_op_eigen.R
#
# Adds $eigen() to the HDF5Matrix R6 class.
# Delegates to rcpp_hdf5dataset_eigen() -> RcppbdEigen_hdf5() (header).


HDF5Matrix$set("public", "eigen",
# @description
# Eigenvalue decomposition of a symmetric (or general) HDF5Matrix via
# Spectra (block-wise for large matrices, direct LAPACK for small ones).
#
# @param k               Integer. Number of eigenvalues/vectors to compute.
#   0 (default) = all.
# @param which           Character. Which eigenvalues to compute:
#   "LM" largest magnitude (default), "SM", "LR", "SR", "LI", "SI".
# @param ncv             Integer. Number of Lanczos basis vectors (0 = auto).
# @param symmetric       Logical. Treat matrix as symmetric (default TRUE).
# @param bcenter         Logical. Center columns before decomposition.
# @param bscale          Logical. Scale columns before decomposition.
# @param tolerance       Numeric. Convergence tolerance (default 1e-10).
# @param max_iter        Integer. Maximum Lanczos iterations (default 1000).
# @param compute_vectors Logical. Compute eigenvectors (default TRUE).
# @param overwrite       Logical. Overwrite existing results (default FALSE).
#   When TRUE, any open HDF5Matrix handles pointing to the output paths
#   (EIGEN/<dataset>/values, EIGEN/<dataset>/vectors) are automatically
#   closed and invalidated before the datasets are overwritten. This
#   prevents use-after-free crashes caused by stale open handles.
# @param threads         Integer. OpenMP threads (-1 = auto).
# @return Named list:
#   values       - numeric vector of eigenvalues (in memory),
#   vectors      - HDF5Matrix of eigenvectors, or NULL,
#   is_symmetric - logical.
function(k               = 0L,
         which           = "LM",
         ncv             = 0L,
         symmetric       = TRUE,
         bcenter         = FALSE,
         bscale          = FALSE,
         tolerance       = 1e-10,
         max_iter        = 1000L,
         compute_vectors = TRUE,
         overwrite       = FALSE,
         threads         = -1L) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")

    # When overwrite=TRUE, close any live HDF5Matrix handles that point to
    # the output paths. This prevents use-after-free crashes caused by a
    # previous $eigen() result (e.g. the $vectors component) still being
    # open when the underlying HDF5 datasets are deleted and recreated.
    if (isTRUE(overwrite)) {
        out_group <- paste0("EIGEN/", private$dataset)
        rcpp_hdf5_close_at_paths(
            private$filename,
            c(paste0(out_group, "/values"),
              paste0(out_group, "/vectors"))
        )
    }

    k_eff <- if (as.integer(k) == 0L) min(dim(self)) else as.integer(k)

    res <- rcpp_hdf5dataset_eigen(
        filename        = private$filename,
        group           = private$group,
        dataset         = private$dataset,
        k               = k_eff,
        which           = which,
        ncv             = as.integer(ncv),
        bcenter         = isTRUE(bcenter),
        bscale          = isTRUE(bscale),
        tolerance       = as.double(tolerance),
        max_iter        = as.integer(max_iter),
        compute_vectors = isTRUE(compute_vectors),
        overwrite       = isTRUE(overwrite),
        threads         = as.integer(threads)
    )

    # Eigenvalues are stored as a 1xk HDF5 dataset - read into memory
    vals_mat <- hdf5_matrix(res$file, res$path_values)
    values   <- as.numeric(as.matrix(vals_mat))
    close(vals_mat)

    vecs <- if (isTRUE(compute_vectors) && nzchar(res$path_vectors))
                hdf5_matrix(res$file, res$path_vectors)
            else NULL

    list(
        values       = values,
        vectors      = vecs,
        is_symmetric = isTRUE(symmetric)
    )
})
