# HDF5Matrix - Aggregation methods
# colSums, rowSums, colMeans, rowMeans, colMins, rowMins,
# colMaxs, rowMaxs, colVars, rowVars, colSds, rowSds,
# sum, mean, min, max, var, sd
# Added via HDF5Matrix$set() following the same pattern as the other
# HDF5Matrix_op_*.R files.


# ---------------------------------------------------------------------------
# Internal helper shared by all aggregation methods
# ---------------------------------------------------------------------------

# Resolve paral/threads from method arguments and package-level options.
# Returns a list(paral, threads) ready to pass to the C++ wrappers.
.agg_options <- function(paral, threads) {
    list(
        paral   = .get_option("paral",   default = NULL, override = paral),
        threads = .get_option("threads", default = NULL, override = threads)
    )
}


# ---------------------------------------------------------------------------
# Column-wise methods  (result length = ncols_R)
# ---------------------------------------------------------------------------

# @description Column sums: equivalent to colSums(X).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length ncol(X).
HDF5Matrix$set("public", "colSums", function(paral   = NULL,
                                              wsize   = NULL,
                                              threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_colSums(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Column means: equivalent to colMeans(X).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length ncol(X).
HDF5Matrix$set("public", "colMeans", function(paral   = NULL,
                                               wsize   = NULL,
                                               threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_colMeans(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Column minimums: equivalent to apply(X, 2, min).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length ncol(X).
HDF5Matrix$set("public", "colMins", function(paral   = NULL,
                                              wsize   = NULL,
                                              threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_colMins(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Column maximums: equivalent to apply(X, 2, max).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length ncol(X).
HDF5Matrix$set("public", "colMaxs", function(paral   = NULL,
                                              wsize   = NULL,
                                              threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_colMaxs(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Column variances: equivalent to apply(X, 2, var).
# Uses Bessel's correction (n-1). Returns NaN for columns with < 2 observations.
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length ncol(X).
HDF5Matrix$set("public", "colVars", function(paral   = NULL,
                                              wsize   = NULL,
                                              threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_colVars(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Column standard deviations: equivalent to apply(X, 2, sd).
# Uses Bessel's correction (n-1).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length ncol(X).
HDF5Matrix$set("public", "colSds", function(paral   = NULL,
                                             wsize   = NULL,
                                             threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_colSds(private$ptr, opt$paral, wsize, opt$threads)
})


# ---------------------------------------------------------------------------
# Row-wise methods  (result length = nrows_R)
# ---------------------------------------------------------------------------

# @description Row sums: equivalent to rowSums(X).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length nrow(X).
HDF5Matrix$set("public", "rowSums", function(paral   = NULL,
                                              wsize   = NULL,
                                              threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_rowSums(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Row means: equivalent to rowMeans(X).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length nrow(X).
HDF5Matrix$set("public", "rowMeans", function(paral   = NULL,
                                               wsize   = NULL,
                                               threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_rowMeans(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Row minimums: equivalent to apply(X, 1, min).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length nrow(X).
HDF5Matrix$set("public", "rowMins", function(paral   = NULL,
                                              wsize   = NULL,
                                              threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_rowMins(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Row maximums: equivalent to apply(X, 1, max).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length nrow(X).
HDF5Matrix$set("public", "rowMaxs", function(paral   = NULL,
                                              wsize   = NULL,
                                              threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_rowMaxs(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Row variances: equivalent to apply(X, 1, var).
# Uses Bessel's correction (n-1). Returns NaN for rows with < 2 observations.
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length nrow(X).
HDF5Matrix$set("public", "rowVars", function(paral   = NULL,
                                              wsize   = NULL,
                                              threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_rowVars(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Row standard deviations: equivalent to apply(X, 1, sd).
# Uses Bessel's correction (n-1).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Named numeric vector of length nrow(X).
HDF5Matrix$set("public", "rowSds", function(paral   = NULL,
                                             wsize   = NULL,
                                             threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_rowSds(private$ptr, opt$paral, wsize, opt$threads)
})


# ---------------------------------------------------------------------------
# Scalar methods  (whole matrix treated as a flat vector)
# ---------------------------------------------------------------------------

# @description Sum of all elements: equivalent to sum(X).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Scalar numeric.
HDF5Matrix$set("public", "sum", function(paral   = NULL,
                                          wsize   = NULL,
                                          threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_scalar_sum(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Mean of all elements: equivalent to mean(X).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Scalar numeric.
HDF5Matrix$set("public", "mean", function(paral   = NULL,
                                           wsize   = NULL,
                                           threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_scalar_mean(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Minimum of all elements: equivalent to min(X).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Scalar numeric.
HDF5Matrix$set("public", "min", function(paral   = NULL,
                                          wsize   = NULL,
                                          threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_scalar_min(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Maximum of all elements: equivalent to max(X).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Scalar numeric.
HDF5Matrix$set("public", "max", function(paral   = NULL,
                                          wsize   = NULL,
                                          threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_scalar_max(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Variance of all elements: equivalent to var(as.vector(X)).
# Uses Bessel's correction (N-1). Returns NaN if total elements < 2.
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Scalar numeric.
HDF5Matrix$set("public", "var", function(paral   = NULL,
                                          wsize   = NULL,
                                          threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_scalar_var(private$ptr, opt$paral, wsize, opt$threads)
})


# @description Standard deviation of all elements: equivalent to sd(as.vector(X)).
# Uses Bessel's correction (N-1).
#
# @param paral   Logical or NULL. Enable OpenMP parallelisation.
# @param wsize   Integer or NULL. Block size for HDF5 reads (NULL = auto).
# @param threads Integer or NULL. Number of OpenMP threads.
# @return Scalar numeric.
HDF5Matrix$set("public", "sd", function(paral   = NULL,
                                         wsize   = NULL,
                                         threads = NULL) {
    if (!self$is_valid()) stop("Dataset is closed or invalid")
    opt <- .agg_options(paral, threads)
    rcpp_hdf5dataset_scalar_sd(private$ptr, opt$paral, wsize, opt$threads)
})
