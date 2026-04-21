# HDF5Matrix - I/O methods
# subset, read_all
# Added via R6Class$set() to allow file partitioning


# @description Read a rectangular subset of the dataset
#
# @param rows Integer vector of row indices (1-based)
# @param cols Integer vector of column indices (1-based)
# @return Numeric matrix with requested data
#
# @details
# Both contiguous ranges (\code{1:100}) and non-contiguous index vectors
# (\code{c(1,5,10)}) are supported. For standard R subsetting syntax
# (\code{X[i, j]}), see \code{\link{[.HDF5Matrix}}.
#
# @examples
# \donttest{
# tmp <- tempfile(fileext = ".h5")
# hdf5_create_matrix(tmp, "data/X", data = matrix(1:100, 10, 10))
# X <- hdf5_matrix(tmp, "data/matrix")
#
# X$subset(1:5, 1:3)
# X$subset(c(1, 3, 7), c(2, 5))
#
# X$close()
# unlink(tmp)
# }
HDF5Matrix$set("public", "subset", function(rows, cols) {
  if (!self$is_valid()) stop("Dataset is closed or invalid")

  if (!is.numeric(rows) && !is.integer(rows)) {
    stop("rows must be a numeric or integer vector")
  }
  if (!is.numeric(cols) && !is.integer(cols)) {
    stop("cols must be a numeric or integer vector")
  }

  rcpp_hdf5dataset_subset(private$ptr, as.integer(rows), as.integer(cols))
})


# @description Read the entire dataset into memory
#
# @return Numeric matrix with all data
#
# @details
# Equivalent to \code{X[, ]}. For very large datasets this may exhaust
# available RAM; prefer \code{$subset()} or \code{X[i, j]} for partial reads.
#
# @examples
# \donttest{
# tmp <- tempfile(fileext = ".h5")
# hdf5_create_matrix(tmp, "data/X", data = matrix(1:100, 10, 10))
# X <- hdf5_matrix(tmp, "data/matrix")
#
# all_data <- X$read_all()
# dim(all_data)  # 10 x 10
#
# X$close()
# unlink(tmp)
# }
HDF5Matrix$set("public", "read_all", function() {
  if (!self$is_valid()) stop("Dataset is closed or invalid")
  rcpp_hdf5dataset_read_all(private$ptr)
})
