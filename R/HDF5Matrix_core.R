# HDF5Matrix - Core methods
# dim, info, is_valid, print, close
# Added via R6Class$set() to allow file partitioning


# @description Get dimensions of dataset
# @return Integer vector \code{c(nrows, ncols)}
#
# @examples
# \donttest{
# tmp <- tempfile(fileext = ".h5")
# hdf5_create_matrix(tmp, "data/X", data = matrix(1:100, 10, 10))
# X <- hdf5_matrix(tmp, "data/matrix")
# X$dim()
# X$close()
# unlink(tmp)
# }
HDF5Matrix$set("public", "dim", function() {
  if (!self$is_valid()) stop("Dataset is closed or invalid")
  rcpp_hdf5dataset_dim(private$ptr)
})


# @description Get dataset information
# @return List with \code{filename}, \code{group}, \code{dataset},
#   \code{type}, and \code{is_open}
#
# @examples
# \donttest{
# tmp <- tempfile(fileext = ".h5")
# hdf5_create_matrix(tmp, "data/X", data = matrix(1:100, 10, 10))
# X <- hdf5_matrix(tmp, "data/matrix")
# X$info()
# X$close()
# unlink(tmp)
# }
HDF5Matrix$set("public", "info", function() {
  if (!self$is_valid()) stop("Dataset is closed or invalid")
  rcpp_hdf5dataset_info(private$ptr)
})


# @description Check if dataset is valid and open
# @return Logical
HDF5Matrix$set("public", "is_valid", function() {
  if (is.null(private$ptr)) return(FALSE)
  rcpp_hdf5dataset_is_valid(private$ptr)
})


# @description Print HDF5Matrix summary
# @return Invisible self
#
# @examples
# \donttest{
# tmp <- tempfile(fileext = ".h5")
# hdf5_create_matrix(tmp, "data/X", data = matrix(1:100, 10, 10))
# X <- hdf5_matrix(tmp, "data/matrix")
# print(X)
# X$close()
# unlink(tmp)
# }
HDF5Matrix$set("public", "print", function() {
  cat("HDF5Matrix object\n")
  cat("  File: ",  private$filename, "\n", sep = "")
  cat("  Path: ",  private$group, "/", private$dataset, "\n", sep = "")
  if (self$is_valid()) {
    dims <- self$dim()
    cat("  Dimensions: ", dims[1], " x ", dims[2], "\n", sep = "")
    info <- self$info()
    cat("  Type: ",   info$type, "\n", sep = "")
    cat("  Status: OPEN\n")
  } else {
    cat("  Status: CLOSED\n")
  }
  invisible(self)
})


# @description Close dataset and release file lock immediately
#
# @details
# Closes the HDF5 dataset without waiting for garbage collection.
# After calling \code{close()}, \code{is_valid()} returns \code{FALSE} and
# any further operations on the object will fail. The file is immediately
# accessible by other tools such as HDFView.
#
# For emergency closure of all open objects, see
# \code{\link{hdf5_close_all}}.
#
# @return Invisible self
#
# @examples
# \donttest{
# tmp <- tempfile(fileext = ".h5")
# hdf5_create_matrix(tmp, "data/X", data = matrix(1:100, 10, 10))
# X <- hdf5_matrix(tmp, "data/matrix")
# data <- X[1:5, 1:5]
#
# X$close()
# X$is_valid()  # FALSE
#
# unlink(tmp)
# }
# ── Public path getters ───────────────────────────────────────────────────────
# These allow S3 methods and user code to retrieve the HDF5 location of an
# object without having to call the heavier $info() which opens the C++ handle.

# @description Return the HDF5 file path for this dataset.
# @return Character string.
HDF5Matrix$set("public", "get_filename", function() private$filename)

# @description Return the HDF5 group path for this dataset.
# @return Character string.
HDF5Matrix$set("public", "get_group", function() private$group)

# @description Return the dataset name within its group.
# @return Character string.
HDF5Matrix$set("public", "get_dataset", function() private$dataset)

# @description Return the full HDF5 path \code{"group/dataset"}.
# @return Character string.
HDF5Matrix$set("public", "get_path", function()
    paste(private$group, private$dataset, sep = "/"))


HDF5Matrix$set("public", "close", function() {
  if (!is.null(private$ptr)) {
      rcpp_hdf5dataset_close(private$ptr)
      private$ptr <- NULL
      invisible(gc())
        # message("HDF5 dataset closed")
  } # else {
  #   message("Dataset already closed")
  # }
  invisible(self)
})
