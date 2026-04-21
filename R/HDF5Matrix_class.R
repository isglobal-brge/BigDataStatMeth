#' HDF5Matrix R6 Class
#'
#' @description
#' R6 class for persistent access to HDF5 datasets with automatic cleanup.
#' Provides a familiar matrix-like interface for large datasets stored in HDF5
#' files, using block-wise algorithms and optional OpenMP parallelization.
#'
#' @details
#' This class wraps a C++ \code{hdf5Dataset} object using an external pointer.
#' When the R object is garbage collected, the C++ destructor is called
#' automatically, ensuring proper cleanup of HDF5 resources.
#'
#' Users should create objects via the \code{\link{hdf5_matrix}} constructor,
#' not directly through \code{HDF5Matrix$new()}.
#'
#' @import R6
#' @keywords internal
#' @noRd
HDF5Matrix <- R6::R6Class("HDF5Matrix",

  private = list(
    ptr      = NULL,  # External pointer to C++ hdf5Dataset
    filename = NULL,  # HDF5 file path
    group    = NULL,  # Group path within HDF5 file
    dataset  = NULL,  # Dataset name within group
    sentinel = NULL,

    # Finalizer - called by R garbage collector
    finalize = function() {
        if (!is.null(private$ptr)) {
            tryCatch(rcpp_hdf5dataset_close(private$ptr), error = function(e) NULL)
            private$ptr <- NULL
        }
    }
    # finalize = function() {
    #     if (!is.null(private$ptr)) {
    #         # rcpp_hdf5dataset_close(private$ptr)
    #         tryCatch(rcpp_hdf5dataset_close(private$ptr), error = function(e) NULL)
    #         private$ptr <- NULL
    #     }
    # }
  ),

  public = list(

    #' @description
    #' Create a new HDF5Matrix object
    #'
    #' @param filename Path to HDF5 file
    #' @param path Full dataset path (e.g., \code{"group/dataset"})
    #'
    #' @examples
    #' \donttest{
    #' tmp <- tempfile(fileext = ".h5")
    #'
    #' X <- HDF5Matrix$new(tmp, "data/matrix")
    #' dim(X)
    #'
    #' X$close()
    #' unlink(tmp)
    #' }
    initialize = function(filename, path) {

      if (!is.character(filename) || length(filename) != 1) {
        stop("filename must be a single string")
      }
      if (!is.character(path) || length(path) != 1) {
        stop("path must be a single string")
      }
      if (!file.exists(filename)) {
        stop("File does not exist: ", filename)
      }

      parts <- strsplit(path, "/")[[1]]
      if (length(parts) < 2) {
        stop("path must be in format 'group/dataset' or 'group/subgroup/dataset'")
      }

      private$filename <- filename
      private$dataset  <- parts[length(parts)]
      private$group    <- paste(parts[-length(parts)], collapse = "/")

      private$ptr <- rcpp_hdf5dataset_open(
        private$filename,
        private$group,
        private$dataset
      )
      
      private$sentinel <- new.env(parent = emptyenv())
      reg.finalizer(private$sentinel, function(e) {
          tryCatch({
              if (!is.null(private$ptr)) {
                  rcpp_hdf5dataset_close(private$ptr)
                  private$ptr <- NULL
              }
          }, error = function(x) NULL)
      }, onexit = TRUE)
      
      
      # reg.finalizer(private$sentinel, function(e) {
      #     if (!is.null(private$ptr)) {
      #         tryCatch(rcpp_hdf5dataset_close(private$ptr), error = function(x) NULL)
      #         private$ptr <- NULL
      #     }
      # }, onexit = TRUE)
    }
  )
)
