#' BigDataStatMeth: Scalable statistical computing with R, C++, and HDF5
#'
#' @description
#' BigDataStatMeth provides statistical and linear algebra operations for
#' matrices stored in HDF5 files. The package is designed for workflows in
#' which matrices may be too large to be held entirely in memory, while
#' still allowing users to work with familiar R functions.
#'
#' The recommended user-facing interface is based on \code{HDF5Matrix}
#' objects and standard R methods. HDF5-backed matrices can be manipulated
#' using calls such as \code{dim()}, \code{[}, \code{\%*\%},
#' \code{crossprod()}, \code{tcrossprod()}, \code{scale()},
#' \code{cor()}, \code{svd()}, \code{prcomp()}, \code{qr()},
#' \code{chol()}, and \code{solve()}.
#'
#' @section Main user-facing functionality:
#' \itemize{
#'   \item Core HDF5 matrix handling: \code{hdf5_create_matrix()},
#'         \code{hdf5_matrix()}, \code{list_datasets()}, \code{is_open()},
#'         \code{close()}, and \code{hdf5_close_all()}.
#'   \item Subsetting and conversion: \code{[}, \code{[<-},
#'         \code{as.matrix()}, and \code{as.data.frame()}.
#'   \item Dimension names: \code{rownames()}, \code{colnames()}, and
#'         \code{dimnames()}.
#'   \item Element-wise arithmetic: \code{+}, \code{-}, \code{*}, and
#'         \code{/} for \code{HDF5Matrix} objects.
#'   \item Matrix algebra: \code{\%*\%}, \code{crossprod()},
#'         \code{tcrossprod()}, \code{cbind()}, and \code{rbind()}.
#'   \item Aggregations and summaries: \code{colSums()},
#'         \code{rowSums()}, \code{colMeans()}, \code{rowMeans()},
#'         \code{colVars()}, \code{rowVars()}, \code{colSds()},
#'         \code{rowSds()}, \code{colMins()}, \code{rowMins()},
#'         \code{colMaxs()}, \code{rowMaxs()}, \code{mean()},
#'         \code{var()}, and \code{sd()}.
#'   \item Statistical transformations: \code{scale()}, \code{sweep()},
#'         and \code{cor()}.
#'   \item Matrix decompositions and factorizations: \code{svd()},
#'         \code{prcomp()}, \code{qr()}, \code{chol()}, \code{solve()},
#'         \code{eigen()}, and \code{pseudoinverse()}.
#'   \item Diagonal, split, reduce, and apply operations:
#'         \code{diag()}, \code{diag_op()}, \code{diag_scale()},
#'         \code{split_dataset()}, \code{reduce()}, and
#'         \code{apply_function()}.
#' }
#'
#' @section Additional high-level utilities:
#' Most user workflows can be expressed through \code{HDF5Matrix} objects
#' and standard R methods. Some functions keep the \code{bd*} prefix
#' because they provide additional utilities that do not map directly to a
#' standard R generic, or because they expose workflows available in earlier
#' versions of the package. Examples include utilities for creating HDF5
#' groups, moving datasets, and writing HDF5-backed dimension names. These
#' functions remain part of the package API and are documented in their
#' corresponding help pages.
#'
#' @section Global options and HDF5 resources:
#' Block-wise operations can be configured with
#' \code{hdf5matrix_options()}, including options for parallel execution,
#' number of threads, block size, and HDF5 compression. Open HDF5 resources
#' can be closed explicitly with \code{close()} for individual objects or
#' \code{hdf5_close_all()} for all handles tracked by the package.
#'
#' @section Architecture and developer interfaces:
#' BigDataStatMeth is organized around a standard R interface backed by
#' a C++ computational infrastructure. The user-facing layer is based on
#' \code{HDF5Matrix} objects and S3 methods, allowing HDF5-backed
#' matrices to be used with familiar R functions.
#'
#' Internally, a lightweight R6 layer connects these R methods with the
#' C++ backend. The C++ infrastructure provides classes for managing
#' HDF5 files, groups, and datasets, together with block-wise routines
#' for linear algebra and statistical operations.
#'
#' This design allows developers to implement new scalable methods from
#' Rcpp-based code while reusing the package machinery for HDF5 file
#' management, block iteration, compression handling, and numerical
#' computation.
#'
#' @section Getting started:
#' See \code{vignette("BigDataStatMeth")} for a practical introduction to
#' HDF5-backed matrices and the main user-facing functionality.
#'
#' @examples
#' h5file <- tempfile(fileext = ".h5")
#'
#' set.seed(1)
#' X <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#'
#' X_h5 <- hdf5_create_matrix(
#'   filename = h5file,
#'   dataset = "data/X",
#'   data = X,
#'   overwrite = TRUE
#' )
#'
#' dim(X_h5)
#' colMeans(X_h5)
#'
#' XtX_h5 <- crossprod(X_h5)
#' dim(XtX_h5)
#'
#' close(X_h5)
#' close(XtX_h5)
#' hdf5_close_all(verbose = FALSE)
#'
#' @name BigDataStatMeth
#' @useDynLib BigDataStatMeth
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom utils download.file unzip untar
#' @importFrom data.table "%like%"
#' @importFrom RCurl url.exists
NULL