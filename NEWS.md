# BigDataStatMeth 2.0.0

This is a major release of BigDataStatMeth. It introduces a new
`HDF5Matrix` user-facing interface, standard S3 methods for
HDF5-backed matrices, and an extended C++ infrastructure for block-wise
statistical computing with HDF5 files.

## Major changes

* Added the `HDF5Matrix` interface for working with matrices stored in
  HDF5 files.
* Added S3 methods that allow HDF5-backed matrices to be used with
  familiar R calls such as `dim()`, `[`, `[<-`, `%*%`, `crossprod()`,
  `tcrossprod()`, `scale()`, `cor()`, `svd()`, `prcomp()`, `qr()`,
  `chol()`, and `solve()`.
* Reorganized the package around a standard R interface backed by a C++
  computational infrastructure for block-wise statistical computing.
* Added global configuration through `hdf5matrix_options()` for common
  settings such as parallel execution, number of threads, block size,
  and HDF5 compression.

## User-facing functionality

* Added or extended support for creating, opening, inspecting, and
  closing HDF5-backed matrices.
* Added support for subsetting, assignment, dimension names, and
  conversion to in-memory R objects.
* Added S3 methods for element-wise arithmetic operations on
  `HDF5Matrix` objects.
* Added S3 support for matrix algebra operations, including `%*%`,
  `crossprod()`, `tcrossprod()`, `cbind()`, and `rbind()`.
* Added aggregation and summary methods, including `colSums()`,
  `rowSums()`, `colMeans()`, `rowMeans()`, `colVars()`, `rowVars()`,
  `colSds()`, `rowSds()`, `colMins()`, `rowMins()`, `colMaxs()`,
  `rowMaxs()`, `mean()`, `var()`, and `sd()`.
* Added support for statistical transformations including `scale()`,
  `sweep()`, and `cor()`.
* Added or extended methods for matrix decompositions and
  factorizations, including `svd()`, `prcomp()`, `qr()`, `chol()`,
  `solve()`, `eigen()`, and `pseudoinverse()`.
* Added diagonal operations and utilities, including `diag()`,
  `diag<-()`, `diag_op()`, and `diag_scale()`.
* Added split, reduce, and apply utilities for HDF5-backed workflows.
* Retained additional high-level `bd*` utilities for specialized
  workflows that do not map directly to standard R generics.

## HDF5 storage and resource management

* Added `list_datasets()` for inspecting datasets stored in HDF5 files.
* Added `hdf5_close_all()` for closing open HDF5 handles managed by the
  package.
* Improved HDF5 handle and pointer management, including safer behavior
  when overwriting datasets and when unloading or reloading the package.
* Added HDF5 compression handling and propagation of compression
  settings across output datasets.
* Added support for HDF5 file space management so that free space
  released by deleted datasets can be reused by subsequent writes in
  files created by the package.

## Performance and numerical methods

* Added support for block-wise and parallel execution settings through
  `hdf5matrix_options()`.
* Added TSQR support for tall-skinny QR decompositions, with automatic
  method selection when appropriate.
* Improved SVD and PCA handling for HDF5-backed matrices, including
  support for truncated outputs and block-wise computation.
* Improved QR, Cholesky, inverse, eigen decomposition, and
  pseudoinverse workflows for HDF5-backed matrices.
* Removed redundant OpenMP critical sections where HDF5 locking already
  provides the required synchronization.
* Improved thread handling for CRAN-compatible execution.

## C++ infrastructure

* Extended the C++ backend with classes and routines for managing HDF5
  files, groups, and datasets.
* Exposed reusable block-wise computational infrastructure through the
  package headers for developers implementing new Rcpp-based methods.
* Improved integration between the R interface and the C++ backend for
  matrix operations, decompositions, and HDF5 resource management.
* Added or improved C++ implementations used by the R/S3 interface for
  matrix algebra, decompositions, transformations, and HDF5-backed
  statistical operations.

## Documentation and examples

* Reworked the main package vignette around the `HDF5Matrix` interface
  and standard R methods.
* Improved the package-level help page opened by
  `help("BigDataStatMeth")`.
* Cleaned examples to use temporary files and avoid writing to the
  user's working directory.
* Removed obsolete vignette material and outdated installation examples.
* Moved example data to the package `extdata` directory.
* Improved documentation consistency for HDF5 terminology and package
  usage.

## Bug fixes and reliability

* Fixed error propagation in matrix multiplication and crossproduct
  paths.
* Fixed signed integer overflow in Cholesky block-size computation:
  `minimumBlockSize` promoted from `int` to `long` in
  `Cholesky_decomposition_intermediate_hdf5` and
  `Inverse_of_Cholesky_decomposition_intermediate_hdf5` to prevent
  overflow before promotion to `double` for `sqrt()`. Detected by
  gcc-ASAN on R-hub.
* Fixed index bug in block SVD nzeros threshold loop
  (`matrixSvdBlock.hpp`): wrong singular value was evaluated in the
  rank-truncation check, causing incorrect nzeros counts.
* Fixed thread-safety crash in multi-level block SVD (`H5SL_insert`
  from `H5I_register`) when the number of hierarchical levels `q >= 2`:
  `Next_level_SvdBlock_decomposition_hdf5` is now sequential; the outer
  block loop already provides sufficient parallelism.
* Fixed `nev` parameter not being applied to per-block truncation in
  `First_level_SvdBlock_decomposition_hdf5`: the parameter was declared
  in the function signature but never used, causing unnecessarily large
  intermediate matrices and a final SVD over more components than
  requested.
* Fixed several HDF5 pointer, handle, and finalizer edge cases.
* Fixed portability issues related to platform-specific memory queries.
* Fixed several edge cases in diagonal operations, matrix-vector
  operations, sparse multiplication internals, and SVD/PCA parameter
  handling.
* Improved cleanup of HDF5 handles during package unloading.
* Improved behavior when repeatedly creating, overwriting, closing, and
  reopening HDF5-backed datasets.


# BigDataStatMeth 1.0.3

* Minor documentation and example updates.
* Version previously available on CRAN before the 2.0.0 redesign.


# BigDataStatMeth 1.0.2

## Major changes
* Reduced example execution time for CRAN compliance
* Fixed bdblockmult_hdf5 example (< 5 seconds now)

## Bug fixes
* Fixed Makevars configuration issue causing compilation warnings


# BigDataStatMeth 1.0.1

## CRAN Resubmission Fixes
* **Thread management**: Implemented respect for OMP_THREAD_LIMIT environment variable
* **Default threads**: Limited to maximum 2 threads on CRAN servers to prevent excessive CPU usage
* **Documentation**: Replaced all Unicode characters with proper LaTeX macros (\eqn{}, \deqn{})
* **HTML validation**: Corrected invalid HTML tags in all documentation files
* **ATLAS compatibility**: Tested with standard BLAS/LAPACK configurations

## Bug Fixes
* Fixed thread safety in parallel HDF5 operations
* Improved file locking mechanisms for concurrent access
* Corrected dimension handling in transposed operations
* Enhanced error messages for invalid inputs

## Documentation Improvements
* Updated all function examples with proper mathematical notation
* Added comprehensive CRAN submission notes
* Improved vignette with clearer explanations
* Fixed formatting issues in Rd files


# BigDataStatMeth 1.0.0

## Major changes
* Complete rewrite of the package after the archived version 0.99.32.
* New block-wise computing framework for large matrices stored in HDF5.
* New C++ backend integrated with R through Rcpp and Rhdf5lib.
* API redesigned; not backwards compatible with previous versions.
* Package now focuses on providing scalable building blocks to develop
  new statistical methods for large datasets.

## New features
* Block-wise matrix multiplication for in-memory and HDF5-backed matrices.
* Block-wise SVD and PCA implementations.
* Block-wise QR decomposition.
* Block-wise crossproduct and matrix operations using HDF5 storage.
* Canonical Correlation Analysis (CCA) for HDF5 matrices via
  `bdCCA_hdf5_rcpp()`.
* Improved HDF5 import utilities, including support for large text files.
* Support for parallel computation (OpenMP) when available.
* Low memory footprint for large matrices exceeding system RAM.

## Improvements
* Substantial performance improvements in matrix operations on HDF5 data.
* More robust dimension handling and HDF5 metadata management.
* Unified interface for in-memory and on-disk data.
* Better error handling and validation on input dimensions and block sizes.

## Removed
* Old implementations from the 0.99.x series have been removed.
* Deprecated functions and APIs have been replaced by the new block-wise
  framework.

