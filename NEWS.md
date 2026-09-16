# BigDataStatMeth 2.0.6

## Bug fixes

- Windows: reopening an HDF5 file shortly after closing it no longer fails
  with "HDF5 file is in use by another process"; file locking is no longer
  forced on Windows. Linux and macOS are unchanged.
- Files are opened with a strong close degree, so an object that outlives
  its file handle can no longer keep the file busy and block a later reopen
  on Windows.
- Fixed the HDF5 handle release order in `bdCreate_hdf5_matrix()` error paths.
- `dim()` on a closed `HDF5Matrix` returns `NULL` and `length()` returns 0
  instead of signalling an error.
- `hdf5_close_all()` is silent when there is nothing to close.

# BigDataStatMeth 2.0.5

## Bug fixes

- Restored the OpenMP serialisation of HDF5 I/O: `impute_snps()` could crash,
  deadlock or silently corrupt its output on matrices with over 4000 columns.
- `hdf5_apply()` with `"CrossProd"` / `"tCrossProd"` returned an all-zero
  result for inputs of 1,048,576 elements or more.
- `svd()` / `prcomp()` with `scale = TRUE` now stop on a constant column
  instead of silently returning all-zero singular values. `scale()` is
  unchanged (`NaN`, as in base R).
- The `threads` argument of `svd()` / `prcomp()` is now honoured; it had no
  effect on the internal matrix products.
- `impute_snps()` and the omics filters now propagate dimnames, and the
  filters report a clear error when they remove every row or column.
- Cross-file `crossprod()`, `%*%` and `tcrossprod()` between `HDF5Matrix`
  objects stored in different files no longer fail.
- Windows memory detection no longer relies on `memory.size()` (defunct since
  R 4.2.0), which inflated the block-size and preload budget.
- Row and column names are no longer truncated to 19 characters: ids up to 63
  bytes round-trip exactly, and longer ids raise an error.

## Portability

- `checkHDF5File()` uses the portable HDF5 C API instead of the deprecated
  `H5::H5File::isHdf5()`, fixing a downstream link failure with
  Rhdf5lib >= 1.32.

## Documentation

- `filter_low_coverage()`: `pcent` counts entries equal to `3`, not `NA`.
- `filter_maf()`: the documented threshold direction was inverted (variants
  with MAF at or below the threshold are removed).

## New features

- `svd()` / `prcomp()` now record whether the decomposition was exact or
  approximate (new attributes plus a `message()`), and the new
  `svd_auto_threshold()` returns the size boundary. See `?svd.HDF5Matrix`.
- `svd()` / `prcomp()` warn when the requested `threads` cannot be honoured.
- `hdf5_remove()` and `HDF5Matrix$remove()` delete a dataset from a file.

# BigDataStatMeth 2.0.4

## Bug fixes

- `bdCorr_matrix()`: fixed a C++ exception (`Expecting a single value:
  [extent=0]`) that could abort the R session under recent Rcpp versions.
  The returned list evaluated an optional `Nullable<bool>` argument in a
  boolean context, which now resolves to `as<bool>(NULL)`; the unwrapped
  logical values are used instead. This resolves the CRAN r-devel example
  check failure.

## Error reporting / CRAN compliance

- Removed the literal `"c++ exception"` prefix from all error messages so
  failures surface as ordinary R errors rather than text that reads like a
  leaked C++ exception (per CRAN request). Error-handling control flow is
  unchanged.

## Documentation

- `cbind()` / `rbind()`: documentation no longer states that plain R matrices
  are accepted and auto-converted (they are not). Convert inputs first with
  `hdf5_matrix()` / `hdf5_create_matrix()`.
- Removed `RSpectra` from the README system-requirements list (it is not a
  package dependency; the bundled C++ Spectra library is used internally).
- Consolidated the duplicate `Ops.HDF5Matrix` definition into a single
  documented method and documented the 1xn / nx1 broadcast behavior.

## Internal

- Removed stale backup and superseded source files from the tree.

# BigDataStatMeth 2.0.3

### Bug fixes

- `qr(..., method = "tsqr")`: fixed incorrect results (loss of orthogonality
  in Q) when the last row-block had fewer rows than columns (`m mod block_size
  < n`). The last block is now merged with the previous one, guaranteeing all
  blocks satisfy the rank condition required by the local Householder QR step.
  Affected configurations: any (m, n, threads) where `m mod max(2n, ceil(m /
  max(4, threads * 4))) < n`. Discovered via a systematic regression test
  covering 202 (m, n, threads, block_size) combinations.
  
- `pseudoinverse()` now correctly zeroes singular values at or below the
  tolerance (1e-9), fixing incorrect results for rank-deficient or
  near-singular matrices.
- Fixed silent failure on matrix conversion errors in `bdCrossprod()`,
  `bdtCrossprod()`, and `bd_wproduct()`; these now raise informative errors.
- Removed unreachable code in the correlation HDF5 backend (no functional
  effect).

## New features

* Added `pseudoinverse.matrix()` for consistent in-memory support alongside
  the existing `HDF5Matrix` method.

## Documentation

- Corrected `pseudoinverse()` and `prcomp(rank. =)` documentation to match
  actual behavior.
- Updated package description and citation metadata.

# BigDataStatMeth 2.0.2

## Performance

- PATH 1 preload strategy added to `crossprod()`, `tcrossprod()`,
  aggregation functions, `scale()`, and `solve()`: when input fits
  within 20% of available RAM, a single HDF5 read + BLAS call + single
  write replaces the block loop.
- `scale()` PATH 2 (block-wise streaming) is now parallelized with OpenMP.
- `hdf5matrix_options(paral = TRUE)` now forces PATH 2 (OMP streaming)
  in `%*%`, `crossprod()`, `tcrossprod()`, aggregations, and `scale()`,
  giving explicit thread control regardless of matrix size.
- `getAvailableMemoryMB()` now uses platform-native detection (macOS:
  `host_statistics64`; Linux: `/proc/meminfo`), replacing the hardcoded
  4 GB fallback.
- TSQR: Step 3 Q-assembly parallelized; block size decoupled from thread
  count (`target_blocks = max(4, nthreads * 4)`).

## Bug fixes

- `rbind()` and `cbind()` with 3 or more `HDF5Matrix` arguments no
  longer produce dataset names exceeding 1400 characters. Auto-generated
  names are now `rbind_N_<uid>` with a random 8-character identifier,
  also preventing "dataset already exists" collisions. Two-argument
  calls are unchanged (`A_rbind_B`).
- Fixed `paral = FALSE` being silently ignored in `multiplication()`
  PATH 2/3.

## Correctness

- Restored `#pragma omp critical(accessFile)` in `crossprod()`,
  `tcrossprod()`, and aggregation PATH 2 (race condition fix).
- Fixed SIGABRT in `Cholesky_decomposition_intermediate_hdf5` on
  macOS ARM64 (removed `throw` from OpenMP parallel region).

## Documentation

- `scale()`: added `paral`, `threads` parameters; PATH 1/PATH 2
  `@details`; `@seealso hdf5matrix_options()`.
- `split()`, `hdf5_apply()`, `hdf5_reduce()`: clarified calling
  conventions and two-level access pattern. Fixed legacy example.
- Added vignette section "Dataset operations: split, reduce, and apply".

# BigDataStatMeth 2.0.1

* Fixed vignette build failure on Debian: replaced external CSV file
  with reproducible synthetic data in the `hdf5_import` example.

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

