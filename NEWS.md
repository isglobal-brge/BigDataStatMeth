# BigDataStatMeth 2.0.26

## Bug fixes

* Fixed `multiply_sparse()`: three bugs in `multiplicationSparse.hpp` —
  incorrect output dataset dimensions (`createDataset(M,N)` → `createDataset(N,M)`),
  accumulation reading uninitialised data on first k-block (`kk==0` guard),
  and incorrect write path via `Rcpp::wrap` (replaced with `std::vector<double>`
  overload, matching `multiplication.hpp` pattern).
* Fixed `diag.default`: removed `names` parameter from `base::diag()` call
  that caused incorrect dispatch returning 1×1 instead of N×N identity.
* Fixed `%*%.default`: added S3 default method for plain R matrices when
  package is loaded; `@rawNamespace` directives prevent malformed NAMESPACE.
* Fixed `~hdf5Dataset()` destructor: `H5Iis_valid(pdataset->getId())` guard
  before `close()` prevents crash when handles are invalidated externally.
* Fixed `rcpp_hdf5dataset_close()`: `R_ExternalPtrAddr(ptr)==nullptr` guard
  prevents double-close after `R_ClearExternalPtr()`.
* Fixed `svd.HDF5Matrix`: `center`/`scale` parameters were silently ignored
  (passed to `...`) — now correctly forwarded to the R6 `$svd()` method.
* Fixed `rcpp_hdf5dataset_diag_scale()`: coordinate mismatch for 1×N vector
  datasets — now uses HDF5-native coords for read/write and R coords for
  `createDataset`.
* Fixed `getAvailableMemoryMB()`: `memory.size()` call wrapped in
  `#ifdef _WIN32` to eliminate 35 cross-platform warnings.

## New features

* `length.HDF5Matrix`: new S3 method returning `prod(dim(x))`, consistent
  with `base::length()` for R matrices.
* `list_datasets()`: new S3 function wrapping `bdgetDatasetsList_hdf5()` —
  lists datasets in an HDF5 group from an `HDF5Matrix` object or file path.
* `rcpp_hdf5_close_all_registry()`: extended to call
  `BigDataStatMeth::closeAllHDF5Handles()` (new function in
  `hdf5CheckClose.hpp`) which closes all HDF5 C-library handles not tracked
  by the R6 registry, equivalent to `rhdf5::h5closeAll()`.
* `.onUnload()`: now calls `rcpp_hdf5_close_all_registry()` before unloading
  the shared library, preventing finalizer crashes on package reload.
* `hdf5_close_all()`: now calls `rcpp_hdf5_close_all_registry()` at the end,
  closing all HDF5 handles including those not tracked as `HDF5Matrix` objects.

## CRAN compliance

* Eliminated 35 `memory.size() is Windows-specific` warnings via
  `#ifdef _WIN32` guard.
* Fixed `mem_tcrossprod.cpp`: two empty `catch{}` blocks now re-throw via
  `Rf_error()`, propagating `non-conformable` errors correctly to R.


# BigDataStatMeth 2.0.25

## Major changes

* New hybrid R6+S3 architecture: all operations are now accessible via a
  standard S3 interface. Users interact exclusively through S3 generics
  (`dim`, `[`, `%*%`, `crossprod`, `svd`, `qr`, `cor`, `scale`, etc.);
  the R6 layer is internal.
* Phase 17: implemented S3/R6 wrappers for `eigen()`, `pseudoinverse()`,
  `diag()`/`diag<-()`, `diag_op()`, `diag_scale()`, `sweep()`,
  matrix-vector product, `multiply_sparse()`, `split_dataset()`,
  `reduce()`, and `apply_function()`.
* Phase 16: automatic compression inheritance — all output datasets inherit
  the compression level of their input when not explicitly specified.
  Global override via `hdf5matrix_options(compression = 0..9)`.
* Phase 18: removed 46 redundant `#pragma omp critical(accessFile)`
  directives (double serialisation with internal `HDF5ThreadSafety::LockGuard`).
  One intentional critical section retained in TSQR (`hdf5_tsqr_read`).
* HDF5 pointer safety: live-pointer registry changed from
  `unordered_set<void*>` to `unordered_map<void*, SEXP>`; `claim_ptr()`
  invalidates XPtr handles via `R_ClearExternalPtr()`; new
  `rcpp_hdf5_close_at_paths()` closes stale handles before overwrite.
  Sentinel-based R6 finalizer ensures deterministic handle closure.
* TSQR (Tall-Skinny QR): new parallel block QR algorithm for tall-skinny
  matrices (`m/n > 5` and `m > 1000`). Selectable via `method = "tsqr"`
  or automatically dispatched with `method = "auto"`.

## Bug fixes

* Fixed double-serialisation in parallel HDF5 I/O by removing redundant
  `omp_critical` sections (Phases 18).
* Fixed `ARM64` heap corruption from Eigen expression templates in
  eigendecomposition (replaced `.reverse()` with explicit index loops).
* Fixed HDF5 file corruption from `writeDataset(Rcpp::wrap(VectorXd))`
  using wrong dataspace rank.
* Fixed `ncomponents` truncation in `prcomp()` for large matrices.
* Fixed `SQUARE_SMALL` block sizing in `multiplication.hpp` (removed
  erroneous `/2` divisor).
* Fixed `diag<-` S3 dispatch via `@rawNamespace` directives.

## CRAN compliance

* Thread compliance: all thread-count selection now respects
  `OMP_THREAD_LIMIT` via `getDTthreads(INT_MAX, false)` (Fixes A–E).
* Removed `exportPattern` from `NAMESPACE`.
* All new examples use `tempfile()`/`tempdir()`.
* Corrected HDF5 capitalisation across all 175 man pages.
* Removed `library(SeqArray)` from vignette.
* Removed dead code (`S3_subset_old.R`, `rcpp_hdf5_count_open_datasets`).


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

