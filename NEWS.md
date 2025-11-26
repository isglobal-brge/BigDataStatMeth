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

