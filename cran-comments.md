## Resubmission

This is a complete rewrite of BigDataStatMeth that was archived in 2022.

### Changes addressing previous CRAN concerns:

* **Thread management**: Now respects OMP_THREAD_LIMIT and defaults to maximum 2 threads on CRAN servers
* **Unicode characters**: All replaced with proper LaTeX macros (\eqn{}, \deqn{})
* **HTML validation**: Corrected invalid tags in documentation files

### Test environments

* local: macOS aarch64-apple-darwin20, R 4.5.1
* win-builder: R-devel and R-release (results pending)

### R CMD check results

0 errors ✓ | 0 warnings ✓ | 1 note x

**NOTE**: 
* RcppEigen compilation warnings (external library code, not package code)

### Additional notes

Package size is 11.6Mb due to HDF5 C++ library dependencies and comprehensive 
documentation. This is necessary for the package's core functionality.