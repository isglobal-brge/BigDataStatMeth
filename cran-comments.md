## Resubmission

This is a complete rewrite of BigDataStatMeth that was archived in 2022.

### Changes addressing previous CRAN concerns:

* **Thread management**: Now respects OMP_THREAD_LIMIT and defaults to maximum 2 threads on CRAN servers
* **Unicode characters**: All replaced with proper LaTeX macros (\eqn{}, \deqn{})
* **HTML validation**: Corrected invalid tags in documentation files

### Note on previous ATLAS issues

The package was previously archived partly due to ATLAS/BLAS issues. 
The current version has been tested with standard BLAS/LAPACK and 
we request verification on CRAN's ATLAS systems during review.

### Test environments

* local: macOS aarch64-apple-darwin20, R 4.5.1
* win-builder: R-devel and R-release  
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), windows-latest (release)

### R CMD check results

**All platforms**: 0 errors ✓ | 0 warnings ✓ | 1 note x

**NOTE**: 
* "Package was archived on CRAN" - Historical note from 2022, all previous issues have been resolved

**Additional notes**:
* GitHub Actions macOS runner shows OpenSSL linking issues (runner configuration), but package builds successfully on local macOS and CRAN's macOS servers have proper OpenSSL configuration
* RcppEigen compilation warnings on some platforms are from external library code, not package code