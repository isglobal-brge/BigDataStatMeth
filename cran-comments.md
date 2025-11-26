## Test environments

* Local macOS arm64 (M1): R 4.5.x
* Local Linux x86_64: R 4.4.x
* Local Windows 10/11: R 4.4.x

* r-hub:
  - Linux R-devel (clang-ASan, clang-UBSan)
  - Linux R-devel with valgrind (clang toolchain)
  - Linux R-release
  - Windows
  - macOS-arm64

* rchk:
  - r-hub/rchk Docker image

## R CMD check results

0 errors | 1 warning | 0 notes

The single WARNING is produced during compilation of external headers
(RcppEigen and BH/Boost) on recent compilers, not from the package’s own
C++ sources.

## Resubmission

This is a resubmission after the archived version 0.99.32.  
The package has been largely rewritten, with a new HDF5-based backend and
a revised API (not backwards compatible).

## Additional checks

I have reviewed the clang-ASan, clang-UBSan, and valgrind logs from
r-hub, as well as the C++ compilation logs on recent toolchains.  
All detected diagnostics originate in external libraries
(HDF5/Rhdf5lib/RcppEigen/BH) and not in the BigDataStatMeth source code.
No memory errors, invalid reads/writes, or undefined behaviour were
reported for the package’s own C++ code.

On some platforms (e.g., macOS arm64 with recent compilers), R CMD check
reports a single WARNING during C++ compilation. These messages
(class-memaccess, Boost concept checks, deprecated std::auto_ptr, etc.)
come from RcppEigen and BH/Boost headers. The C++ code in
BigDataStatMeth itself compiles cleanly.

The rchk report shows warnings mainly related to template instantiation
within Rcpp dependencies. The package’s code has been reviewed for
correct PROTECT/UNPROTECT usage, and OpenMP support is optional and
properly detected at compile time.
