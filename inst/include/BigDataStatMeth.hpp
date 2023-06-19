#ifndef BigDataStatMeth_BigDataStatMeth_H
#define BigDataStatMeth_BigDataStatMeth_H

// Include the Rcpp Header
#include <RcppEigen.h>

// Common headers
#include "H5Cpp.h" //1.8.18
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sys/stat.h>
#include <string>


#include "Utilities/pkg_omp.h" // first for clang-13-omp, #5122
#include <R.h>
#include <Rversion.h>
#if !defined(R_VERSION) || R_VERSION < R_Version(3, 5, 0)  // R-exts$6.14
#  define ALTREP(x) 0     // #2866
#  define USE_RINTERNALS  // #3301
#  define DATAPTR_RO(x) ((const void *)DATAPTR(x))
#endif
#include <Rinternals.h>
#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped STRING_ELT and VECTOR_ELT
#include <stdint.h>    // for uint64_t rather than unsigned long long
#include <stdbool.h>
#ifdef WIN32  // positional specifiers (%n$) used in translations; #4402
#  define snprintf dt_win_snprintf  // see our snprintf.c; tried and failed to link to _sprintf_p on Windows
#endif
#ifdef sprintf
#undef sprintf
#endif

// #include <signal.h> // the debugging machinery + breakpoint aidee
// raise(SIGINT);
#define IS_TRUE_OR_FALSE(x) (TYPEOF(x)==LGLSXP && LENGTH(x)==1 && LOGICAL(x)[0]!=NA_LOGICAL)

// for use with bit64::integer64
#define NA_INTEGER64  INT64_MIN
#define MAX_INTEGER64 INT64_MAX

// openme-utils.cpp functions
void initDTthreads();
int getDTthreads(const int64_t n, const bool throttle);
static const char *mygetenv(const char *name, const char *unset);
SEXP getDTthreads_R(SEXP verbose);
void when_fork();
void after_fork();
void avoid_openmp_hang_within_fork();



// Load class headers from BigDataStatMeth
#include "hdf5Utilities/hdf5Files.hpp"
#include "hdf5Utilities/hdf5Groups.hpp"
#include "hdf5Utilities/hdf5Datasets.hpp"



// Load function definition from BigDataStatMeth



#endif