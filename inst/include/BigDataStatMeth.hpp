#ifndef BigDataStatMeth_BigDataStatMeth_H
#define BigDataStatMeth_BigDataStatMeth_H

#ifndef _COMMON_HEADERS_H
#define _COMMON_HEADERS_H
    // Include the Rcpp Header
    #include <RcppEigen.h>
    
    // Common headers
    #include "H5Cpp.h" //1.8.18
    #include <iostream>
    #include <boost/algorithm/string.hpp>
    #include <fstream>
    #include <sys/stat.h>
    #include <string>
    #include <regex>
    
    // Common headers - multi-threading
    #include <thread>
#endif


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

// Constants

#define MAX_NAME 1024
// #define MAXSVDBLOCK 1500 //. ORIGINAL - CORRECTE.//

const int RANK1 = 1;
const int RANK2 = 2;
const int RANK3 = 3;
const int DIM1 = 1;
const int DIM2 = 2;
const int DIM3 = 3;

const int MAXSTRING = 20;
const hsize_t MAXSTRBLOCK = 1 << 5;

const hsize_t MAXELEMSINBLOCK = ((2 << 29) - 1);
// const hsize_t MAXELEMSINBLOCK = ((2 << 5) - 1);
const hsize_t MAXBLOCKSIZE = std::floor(std::sqrt(MAXELEMSINBLOCK));
const hsize_t MAXMULTBLOCKSIZE = 1 << 13;

// const int maxElemBlock = 3000000;
//.Only test.// const int maxElemBlock = 30;
const int EXEC_OK = 0;
const int EXEC_ERROR = 1;
const int EXEC_WARNING = 2;

// openme-utils.cpp functions
void initDTthreads();
int getDTthreads(const int64_t n, const bool throttle);
static const char *mygetenv(const char *name, const char *unset);
SEXP getDTthreads_R(SEXP verbose);


// Load class headers from BigDataStatMeth


#ifndef _BIGDATASTATMETH_BASIC_HEADERS_HPP
#define _BIGDATASTATMETH_BASIC_HEADERS_HPP

    // Classes
    #include "hdf5Utilities/hdf5Files.hpp"
    #include "hdf5Utilities/hdf5Groups.hpp"
    #include "hdf5Utilities/hdf5Datasets.hpp"
    #include "hdf5Utilities/hdf5DatasetsInternal.hpp"
    #include "hdf5Utilities/hdf5Dims.hpp"
    
    #include "hdf5Utilities/hdf5CheckClose.hpp"

#endif

// Load function definition from BigDataStatMeth
#ifndef _BIGDATASTATMETH_FUNCTION_HEADERS_HPP
#define _BIGDATASTATMETH_FUNCTION_HEADERS_HPP
    #include "hdf5Algebra.hpp"  // Algebra (hdf5 data files)
    #include "hdf5Utilities.hpp" // Utilities (hdf5 data files)
    #include "memAlgebra.hpp"  // Algebra (on-memort)
#endif




// // Algebra related functions
    // #include "hdf5Algebra/matrixNormalization.hpp"
    // #include "hdf5Algebra/matrixSdMean.hpp"
    // #include "hdf5Algebra/multiplication.hpp"
    // #include "hdf5Algebra/matrixPCA.hpp"
    // #include "hdf5Algebra/matrixEquationSolver.hpp"
    // #include "hdf5Algebra/matrixSvd.hpp"
    // #include "hdf5Algebra/matrixNormalization.hpp"
    // #include "hdf5Algebra/matrixSvdBlock.hpp"
    // #include "hdf5Algebra/multiplicationSparse.hpp"
    // #include "hdf5Algebra/matrixQR.hpp"
    // #include "hdf5Algebra/tcrossprod.hpp"
    // #include "hdf5Algebra/matrixInvCholesky.hpp"
    // #include "hdf5Algebra/matrixSum.hpp"
    // #include "hdf5Algebra/matrixDiagonal.hpp"
    // #include "hdf5Algebra/vectormatrix.hpp"
    // #include "hdf5Algebra/matrixSubstract.hpp"
    // #include "hdf5Algebra/matrixPseudoinverse.hpp"
    // #include "hdf5Algebra/matrixTriangular.hpp"
    // #include "hdf5Algebra/crossprod.hpp"

// Omics related functions
    // #include "hdf5Omics/hdf5RemoveMAF.hpp"
    // #include "hdf5Omics/hdf5OmicsUtils.hpp"


// Other hdf5 functions
    // #include "hdf5Utilities/hdf5ReduceDataset.hpp"
    // #include "hdf5Utilities/hdf5SplitDataset.hpp"
    // #include "hdf5Utilities/hdf5ImportFiles.hpp"
    // #include "hdf5Utilities/hdf5BindDatasets.hpp"
    // #include "hdf5Utilities/hdf5Methods.hpp"
    // #include "hdf5Utilities/hdf5RemoveElements.hpp"
    // #include "hdf5Utilities/hdf5RemoveLowData.hpp"
    // #include "hdf5Utilities/hdf5Dims.hpp"
    // #include "hdf5Utilities/hdf5SortDataset.hpp"
    // #include "hdf5Utilities/hdf5ImputeData.hpp"

// On-memory related functions
    // #include "memAlgebra/memOtherFunctions.hpp"
    // #include "memAlgebra/memOptimizedProducts.hpp"
    // #include "memAlgebra/memMultiplication.hpp"
    // #include "memAlgebra/memSubstract.hpp"
    // #include "memAlgebra/memSum.hpp"


#endif