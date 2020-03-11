#ifndef pca
#define pca

  #include <RcppEigen.h>
  #include "matrix_utilities.h"
  #include "svdDecomposition.h"
  #include "optimizedproduct.h"
  #include "parallelBlockMult.h"

  Rcpp::RObject bdPCA(const Rcpp::RObject & x, int k);



#endif