#ifndef parallel_CrossProd
#define parallel_CrossProd

  #include <RcppEigen.h>
  #include <cmath>
  #include <algorithm>
  #include <RcppParallel.h>
  #include "matrix_utilities.h"


  Rcpp::NumericMatrix rcpp_parallel_tCrossProd(Rcpp::NumericMatrix mat);
  Rcpp::NumericMatrix rcpp_parallel_tCrossProdblock(Rcpp::NumericMatrix mat);
  Rcpp::NumericMatrix rcpp_parallel_CrossProd(Rcpp::NumericMatrix mat);
  Rcpp::NumericMatrix rcpp_parallel_CrossProdblock(Rcpp::NumericMatrix mat);
  Rcpp::RObject partCrossProd(Rcpp::RObject X);
  Rcpp::RObject parCrossProd(Rcpp::RObject X);

#endif

