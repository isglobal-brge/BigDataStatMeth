#ifndef matrixvector_mult
#define matrixvector_mult


#include <RcppEigen.h>
#include "ReadDelayedData.h"


  Rcpp::RObject bdwXw(Rcpp::RObject X, Rcpp::RObject W, std::string op);

#endif