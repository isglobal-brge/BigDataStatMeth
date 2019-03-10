#ifndef parallel_CrossProd
#define parallel_CrossProd

  #include <RcppEigen.h>
  #include <cmath>
  #include <algorithm>
  #include <RcppParallel.h>
  #include "matrix_utilities.h"

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ColumnVector;
  typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> InversMat;
  typedef Eigen::Map<InversMat> MapInversMat;

  Rcpp::NumericMatrix rcpp_parallel_tCrossProd(const Rcpp::NumericMatrix& mat);
  Eigen::Map<Eigen::MatrixXd> rcpp_parallel_tCrossProd_eigen(const Eigen::Map<Eigen::MatrixXd>& mat);
  Rcpp::NumericMatrix rcpp_parallel_tCrossProdblock(Rcpp::NumericMatrix mat);
  Rcpp::NumericMatrix rcpp_parallel_CrossProd(Rcpp::NumericMatrix mat);
  Rcpp::NumericMatrix rcpp_parallel_CrossProdblock(Rcpp::NumericMatrix mat);
  Rcpp::RObject partCrossProd(Rcpp::RObject X);
  Rcpp::RObject parCrossProd(Rcpp::RObject X);

#endif

