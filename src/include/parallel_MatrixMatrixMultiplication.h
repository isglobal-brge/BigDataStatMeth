#ifndef prallel_MatrixMatrixMultiplication
#define prallel_MatrixMatrixMultiplication

  #include <RcppEigen.h>
  #include <RcppParallel.h>
  #include "matrix_utilities.h"

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ColumnVector;
  
  Rcpp::NumericMatrix rcpp_parallel_XYProd(Rcpp::NumericMatrix matx, Rcpp::NumericMatrix maty);
  Eigen::MatrixXd rcpp_parallel_XYProd_eigen(Eigen::MatrixXd matX, Eigen::MatrixXd matY);
  Rcpp::NumericMatrix rcpp_parallel_XYProdBlock(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY);
  Rcpp::NumericMatrix rcpp_parallel_XtYProd(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY);
  Rcpp::NumericMatrix rcpp_parallel_XYtProd(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY);
  


#endif
