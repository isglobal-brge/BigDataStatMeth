#ifndef parallel_MatrixVectorMultiplication
#define parallel_MatrixVectorMultiplication

  #include <RcppEigen.h>
  #include <RcppParallel.h>
  #include <cmath>

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ColumnVector;

  Rcpp::NumericMatrix rcpp_xwxt(Rcpp::NumericMatrix mat, Rcpp::NumericVector w);
  //..// Eigen::Map<Eigen::MatrixXd> rcpp_xwxt_eig(const Eigen::Map<Eigen::MatrixXd>* mat, const Eigen::Map<Eigen::VectorXd>* w, Eigen::Map<Eigen::VectorXd>* rmat);
  void rcpp_xwxt_eig(const Eigen::Map<Eigen::MatrixXd>* mat, const Eigen::Map<Eigen::VectorXd>* w, Eigen::Map<Eigen::MatrixXd>* rmat);
  
  Rcpp::NumericMatrix rcpp_xtwx(Rcpp::NumericMatrix mat, Rcpp::NumericVector w);
  Rcpp::NumericMatrix rcpp_parallel_Xy(Rcpp::NumericMatrix mat, Rcpp::NumericVector y);
  //..// Eigen::VectorXd rcpp_parallel_Xy_eigen( Eigen::MatrixXd mat, Eigen::VectorXd y);
  void rcpp_parallel_Xy_eigen( Eigen::Map<Eigen::MatrixXd>* mat, const Eigen::Map<Eigen::VectorXd>* y, Eigen::Map<Eigen::VectorXd>* rmat);
  Rcpp::RObject xwxt(Rcpp::RObject mat, Rcpp::RObject w);
  Rcpp::RObject xtwx(Rcpp::RObject mat, Rcpp::RObject w);

#endif
