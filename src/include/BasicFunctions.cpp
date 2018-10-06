#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include "beachmat/numeric_matrix.h"  // To access numeric matrix
#include "beachmat/integer_matrix.h"  // To access numeric matrix
#include<cmath>


//' Normalize Data
//' 
//' This function returns normalized matrix ...
//'
//' @param RObject x
//' @export
// [[Rcpp::export]]
Rcpp::RObject Normalize_Data (const Rcpp::RObject & x )
{
  Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd>(x);
  Eigen::RowVectorXd mean = X.colwise().mean();
  Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
  X = (X.rowwise() - mean).array().rowwise() / std.array();
  return Rcpp::wrap(X);
}

//' Cros-Product Eigen
//' 
//' Return a matrix cross-product
//'
//' @param Eigen::MatrixXd x
//' @param bool tpc if true returns x %*% t(x) if tpc=false returns t(x) %*% x
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd CrossProduct_eig (const Eigen::MatrixXd x, bool tpc )
{
  if ( tpc == true) {
    return (x * x.transpose());
  }else {
    return (x.transpose() * x);
  }
}

//' Cros-Product
//' 
//' Return a matrix cross-product
//'
//' @param RObject x
//' @param bool tpc if true returns x %*% t(x) if tpc=false returns t(x) %*% x
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd CrossProduct (const Rcpp::RObject & x, bool tpc )
{
  Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd>(x);
  if ( tpc == true) {
    return (X * X.transpose());
  }else {
    return (X.transpose() * X);
  }
}


/*** R
timesTwo(42)
*/
