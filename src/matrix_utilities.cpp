#include "include/matrix_utilities.h"
#include "include/parallel_CrossProd.h"



Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x)
{
  double xcol = x.ncol();
  double xrow = x.nrow();
  double dim = xrow * xcol;
  
  Rcpp::NumericVector fr( dim );
  
  for( double i=0; i< xrow; i++)
  {
    for( double j=0; j<xcol; j++)
    {
      fr[i*xcol+j] = x(i,j);
    }
  }
  
  return(fr);
}


Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x)
{
  double xcol = x.ncol();
  double xrow = x.nrow();
  double dim = xrow * xcol;
  
  Rcpp::NumericVector fc( dim );
  
  for( double i=0; i< xrow; i++)
  {
    for( double j=0; j<xcol; j++)
    {
      fc[i+j*xrow] = x(i,j);
    }
  }
  return(fc);
}


Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X )
{
  Eigen::RowVectorXd mean = X.colwise().mean();
  Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
  return (X.rowwise() - mean).array().rowwise() / std.array();
}



Rcpp::NumericMatrix RcppNormalize_Data_r ( Rcpp::NumericMatrix  x )
{
  Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd> (x);
  Eigen::RowVectorXd mean = X.colwise().mean();
  Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
  return Rcpp::wrap((X.rowwise() - mean).array().rowwise() / std.array());
}



// [[Rcpp::export]]
Rcpp::RObject Normalize_Data ( Rcpp::RObject & x )
{
  
  Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd>(x);
  Eigen::RowVectorXd mean = X.colwise().mean();
  Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
  X = (X.rowwise() - mean).array().rowwise() / std.array();
  return Rcpp::wrap(X);
}



