#include "include/ReadDelayedData.h"


// [[Rcpp::depends(RcppEigen)]]


Eigen::MatrixXd read_DelayedArray_int( Rcpp::RObject A )
{
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  
  auto dmat = beachmat::create_integer_matrix(A);
  ncols = dmat->get_ncol();
  nrows = dmat->get_nrow();
  
  eX.resize(nrows, ncols);
  
  if ( nrows < ncols )
  {
    Rcpp::IntegerVector output(dmat->get_nrow());
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      dmat->get_col(ncol, output.begin());
      eX.col(ncol) = Rcpp::as<Eigen::VectorXd >(output);
    }
  }else {
    Rcpp::IntegerVector output(dmat->get_ncol());
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      eX.row(nrow) = Rcpp::as<Eigen::VectorXd >(output);
    }
  }
  return eX;
}


Eigen::MatrixXd read_DelayedArray_real( Rcpp::RObject A )
{
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  
  auto dmat = beachmat::create_numeric_matrix(A);

  ncols = dmat->get_ncol();
  nrows = dmat->get_nrow();
  eX.resize(nrows, ncols);
  
  
  if ( ncols < nrows )
  {
    Rcpp::NumericVector output(nrows);
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      dmat->get_col(ncol, output.begin());
      eX.col(ncol) = Rcpp::as<Eigen::VectorXd >(output);
    }
  }else {
    Rcpp::NumericVector output(ncols);
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      eX.row(nrow) = Rcpp::as<Eigen::VectorXd >(output);
    }
  }
  return eX;
}


Eigen::MatrixXd read_DelayedArray( Rcpp::RObject A  )
{
  auto dmtypex = beachmat::find_sexp_type(A);
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  
  if ( dmtypex == INTSXP )
  {
    return (read_DelayedArray_int(A));
  } else if (dmtypex==REALSXP)
  {
    return (read_DelayedArray_real(A));
  }else
  {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}

bool isDelayedData()
{
  
}

