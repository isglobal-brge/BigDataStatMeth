#include "include/ReadDelayedData.h"


// [[Rcpp::depends(RcppEigen)]]


Eigen::MatrixXd read_DelayedArray_int( Rcpp::RObject A )
{

  auto dmat = beachmat::create_integer_matrix(A);
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  Eigen::MatrixXd eX(nrows, ncols);
  
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
  // size_t ncols = 0, nrows=0;
  
  auto dmat = beachmat::create_numeric_matrix(A);

  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  Eigen::MatrixXd eX(nrows, ncols);
  // eX.resize(nrows, ncols);
  
  
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







Rcpp::NumericMatrix read_DelayedArray_int_r( Rcpp::RObject A )
{

  auto dmat = beachmat::create_integer_matrix(A);
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  Rcpp::NumericMatrix eX(nrows, ncols);

  if ( nrows < ncols )
  {
    Rcpp::IntegerVector output(dmat->get_nrow());
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      dmat->get_col(ncol, output.begin());
      eX.column(ncol) = output;
    }
  }else {
    Rcpp::IntegerVector output(dmat->get_ncol());
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      eX.row(nrow) = output;
    }
  }
  return eX;
}


Rcpp::NumericMatrix read_DelayedArray_real_r( Rcpp::RObject A )
{
  auto dmat = beachmat::create_numeric_matrix(A);
  
  size_t ncols = dmat->get_ncol();
  size_t nrows = dmat->get_nrow();
  
  Rcpp::NumericMatrix eX(nrows, ncols);
  
  if ( ncols < nrows )
  {
    Rcpp::NumericVector output(nrows);
    for (size_t ncol=0; ncol<ncols; ++ncol) {
      dmat->get_col(ncol, output.begin());
      eX.column(ncol) = output;
    }
  }else {
    Rcpp::NumericVector output(ncols);
    for (size_t nrow=0; nrow<nrows; ++nrow) {
      dmat->get_row(nrow, output.begin());
      eX.row(nrow) = output;
    }
  }
  return eX;
}




Rcpp::NumericMatrix read_DelayedArray_rcpp( Rcpp::RObject A  )
{
  auto dmtypex = beachmat::find_sexp_type(A);
  // size_t ncols = 0, nrows=0;
  
  if ( dmtypex == INTSXP )
  {
    return (read_DelayedArray_int_r(A));
  } else if (dmtypex==REALSXP)
  {
    return (read_DelayedArray_real_r(A));
  }else
  {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}



