#include <RcppEigen.h>
#include "beachmat/numeric_matrix.h"   // To access numeric matrix
#include "beachmat/integer_matrix.h"   // To access integer matrix
#include "beachmat/character_matrix.h" // To access character matrix

using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::RObject BDtrans_numeric (const Rcpp::RObject & x)
{
  
  auto dmtype = beachmat::find_sexp_type(x);
  
  if ( dmtype == INTSXP ) {
    auto dmat = beachmat::create_integer_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    size_t ccol;//, crow;
    Eigen::MatrixXi proves(nrows,ncols);
    
    Rcpp::IntegerVector output(dmat->get_nrow());
    
    for (ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, output.begin());
      proves.col(ccol) = Rcpp::as<Eigen::VectorXi >(output);
    }
    
    proves.transposeInPlace() ;
    
    return wrap(proves);
    
  } else if (dmtype==REALSXP) {
    auto dmat = beachmat::create_numeric_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    size_t ccol;//, crow;
    Eigen::MatrixXd proves(nrows,ncols);
    
    Rcpp::NumericVector output(dmat->get_nrow());
    
    for (ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, output.begin());
      proves.col(ccol) = Rcpp::as<Eigen::VectorXd >(output);
    }
    
    proves.transposeInPlace() ;
    
    return wrap(proves);
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}


// [[Rcpp::export]]
Rcpp::RObject BDtrans_hdf5(const Rcpp::RObject & x)
{
  
  auto dmtype = beachmat::find_sexp_type(x);
  
  
  if ( dmtype == INTSXP ) {
    auto dmat = beachmat::create_integer_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    beachmat::output_param oparam(dmat->get_matrix_type(), FALSE, TRUE);  
    auto out_dmat = beachmat::create_integer_output(ncols, nrows, oparam);
    
    Rcpp::IntegerVector v_in(dmat->get_nrow());
    Rcpp::IntegerVector v_out(out_dmat->get_nrow());
    
    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, v_in.begin());
      out_dmat ->set_row(ccol,v_in.begin());
    }
    
    return (out_dmat->yield());
    
  } else if (dmtype==REALSXP) {
    
    auto dmat = beachmat::create_numeric_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    beachmat::output_param oparam(dmat->get_matrix_type(), FALSE, TRUE);  
    auto out_dmat = beachmat::create_numeric_output(ncols, nrows, oparam);
    
    Rcpp::NumericVector v_in(dmat->get_nrow());
    Rcpp::NumericVector v_out(out_dmat->get_nrow());
    
    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, v_in.begin());
      out_dmat ->set_row(ccol,v_in.begin());
    }
    
    return (out_dmat->yield());
    
  } else if (dmtype==STRSXP) {
    
    auto dmat = beachmat::create_character_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    beachmat::output_param oparam(dmat->get_matrix_type(), FALSE, TRUE);  
    auto out_dmat = beachmat::create_character_output(ncols, nrows, oparam);
    
    Rcpp::CharacterVector v_in(dmat->get_nrow());
    Rcpp::CharacterVector v_out(out_dmat->get_nrow());
    
    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, v_in.begin());
      out_dmat ->set_row(ccol,v_in.begin());
    }
    
    return (out_dmat->yield());
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

library(Rcpp)
library(DelayedArray)
library("microbenchmark")
  
  blah <- matrix(runif(10000), ncol=200, nrow=50)
  X <- DelayedArray(blah)

  results <- microbenchmark(transR = t(X),
                            transCpp = BDtrans_numeric(X),
                            transCpp = BDtrans_hdf5(X),
                            times = 50L)

*/
