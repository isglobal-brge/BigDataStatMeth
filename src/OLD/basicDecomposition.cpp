// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include "beachmat/numeric_matrix.h"  // To access numeric matrix
#include "beachmat/integer_matrix.h"  // To access numeric matrix
#include "include/BasicFunctions.h"
#include "include/spectra/SymEigsSolver.h"  // To access symmetric matrix


// FUNCIONS BEACHMAT
// [[Rcpp::export]]
arma::mat BDCrossProduct_benchmat_nodelayed ( Rcpp::NumericMatrix x, bool tpc )
{
  arma::mat mx = Rcpp::as<arma::mat >(x);
  if (tpc == true) {
    return (mx * mx.t()); // X'*X
  }else {
    return (mx.t() * mx); // X*X'
  }
  
}

// FI FUNCIÓ PROVA BENCHMAT


// [[Rcpp::export]]
Rcpp::RObject BDCrossProduct (const Rcpp::RObject & x, bool tpc )
{
  auto dmtype = beachmat::find_sexp_type(x);
  
  if ( dmtype == INTSXP ) {
    
    auto dmat = beachmat::create_integer_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    Rcpp::IntegerVector output(dmat->get_nrow());
    arma::mat amx(nrows, ncols);
    
    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, output.begin());
      amx.col(ccol) = Rcpp::as<arma::colvec >(output);
    }
    
    if ( tpc == true) {
      return Rcpp::wrap(amx * amx.t());
    }else {
      
      return Rcpp::wrap(amx.t() * amx);
    }
    
  } else if (dmtype==REALSXP) {
    
    auto dmat = beachmat::create_numeric_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    Rcpp::NumericVector output(dmat->get_nrow());
    arma::mat amx(nrows, ncols);
    
    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, output.begin());
      amx.col(ccol) = Rcpp::as<arma::colvec >(output);
    }
    
    if ( tpc == true ) {
      return Rcpp::wrap(amx * amx.t());
    }else {
      return Rcpp::wrap(amx.t() * amx);
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}



// [[Rcpp::export]]

Rcpp::RObject BDCInverse_Cholesky (const Rcpp::RObject & x )
{
  auto dmtype = beachmat::find_sexp_type(x);
  
  if ( dmtype == INTSXP ) {
    
    auto dmat = beachmat::create_integer_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    Rcpp::IntegerVector output(dmat->get_nrow());
    arma::mat mx(nrows, ncols);
    
    // bool tcprod = true;
    
    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, output.begin());
      mx.col(ccol) = Rcpp::as<arma::colvec >(output);
    }
    if(mx.is_hermitian()) 
    {
      arma::mat minvC(nrows, ncols), R;  
      
      chol( R, mx );
      minvC = mx.t()* inv_sympd((mx*mx.t()));
      
      return Rcpp::wrap(minvC); 
      
    }else {
      throw std::runtime_error("Matrix must be symmetric");
    }

  } else if (dmtype==REALSXP) {
    
    auto dmat = beachmat::create_numeric_matrix(x);
    const size_t ncols = dmat->get_ncol();
    const size_t nrows = dmat->get_nrow();
    
    Rcpp::NumericVector output(dmat->get_nrow());
    arma::mat mx(nrows, ncols);

    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, output.begin());
      mx.col(ccol) = Rcpp::as<arma::colvec >(output);
    }

    if(mx.is_hermitian()) 
    {
      arma::mat minvC(nrows, ncols), R(nrows, ncols);  
  
      arma::chol(R,mx);
      minvC = mx.t()* inv_sympd((mx*mx.t()));

      return Rcpp::wrap(minvC); 
      
    }else {
      throw std::runtime_error("Matrix must be symmetric");
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
}


// [[Rcpp::export]]
Rcpp::RObject bdSVD (const Rcpp::RObject & x, int k, int nev )
{
  auto dmtype = beachmat::find_sexp_type(x);
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd X;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP ) {
    auto dmat = beachmat::create_integer_matrix(x);
    ncols = dmat->get_ncol();
    nrows = dmat->get_nrow();
    
    X.resize(nrows, ncols);
    Rcpp::IntegerVector output(dmat->get_nrow());
    
    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, output.begin());
      X.col(ccol) = Rcpp::as<Eigen::VectorXd >(output);
    }

  } else if (dmtype==REALSXP) {
    
    auto dmat = beachmat::create_numeric_matrix(x);
    ncols = dmat->get_ncol();
    nrows = dmat->get_nrow();
    X.resize(nrows, ncols);

    Rcpp::NumericVector output(dmat->get_nrow());
    //Eigen::MatrixXd X(nrows,ncols);
    
    for (size_t ccol=0; ccol<ncols; ++ccol) {
      dmat->get_col(ccol, output.begin());
      X.col(ccol) = Rcpp::as<Eigen::VectorXd >(output);
    }

  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  Eigen::MatrixXd Xtcp = CrossProduct(Normalize_Data(Rcpp::wrap(X)),true);
  Spectra::DenseSymMatProd<double> op(Xtcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, nev);
  
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  
  if(eigs.info() == Spectra::SUCCESSFUL)
  {
    ret["d$"] = eigs.eigenvalues().cwiseSqrt();
    ret["u$"] = eigs.eigenvectors();
  }
  
  Eigen::MatrixXd Xcp = CrossProduct(Normalize_Data(Rcpp::wrap(X)),false); // Normalitzem matriu
  // Eigen::MatrixXd Xcp = CrossProduct_eig(X,false); // No normalitzem matriu - matriu normalitzada prèviament
  
  Spectra::DenseSymMatProd<double> opv(Xcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigsv(&opv, k, nev);
  
  // Initialize and compute
  eigsv.init();
  nconv = eigsv.compute();
  
  // Retrieve results
  if(eigsv.info() == Spectra::SUCCESSFUL)
  {
    ret["v$"] = eigsv.eigenvectors();
  }
  
  return Rcpp::wrap(ret);
  
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

  val <- matrix(sample(1:100,1000, replace = TRUE), ncol = 100)
  da_mat <- DelayedArray(seed = val) ; da_mat
  Z <- da_mat # DelayedArray
  X <- val    # NO DelayedArray
  
  X <- matrix(runif(10000), ncol=100)
  Z <- DelayedArray(X)  # DelayedArray

  library(microbenchmark)
  
  # Proves : X*X' i X'*X
  X <- matrix(runif(10000), ncol=100)
  Z <- DelayedArray(X)  # DelayedArray
  
  results <- microbenchmark(crosR <- crossprod(X), # Funció R
                            crosArma <- BDCrossProduct(Z,FALSE), # DelayedArray
                            crosArmaND <- BDCrossProduct_benchmat_nodelayed(X,FALSE), # NO DelayedArray
                            tcrosR <- tcrossprod(X), # Funció R
                            tcrosArma <- BDCrossProduct(Z,TRUE), # DelayedArray
                            tcrosArmaND <- BDCrossProduct_benchmat_nodelayed(X,TRUE), # NO DelayedArray
                            times=50L)
  
  print(summary(results)[, c(1:7)],digits=1)
  
  all.equal(crossprod(X), crosArmaND)
  all.equal(crossprod(X), crosArma)
  all.equal(tcrossprod(X), tcrosArmaND)
  all.equal(tcrossprod(X), tcrosArma)
  
  
  # Proves : Matriu Inversa - Cholesky
  
  Z <- matrix(sample(1:100, 10000, replace=TRUE), nrow=100)
  Z[lower.tri(Z)] <- t(Z)[lower.tri(Z)]
  Z <- DelayedArray(Z);Z
  
  results <- microbenchmark(invR <- solve(Z),  # Inversa amb Solve
                              invCpp <- BDCInverse_Cholesky(Z))  # Inversa amb Cholesky
  print(summary(results)[, c(1:7)],digits=3)

  all.equal(invR, invCpp)  
  
  # Proves : Descomposició SVD Spectra
  
  library(RSpectra)
  
  set.seed(5283)
  A <- matrix(sample(1:100, 25, replace=TRUE), nrow=5);A
  Ad <- DelayedArray(seed = A) ; Ad
  results <- microbenchmark(svdR <- svds(scale(Ad,center = TRUE, scale = TRUE ), 2),  # SVD amb RSpectra - Passem matriu normalitzada
                            svdCpp <- bdSVD(Ad, 2, 5), 
                            times = 10)  # SVD amb Spectra - passem matriu sense normalitzar
  print(summary(results)[, c(1:7)],digits=3)
  
  set.seed(5283)
  A <- matrix(sample(1:100, 1000000, replace=TRUE), ncol = 125)
  Ad <- DelayedArray(seed = A) ; Ad
  results <- microbenchmark(svdR <- svds(scale(A,center = TRUE, scale=TRUE), 2, ncv = 100),  # SVD amb RSpectra - Passem matriu normalitzada
                            svdCpp <- bdSVD(Ad, 2, 100),
                            times=10L)  # SVD amb Spectra - passem matriu sense normalitzar
  print(summary(results)[, c(1:7)],digits=3)
  

  A <- matrix(runif(1000000), ncol=500)
  Ad <- DelayedArray(A)  # DelayedArray
  results <- microbenchmark(svdR <- svds(scale(Ad,center = TRUE, scale=TRUE), 2, ncv = 100),  # SVD amb RSpectra - Passem matriu normalitzada
                            svdCpp <- bdSVD(Ad, 2, 100),
                                    times=10L)  # SVD amb Spectra - passem matriu sense normalitzar
  print(summary(results)[, c(1:7)],digits=3)
  
  svdR$d
  svdCpp$d
  
  head(svdR$u)
  head(svdCpp$u)
  
  head(svdR$v)
  head(svdCpp$v)
  
*/
