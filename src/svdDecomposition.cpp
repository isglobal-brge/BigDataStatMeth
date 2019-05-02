#include "include/svdDecomposition.h"
#include "include/ReadDelayedData.h"



svdeig RcppBDsvd_eig( Eigen::MatrixXd& X, int k, int nev, bool normalize )
{

  svdeig retsvd;
  Eigen::MatrixXd nX;
  
  if( k==0 )    k = (std::min(X.rows(), X.cols()))-1;
  if(nev == 0)  nev = k + 1 ;
  
  Eigen::MatrixXd Xtcp;
  if(normalize ==true )  {
    // Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(RcppNormalize_Data(X))));
    nX = RcppNormalize_Data(X);
    Xtcp =  bdtcrossproduct(nX);
  }else {
    //Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(X)));
    Xtcp =  bdtcrossproduct(X);

  }
  
  Spectra::DenseSymMatProd<double> op(Xtcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, nev);
  
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  
  if(eigs.info() == Spectra::SUCCESSFUL)
  {
    retsvd.d = eigs.eigenvalues().cwiseSqrt();
    retsvd.u = eigs.eigenvectors();
  }
  
  Eigen::MatrixXd Xcp;
  if(normalize ==true )  {
    // Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(RcppNormalize_Data(X))));  
    Xcp =  bdcrossproduct(nX);  
  }else {
    // Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(X)));
    Xcp =  bdcrossproduct(X);
  }
  
  Spectra::DenseSymMatProd<double> opv(Xcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigsv(&opv, k, nev);
  
  // Initialize and compute
  eigsv.init();
  nconv = eigsv.compute();

  // Retrieve results
  if(eigsv.info() == Spectra::SUCCESSFUL)
  {
    retsvd.v = eigsv.eigenvectors();
  }
  
  return retsvd;
  
}




svd RcppBDsvd ( const Rcpp::NumericMatrix& X, int k, int nev, bool normalize )
{
  
  svd retsvd;
  if( k==0 )    k = (std::min(X.rows(), X.cols()))-1;
  if(nev == 0)  nev = k + 1 ;
  
  Eigen::MatrixXd Xtcp;
  if(normalize ==true )  {
    Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( RcppNormalize_Data_r(X)));
  }else {
    Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( X));
  }
  
  Spectra::DenseSymMatProd<double> op(Xtcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, nev);
  
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  
  if(eigs.info() == Spectra::SUCCESSFUL)
  {
    retsvd.d = Rcpp::wrap(eigs.eigenvalues().cwiseSqrt());
    retsvd.u = Rcpp::wrap(eigs.eigenvectors());
  }
  
  Eigen::MatrixXd Xcp;
  if(normalize ==true )  {
    Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( RcppNormalize_Data_r(X)));  
  }else {
    Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( X));
  }
  
  Spectra::DenseSymMatProd<double> opv(Xcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigsv(&opv, k, nev);
  
  // Initialize and compute
  eigsv.init();
  nconv = eigsv.compute();
  
  // Retrieve results
  if(eigsv.info() == Spectra::SUCCESSFUL)
  {
    retsvd.v = Rcpp::wrap(eigsv.eigenvectors());
  }
  
  return retsvd;
  
}


svdeig RcppCholDec(const Eigen::MatrixXd& X)
{
  Eigen::MatrixXd mX = X;
  svdeig decomp;
  
  Eigen::LDLT<Eigen::MatrixXd> cholSolv = mX.ldlt();

  // if symetric + positive definite -> info = Success
  // else no Cholesky Decomposition --> svd decomposition with Expectra
  if(cholSolv.info()==Eigen::Success)
  {
    size_t n = cholSolv.cols();
    decomp.d = cholSolv.vectorD();
    Eigen::MatrixXd preinv = Eigen::MatrixXd::Identity(n, n);
    decomp.v = cholSolv.solve(preinv);
  } else {
    Rcpp::Rcout<<"No symetric positive matrix, Cholesky decomposition not viable.";
    decomp = RcppBDsvd_eig(mX, int(), int(), false);
  }
  return(decomp);
}








// [[Rcpp::export]]
Eigen::MatrixXd bdInvCholesky (const Rcpp::RObject & x )
{
  
  svdeig result;
  Eigen::MatrixXd X;
  
  if ( x.isS4() == true)    
  {
    X = read_DelayedArray(x);
  } else {
    try{  
      X = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(x);
    }
    catch(std::exception &ex) { }
  }
  
  result = RcppCholDec(X);
  
  /*
  if(x.sexp_type()==13)
  {
    Eigen::MatrixXi eX(Rcpp::as<Eigen::Map<Eigen::MatrixXi> >(x));
    result = RcppCholDec(eX.cast<double>());
  }
  else if (x.sexp_type()==14)
  {
    Eigen::MatrixXd eX(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(x));
    result = RcppCholDec(eX);
  }
  */
  
  return (result.v);
  
}

/*
// [[Rcpp::export]]
Rcpp::RObject BDsvd (const Rcpp::RObject & x, int k, int nev, bool normalize )
{

  auto dmtype = beachmat::find_sexp_type(x);

  // size_t ncols = 0, nrows=0;
  Eigen::MatrixXd X;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( x.isS4() == true){
      X = read_DelayedArray(x);
    }else {
      try{
        X = Rcpp::as<Eigen::MatrixXd >(x);
      }catch(std::exception &ex) {
        X = Rcpp::as<Eigen::VectorXd >(x);
      }
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  Eigen::MatrixXd Xtcp;

  if( normalize == true)  {
    Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(RcppNormalize_Data(X))));  
  } else {
    Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(X)));  
  }
  
  
  Spectra::DenseSymMatProd<double> op(Xtcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, nev);
  
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  

  if(eigs.info() == Spectra::SUCCESSFUL)
  {
    ret["d$"] = eigs.eigenvalues();
    ret["u$"] = eigs.eigenvectors();
  }
  
  Eigen::MatrixXd Xcp;
  
  if(normalize==true )  {
    Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(RcppNormalize_Data(X))));
  } else
  {
    Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(X)));
  }

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

*/




// [[Rcpp::export]]
Rcpp::RObject BDsvd (const Rcpp::RObject & x, Rcpp::Nullable<int> k=0, Rcpp::Nullable<int> nev=0, Rcpp::Nullable<bool> normalize=true )
{
  
  auto dmtype = beachmat::find_sexp_type(x);
  int ks, nvs;
  bool bnorm;
  
  if(k.isNull())  ks = 0 ;
  else    ks = Rcpp::as<int>(k);
  
  if(nev.isNull())  nvs = 0 ;
  else    nvs = Rcpp::as<int>(nev);
  
  if(normalize.isNull())  bnorm = true ;
  else    bnorm = Rcpp::as<bool>(normalize);
  
  // size_t ncols = 0, nrows=0;
  Eigen::MatrixXd X;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( x.isS4() == true){
      X = read_DelayedArray(x);
    }else {
      try{
        X = Rcpp::as<Eigen::MatrixXd >(x);
      }catch(std::exception &ex) {
        X = Rcpp::as<Eigen::VectorXd >(x);
      }
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  svdeig retsvd;
  retsvd = RcppBDsvd_eig(X,ks,nvs,bnorm);
  
  ret["u"] = retsvd.u;
  ret["v"] = retsvd.v;
  ret["d"] = retsvd.d;
  
  return Rcpp::wrap(ret);
  
}


/***R

library(microbenchmark)
library(DelayedArray)

# Proves : Matriu Inversa - Cholesky
set.seed(12)
n <- 10
p <- 10
Z <- matrix(rnorm(n*p), nrow=n, ncol=p)
X <- Z
Z[lower.tri(Z)] <- t(Z)[lower.tri(Z)]

# Z <- DelayedArray(Z);Z

A <- matrix(c(5,0,2,5,0,5,-1,-1,2,-1,5,-1,5,-1,-1,5),byrow = TRUE, nrow = 4)

A <- Posdef(n=100, ev=1:100)
AD <- DelayedArray(A)
invCpp <- bdInvCholesky_LDL_eigen(A); invCpp[1:5,1:5]
invR <- solve(A); invR[1:5,1:5]
invP <- inversechol_par(A);invP[1:5,1:5]

stopifnot(all.equal(invCpp,invR ))
stopifnot(all.equal(invCpp,invR ))

invCpp <- bdInvCholesky_LDL_eigen(AD); invCpp[1:5,1:5]
invP <- inversechol_par(AD);invP[1:5,1:5]



invCpp <- bdInvCholesky_LDL(A); invCpp$v[1:4,1:4]
invR <- solve(A); invR[1:4,1:4]
invP <- inversechol(A);invP[1:4,1:4]

stopifnot(all.equal(solve(Z),bdInvCholesky_LDL(Z)$v ))






  results <- microbenchmark(invR <- solve(Z),  # Inversa amb Solve
                            invCpp <- bdInvCholesky_LDL(Z),
                            invCpppar <- inversechol(Z),
                            times = 5L)  # Inversa amb Cholesky
                            
  print(summary(results)[, c(1:7)],digits=3)
  
  all.equal(invR, invCpppar)  
  
  
  cc <- eigen(tcrossprod(X))
  sqrt(cc$values[1:10])
  sqrt(invCpp$d)
  invCpp$v  
  cc$vectors
  
  
  
  A <- matrix(c(5,-3,4,-3,3,-4,4,-4,6), byrow = TRUE, ncol = 3)
  A <- matrix(c(2,-1,0,-1,2,-1,0,-1,1), byrow = TRUE, ncol = 3)
  A <- matrix(c(3,4,3,4,8,6,3,6,9), byrow = TRUE, ncol = 3)
  
  invR <- solve(A) ; invR
  eigen(tcrossprod(A))
  invCpp <- bdInvCholesky(A); invCpp
  
  LOOE_BLAST()
  
  n <- 10
  p <- 10
  Z <- matrix(rnorm(n*p), nrow=n, ncol=p)
  
  a <- BDsvd(Z,9,10,FALSE )
  b <- eigen(tcrossprod(Z))
  a$d^2
  b$values
  ;
  
  n <- 500
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  AD <- DelayedArray(A)
  
  dim(A)
  
  
  res <- microbenchmark( bdsvd <- BDsvd( A, n-1, n, FALSE), # No normalitza la matriu
                         bdsvdD <- BDsvd( AD, n-1, n, FALSE), # No normalitza la matriu
                         sbd <- svd(tcrossprod(A)),
                         times = 5, unit = "s")
  
  print(summary(res)[, c(1:7)],digits=3)
  
  sqrt(sbd$d[1:10])
  bdsvd$d[1:10]
  rsv$d[1:10]
  rsv <- rsvd::rsvd(A)
  
  bdsvd$u[1:5,1:5]
  
  svd(tcrossprod(A))$d[1:10]
  
*/