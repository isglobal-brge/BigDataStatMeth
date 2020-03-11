#include "include/svdDecomposition.h"
#include "include/ReadDelayedData.h"

// SVD decomposition 
svdeig RcppbdSVD( Eigen::MatrixXd& X, int k, int nev, bool bcenter, bool bscale )
{
  
  svdeig retsvd;
  Eigen::MatrixXd nX;
  int nconv;
  
  if( k==0 )    k = (std::min(X.rows(), X.cols()))-1;
  else if (k > (std::min(X.rows(), X.cols()))-1 ) k = (std::min(X.rows(), X.cols()))-1;
  
  if(nev == 0)  nev = k + 1 ;
  if(nev<k) nev = k + 1;
  
  {
    Eigen::MatrixXd Xtcp;
    //..//if(normalize ==true )  {
    if(bcenter ==true || bscale == true)  {
      // Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(RcppNormalize_Data(X))));
      nX = RcppNormalize_Data(X, bcenter, bscale);
      Xtcp =  bdtcrossproduct(nX);
    }else {
      //Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(X)));
      Xtcp =  bdtcrossproduct(X);
      
    }
    
    Spectra::DenseSymMatProd<double> op(Xtcp);
    Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, nev);
    
    // Initialize and compute
    eigs.init();
    nconv = eigs.compute();
    
    if(eigs.info() == Spectra::SUCCESSFUL)
    {
      retsvd.d = eigs.eigenvalues().cwiseSqrt();
      retsvd.u = eigs.eigenvectors();
      retsvd.bokuv = true;
    } else {
      retsvd.bokuv = false;
    }
    
  }
  
  if(retsvd.bokuv == true)
  {
    Eigen::MatrixXd Xcp;
    if(bcenter ==true || bscale==true )  {
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
    } else {
      retsvd.bokd = false;
    }
    
  }
  
  return retsvd;
}



svdeig RcppCholDec(const Eigen::MatrixXd& X)
{
  Eigen::MatrixXd mX = X;
  svdeig decomp;
  
  Eigen::LDLT<Eigen::MatrixXd> cholSolv = mX.ldlt();

  // if symetric + positive definite -> info = Success
  // else no Cholesky Decomposition --> svd decomposition with Spectra
  if(cholSolv.info()==Eigen::Success)
  {
    size_t n = cholSolv.cols();
    decomp.d = cholSolv.vectorD();
    Eigen::MatrixXd preinv = Eigen::MatrixXd::Identity(n, n);
    decomp.v = cholSolv.solve(preinv);
  } else {
    Rcpp::Rcout<<"No symetric positive matrix, Cholesky decomposition not viable.";
    decomp = RcppbdSVD(mX, int(), int(), false);
  }
  return(decomp);
}







//' Inverse Cholesky of Delayed Array
//' 
//' This function get the inverse of a numerical or Delayed Array matrix. If x is hermitian and positive-definite matrix then 
//' performs get the inverse using Cholesky decomposition
//' 
//' 
//' @param x numerical or Delayed Array matrix. If x is Hermitian and positive-definite performs
//' @return inverse matrix of d 
//' @examples
//' 
//' A <- matrix(c(3,4,3,4,8,6,3,6,9), byrow = TRUE, ncol = 3)
//' bdInvCholesky(A)
//' 
//' # with Delayed Array
//' DA <- DelayedArray(A)
//' bdInvCholesky(DA)
//' 
//' @export
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
  
  return (result.v);
  
}



//' SVD of DelayedArray 
//' 
//' This function performs a svd decomposition of numerical matrix or Delayed Array
//' 
//' @param x numerical or Delayed Array matrix
//' @param k number of eigen values , this should satisfy k = min(n, m) - 1
//' @param nev (optional, default nev = n-1) Number of eigenvalues requested. This should satisfy 1≤ nev ≤ n, where n is the size of matrix. 
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
//' @return u eigenvectors of AA^t, mxn and column orthogonal matrix
//' @return v eigenvectors of A^tA, nxn orthogonal matrix
//' @return d singular values, nxn diagonal matrix (non-negative real values)
//' @examples
//' n <- 500
//' A <- matrix(rnorm(n*n), nrow=n, ncol=n)
//' AD <- DelayedArray(A)
//' 
//' # svd without normalization
//' bdSVD( A, bscale = FALSE, bcenter = FALSE ), # No matrix normalization
//' decsvd$d
//' decsvd$u
//' 
//' # svd with normalization
//' decvsd <- bdSVD( A, bscale = TRUE, bcenter = TRUE), # Matrix normalization
//' 
//' decsvd$d
//' decsvd$u
//' 
//' # svd with scaled matrix (sd)
//' decvsd <- bdSVD( A, bscale = TRUE, bcenter = FALSE), # Scaled matrix
//' 
//' decsvd$d
//' decsvd$u

//' # svd with centered matrix (sd)
//' decvsd <- bdSVD( A, bscale = FALSE, bcenter = TRUE), # Centered matrix
//' decsvd$d
//' decsvd$u

//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSVD (const Rcpp::RObject & x, Rcpp::Nullable<int> k=0, Rcpp::Nullable<int> nev=0,
                     Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true)
{
  
  auto dmtype = beachmat::find_sexp_type(x);
  int ks, nvs;
  bool bnorm, bcent, bscal;
  
  if(k.isNull())  ks = 0 ;
  else    ks = Rcpp::as<int>(k);
  
  if(nev.isNull())  nvs = 0 ;
  else    nvs = Rcpp::as<int>(nev);
  
  
  if(bcenter.isNull())  bcent = true ;
  else    bcent = Rcpp::as<bool>(bcenter);
  
  if(bscale.isNull())  bscal = true ;
  else    bscal = Rcpp::as<bool>(bscale);
  
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
  retsvd = RcppbdSVD(X,ks,nvs, bcent, bscal);
  
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
  p <- 20
  Z <- matrix(rnorm(n*p), nrow=n, ncol=p)
  Zn <- scale(Z,center = TRUE, scale = TRUE)
  
  library(BigDataStatMeth)
  a <- bdSVD(Z,8,10, TRUE, TRUE)
  b <- eigen(tcrossprod(Z))
  
  svd(Zn)$d
  a$d
  
  svd(Z)$d^2
  a$d^2
  b$values
  
  
  
  a$d
  a$okd
  
  a$d^2
  b$values
  
  ;
  
  n <- 1000
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  AD <- DelayedArray(A)
  
  dim(A)
  
  
  res <- microbenchmark( bdsvd <- bdSVD( A, n-1, n, FALSE), # No normalitza la matriu
                         bdsvdD <- bdSVD( AD, n-1, n, FALSE), # No normalitza la matriu
                         sbd <- svd(tcrossprod(A)),
                         times = 5, unit = "s")
  
  print(summary(res)[, c(1:7)],digits=3)
  
  sqrt(sbd$d[1:10])
  bdsvd$d[1:10]
  rsv$d[1:10]
  rsv <- rsvd::rsvd(A)
  
  bdsvd$u[1:5,1:5]
  
  svd(tcrossprod(A))$d[1:10]
  
  A <- matrix(c(5,-3,4,-3,3,-4,4,-4,6,7,2,3), byrow = TRUE, ncol = 3)
  svd(A)$u
  
  bdSVD(A)$u
  
  
*/