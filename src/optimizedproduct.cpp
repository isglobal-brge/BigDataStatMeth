#include "include/optimizedproduct.h"


// Computes CrossProduct X'X
Eigen::MatrixXd bdcrossproduct (Eigen::MatrixXd& X)
{
  size_t nc(X.cols());
  Eigen::MatrixXd XtX(Eigen::MatrixXd(nc, nc).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint()));
  return(XtX);
}

// Compute tCrossProduct XX'
Eigen::MatrixXd bdtcrossproduct (Eigen::MatrixXd& X)
{
  
  size_t nr(X.rows());
  Eigen::MatrixXd XXt(Eigen::MatrixXd(nr, nr).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X));
  return(XXt);
}


// Compute weighted crossproduct XwX'
Eigen::MatrixXd xwxt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w) 
{
  const int n(X.rows());
  Eigen::MatrixXd XwXt(Eigen::MatrixXd(n, n).setZero().
                         selfadjointView<Eigen::Lower>().rankUpdate(X * w.array().sqrt().matrix().asDiagonal()));
  return (XwXt);
}


// Compute transposed weighted crossproduct X'wX 
Eigen::MatrixXd xtwx(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w)
{
  const int n(X.cols());
  Eigen::MatrixXd XtwX(Eigen::MatrixXd(n, n).setZero().
                         selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint() * w.array().sqrt().matrix().asDiagonal()));
  return (XtwX);
}


// Compute weighted crossproduct Xw
Eigen::MatrixXd Xw(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w) 
{
  Eigen::MatrixXd Xw = X * w.array().matrix().asDiagonal();
  return (Xw);
}


// Compute weighted crossproduct Xw
Eigen::MatrixXd wX(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w) 
{
  Eigen::MatrixXd wX = w.array().matrix().asDiagonal()*X;
  return (wX);
}



// Matirx - vector as diagonal matrix Multiplication
Eigen::MatrixXd Xwd(const Eigen::MatrixXd& X, const Eigen::VectorXd& w)
{
  int n = X.rows();
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n,X.cols()) ; 

  for (int i=0; i<n; i++) {
    C.col(i) = X.col(i)*w(i);
  }
  return(C);
}


//vector as diagonal matrix -  Matirx Multiplication
Eigen::MatrixXd wdX(const Eigen::MatrixXd& X, const Eigen::VectorXd& w)
{
  int n = X.cols();
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(X.rows(),n) ; 
  
  for (int i=0; i<n; i++) {
    C.row(i) = w(i)*X.row(i);
  }
  return(C);
}


// Matirx - vector as diagonal matrix Multiplication (parallel)
Eigen::MatrixXd Xwd_parallel(const Eigen::MatrixXd& X, const Eigen::VectorXd& w, Rcpp::Nullable<int> threads = R_NilValue)
{
  int n = X.rows();
  unsigned int ithreads;
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n,X.cols()) ; 
  
  if(threads.isNotNull()) 
  {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
      ithreads = Rcpp::as<int> (threads);
    else 
      ithreads = std::thread::hardware_concurrency()/2;
  }
  else    ithreads = std::thread::hardware_concurrency() / 2; //omp_get_max_threads();
  
  omp_set_num_threads(ithreads);

#pragma omp parallel shared(X, w, C) 
{
#pragma omp for schedule (dynamic)
  for (int i=0; i<n; i++)
  {
    C.col(i) = X.col(i)*w(i);
  }  
}
  return(C);
}


//vector as diagonal matrix -  Matirx Multiplication (parallel)
Eigen::MatrixXd wdX_parallel(const Eigen::MatrixXd& X, const Eigen::VectorXd& w, Rcpp::Nullable<int> threads = R_NilValue)
{
  int n = X.cols();
  unsigned int ithreads;
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(X.rows(),n);
  
  if(threads.isNotNull()) 
  {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
      ithreads = Rcpp::as<int> (threads);
    else 
      ithreads = std::thread::hardware_concurrency()/2;
  }
  else    ithreads = std::thread::hardware_concurrency() /2; //omp_get_max_threads();
  
  omp_set_num_threads(ithreads);
  
#pragma omp parallel shared(X, w, C) 
{
#pragma omp for schedule (dynamic)
  for (int i=0; i<n; i++)
  {
    C.row(i) = w(i)*X.row(i);
  }
}
  return(C);
}



//' Crossproduct and transposed crossproduct of DelayedArray
//' 
//' This function performs a crossproduct or transposed crossproduct of numerical or DelayedArray matrix.
//' 
//' @param a numerical or Delayed Array matrix
//' @param transposed (optional, default = false) boolean indicating if we have to perform a crossproduct (transposed=false) or transposed crossproduct (transposed = true)
//' @return numerical matrix with crossproduct or transposed crossproduct 
//' @examples
//' 
//' library(DelayedArray)
//' 
//' n <- 100
//' p <- 60
//' 
//' X <- matrix(rnorm(n*p), nrow=n, ncol=p)
//' 
//' # without DelayedArray
//' bdcrossprod(X)
//' bdcrossprod(X, transposed = TRUE)
//' 
//' # with DelayedArray
//' XD <- DelayedArray(X)
//' bdcrossprod(XD)
//' bdcrossprod(XD, transposed = TRUE)
//' 
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd bdcrossprod(Rcpp::RObject a, Rcpp::Nullable<Rcpp::Function> transposed = R_NilValue)
{
  
  Eigen::MatrixXd A;
  bool btrans;
  
  if(transposed.isNotNull())
    btrans = Rcpp::as<bool> (transposed);
  else
    btrans = FALSE;

  // Read DelayedArray's a and b
  if ( a.isS4() == true)    
  {
    A = read_DelayedArray(a);
  } else {
    try{  
      A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
    }
    catch(std::exception &ex) { }
  }

  if(btrans == true) {
    return(bdtcrossproduct(A)) ;
  }else {
    return(bdcrossproduct(A));
  }
}


// Weighted product
//' Matrix - Weighted vector Multiplication with numerical or DelayedArray data
//' 
//' This function performs a weighted product of a matrix(X) with a weighted diagonal matrix (w)
//' 
//' @param a numerical or Delayed Array matrix
//' @param w vector with weights
//' @param op string indicating if operation 'xtwx' and 'xwxt' for weighted cross product (Matrix - Vector - Matrix) or 'Xw' and 'wX' for weighted product (Matrix - Vector)
//' @return numerical matrix 
//' @examples
//' 
//' library(DelayedArray)
//' 
//' n <- 100
//' p <- 60
//' 
//' X <- matrix(rnorm(n*p), nrow=n, ncol=p)
//' 
//' u <- runif(n)
//' w <- u * (1 - u)
//' ans <- bdwproduct(X, w,"xtwx")
//' 
//' # with Delayed Array
//' 
//' DX <- DelayedArray(X)
//' 
//' ans <- bdwproduct(DX, w,"xtwx")
//' 
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd bdwproduct(Rcpp::RObject a, Rcpp::RObject w, std::string op)
{
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd W;
  
  
  // Read DelayedArray's a and b
  if ( a.isS4() == true)    
  {
    A = read_DelayedArray(a);
  } else {
    try{  
      A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
    }
    catch(std::exception &ex) { }
  }
  
  if ( w.isS4() == true)   
  {
    W = read_DelayedArray(w);
  }  else {
    W = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(w);
  } 
  
  if(op == "xwxt") {
    return(xwxt(A,W)) ;
  }else if (op == "xtwx") {
    return(xtwx(A,W));
  }else if (op == "Xw") {
    return(Xw(A,W));
  }else if (op == "wX") {
    return(wX(A,W));
  } else
  {
    throw("Invalid option, valid options : 'xtwx' and 'xwxt' for weighted cross product or 'Xw' and 'wX' for weighted product");
  }
}

// Weighted product
//' Matrix - Weighted Scalar Multiplication with numerical or DelayedArray data
//' 
//' This function performs a weighted product of a matrix(X) with a weighted diagonal matrix (w)
//' 
//' @param a numerical or Delayed Array matrix
//' @param w scalar, weight
//' @param op string indicating if operation  "Xw" or "wX"
//' @return numerical matrix 
//' @examples
//' 
//' library(DelayedArray)
//' 
//' n <- 100
//' p <- 60
//' 
//' X <- matrix(rnorm(n*p), nrow=n, ncol=p)
//' w <- 0.75
//' 
//' bdScalarwproduct(X, w,"Xw")
//' bdScalarwproduct(X, w,"wX")
//' 
//' # with Delayed Array
//' 
//' DX <- DelayedArray(X)
//' 
//' bdScalarwproduct(DX, w,"Xw")
//' bdScalarwproduct(DX, w,"wX")
//' 
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd bdScalarwproduct(Rcpp::RObject a, double w, std::string op)
{
  
  Eigen::MatrixXd A;

  // Read DelayedArray's a and b
  if ( a.isS4() == true)    
  {
    A = read_DelayedArray(a);
  } else {
    try{  
      A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
    }
    catch(std::exception &ex) { }
  }
  
  std::transform(op.begin(), op.end(), op.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  
  if (op == "xw" || op == "wx") {
    return(w*A);
  } else  {
    throw("Invalid option (valid options : xtwx or xwxt)");
  }
}




/***R
n <- 10
p <- 6
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  
u <- runif(n)
w <- u * (1 - u)
# OK

bdwproduct(X, w,"xtwx",)
crossprod(X, w*X)
t(X)%*%diag(w)%*%X

stopifnot(all.equal(bdwproduct(X, w,"xtwx"),crossprod(X, w*X)),
          all.equal(bdwproduct(X, w,"xtwx"),t(X)%*%diag(w)%*%X))



n <- 10
p <- 6
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
XD <- DelayedArray(X)
ru <- runif(p)
w <- ru * (1 - ru)
wD <- DelayedArray(as.matrix(w))
  
bdwproduct(X, w,"xwxt")
bdwproduct(XD, w,"xwxt")
X%*%diag(w)%*%t(X)

stopifnot(all.equal(bdwproduct(X, w,"xwxt"), X%*%diag(w)%*%t(X) ),
          all.equal(bdwproduct(XD, w,"xwxt"), X%*%diag(w)%*%t(X) ),
          all.equal(bdwproduct(XD, wD,"xwxt"), X%*%diag(w)%*%t(X) )) 



XD <- DelayedArray(X)
bdcrossprod(X)
bdcrossprod(XD)

stopifnot(all.equal(bdcrossprod(X),crossprod(X)),
          all.equal(bdcrossprod(XD),crossprod(X)))

stopifnot(all.equal(bdcrossprod(X, transposed = TRUE),tcrossprod(X)),
          all.equal(bdcrossprod(XD, transposed = TRUE),tcrossprod(X)))

bdcrossprod(X,TRUE)
bdcrossprod(XD,TRUE)

n <- 500
A <- matrix(rnorm(n*n), nrow=n, ncol=n)
DA <- DelayedArray(A)

A[1:5,1:5]
DA[1:5,1:5]
  
cpA <- bdcrossprod(A, transposed = FALSE)  

cpDA <- bdcrossprod(DA, transposed = TRUE)  # With DelayedArray data type


cpDA <- bdcrossprod(DA, transposed = FALSE)  # With DelayedArray data type
  
all.equal(cpDA, tcrossprod(A)) #
all.equal(cpDA, crossprod(A))



all.equal(cpDA, A%*%t(A))
all.equal(cpA, A%*%t(A))




*/