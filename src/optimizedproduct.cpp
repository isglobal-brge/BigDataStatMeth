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


// Compute weighted crossproduct XwXi
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





// [[Rcpp::export]]
Eigen::MatrixXd bdcrossprod(Rcpp::RObject a, Rcpp::Nullable<Rcpp::Function> transposed = R_NilValue)
{
  
  Eigen::MatrixXd A;
  bool btrans;
  
  if(transposed.isNotNull())
    btrans = transposed;
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
  } else
  {
    throw("Invalid option (valid options : xtwx or xwxt)");
  }
}



/***R
n <- 100
p <- 60
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  
u <- runif(n)
w <- u * (1 - u)
# OK

bdwproduct(X, w,"xtwx")
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

*/