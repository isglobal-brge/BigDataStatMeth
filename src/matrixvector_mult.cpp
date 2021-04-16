#include "include/matrixvector_mult.h"

using namespace std;


//' Matrix - Weighted vector Multiplication with numerical or DelayedArray data
//' 
//' This function performs a weighted product of a matrix(X) with a weighted diagonal matrix (w)
//' 
//' @param X numerical or Delayed Array matrix
//' @param W vector with weights
//' @param op string indicating if operation  "Xw" , "wX" or "wXw"
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
//' 
//' # by columnes
//' bdwXw(X, w,"wX") 
//' 
//' # with Delayed Array
//' 
//' DX <- DelayedArray(X)
//' 
//' bdwXw(DX, w,"wX")
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdwXw(Rcpp::RObject X, Rcpp::RObject W, std::string op)
{
  
  //.commented 20201120 - warning check().// char Nchar='N';
  //.commented 20201120 - warning check().// char Uchar='T';
  //.commented 20201120 - warning check().// char Lchar='C';
  //.commented 20201120 - warning check().// double done = 1.0;
  //.commented 20201120 - warning check().// int ione = 1;
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd w;  // column or row vector
  
  //.commented 20201120 - warning check().// bool bpar=true;
  
  // Read DelayedArray's X and b
  
  if ( X.isS4() == true)
  {
    A = read_DelayedArray(X);
  } else {
    try{  
      A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
    }
    catch(std::exception &ex) { }
  }
  
  // Read DelayedArray's X and b
  if ( W.isS4() == true)
  {
    w = read_DelayedArray(W);
  } else {
    try{  
      w = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(W);
    }
    catch(std::exception &ex) { }
  }
  
  if(op=="wX")
  {
    if(w.rows() == A.rows())
      return(Rcpp::wrap(w.asDiagonal()*A ));
    else
      throw std::runtime_error("Dimension error");
  }
  else if(op=="Xw")  
  {
    if(w.rows() == A.cols())
      return(Rcpp::wrap(A*w.asDiagonal()));
    else
      throw std::runtime_error("Dimension error");
    
  }
  else if(op == "wXw")
  {
    if(w.rows() == A.rows())
      return(Rcpp::wrap(w.asDiagonal()*A*w.asDiagonal()));
    else
      throw std::runtime_error("Dimension error");
  }else
  {
    throw std::runtime_error("Operation type error");
  }
  
}

