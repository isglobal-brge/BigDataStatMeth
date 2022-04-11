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
  
  if(threads.isNotNull()) {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
      ithreads = Rcpp::as<int> (threads);
    } else {
      ithreads = getDTthreads(0, true);
      //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;}
    }
  } else {
    ithreads = getDTthreads(0, true);
    //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;
  }
  
  //.OpenMP.//omp_set_num_threads(ithreads);

  //.OpenMP.//#pragma omp parallel shared(X, w, C) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(X, w, C) 
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
  
  if(threads.isNotNull()) {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
      ithreads = Rcpp::as<int> (threads);
    } else {
      ithreads = getDTthreads(0, true);
      //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;}
    }
  } else {
    ithreads = getDTthreads(0, true);
    //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;
  }
  
  //.OpenMP.// omp_set_num_threads(ithreads);
  
  //.OpenMP.//#pragma omp parallel shared(X, w, C) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(X, w, C) 
{
#pragma omp for schedule (dynamic)
  for (int i=0; i<n; i++)
  {
    C.row(i) = w(i)*X.row(i);
  }
}
  return(C);
}




// [[Rcpp::export(.bdCrossprod_generic)]]
Eigen::MatrixXd bdCrossprod_generic( Rcpp::RObject A, Rcpp::Nullable<Rcpp::RObject> B=  R_NilValue, 
                                     Rcpp::Nullable<bool> transposed = R_NilValue,
                                     Rcpp::Nullable<int> block_size = R_NilValue, 
                                     Rcpp::Nullable<bool> paral = R_NilValue,
                                     Rcpp::Nullable<int> threads = R_NilValue )
{
  
  Eigen::MatrixXd mA;
  Eigen::MatrixXd mB;
  Eigen::MatrixXd C;
  
  bool btrans, bparal;
  int iblock_size;
  
  if(transposed.isNotNull())
    btrans = Rcpp::as<bool> (transposed);
  else
    btrans = FALSE;
  
  if( paral.isNull()) {
    bparal = false;
  } else {
    bparal = Rcpp::as<bool> (paral);
  }
  
  
  
  if(block_size.isNotNull())
  {
    iblock_size = Rcpp::as<int> (block_size);
  } else {
      iblock_size = 128;
  }
  
  
  
  if( paral.isNull()) {
    bparal = false;
  } else {
    bparal = Rcpp::as<bool> (paral);
  }
  
  
  // Read DelayedmArray's A and b
  if ( A.isS4() == true)    
  {
    mA = read_DelayedArray(A);
  } else {
    try{  
      mA = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(A);
    }
    catch(std::exception &ex) { }
  }
  
  if(B.isNull())
  {
    
    if(btrans == true) {
      
       C = bdtcrossproduct(mA) ;
    
    }else {
    
      C = bdcrossproduct(mA);
    }
    
    
  } else {
    
    // Read DelayedArray's a and b
    if ( as<Rcpp::RObject>(B).isS4() == true) {
      mB = read_DelayedArray(as<Rcpp::RObject>(B));
    } else {
      try{  
        mB = Rcpp::as<Eigen::MatrixXd>(B); 
      }
      catch(std::exception &ex) { }
    }
    
    if( btrans == true )
    {
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mTrans(mB.data(), mB.cols(), mB.rows());
      if(bparal == true) {
        C = Bblock_matrix_mul_parallel(mA, mTrans, iblock_size, threads);
        
      } else if (bparal == false)  {
        C = Bblock_matrix_mul(mA, mTrans, iblock_size);
      }
      
      
    } else {
      
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mTrans(mA.data(), mA.cols(), mA.rows());
      
      if(bparal == true) {
        C = Bblock_matrix_mul_parallel(mTrans, mB, iblock_size, threads);
        
      } else if (bparal == false)  {
        C = Bblock_matrix_mul(mTrans, mB, iblock_size);
      }
    }

    
  }
  
  return(C);
  
  
}







// Weighted product
//' Matrix - Weighted vector Multiplication with numerical or DelayedArray data
//' 
//' This function performs a weighted product of a matrix(X) with a weighted diagonal matrix (w)
//' 
//' @param X numerical or Delayed Array matrix
//' @param w vector with weights
//' @param op string indicating if operation 'XtwX' and 'XwXt' for weighted cross product (Matrix - Vector - Matrix) or 'Xw' and 'wX' for weighted product (Matrix - Vector)
//' @return numerical matrix 
//' @examples
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
//' 
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd bdwproduct(Rcpp::RObject X, Rcpp::RObject w, std::string op)
{
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd W;
  
  
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
  
  if ( w.isS4() == true)   
  {
    W = read_DelayedArray(w);
  }  else {
    W = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(w);
  } 
  
  if(op == "XwXt" || op == "xwxt") {
    return(xwxt(A,W)) ;
  }else if (op == "XtwX" || op == "xtwx") {
    return(xtwx(A,W));
  }else if (op == "Xw" || op == "xw") {
    return(Xw(A,W));
  }else if (op == "wX" || op == "wx" ) {
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
//' @param A numerical or Delayed Array matrix
//' @param w scalar, weight
//' @param op string indicating if operation  "Xw" or "wX"
//' @return numerical matrix 
//' @examples
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
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd bdScalarwproduct(Rcpp::RObject A, double w, std::string op)
{
  
  Eigen::MatrixXd mA;

  // Read DelayedArray's A and b
  if ( A.isS4() == true)    
  {
    mA = read_DelayedArray(A);
  } else {
    try{  
      mA = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(A);
    }
    catch(std::exception &ex) { }
  }
  
  std::transform(op.begin(), op.end(), op.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  
  if (op == "xw" || op == "wx") {
    return(w*mA);
  } else  {
    throw("Invalid option (valid options : xtwx or xwxt)");
  }
}




/***R

*/
