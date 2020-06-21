#include "include/qrdecomp.h"


using namespace std;


strQR rcpp_bdQR( Eigen::MatrixXd & A, bool bthin)
{
  
  int m = A.rows(), n = A.cols();
  int irank;
  strQR vQR;
  
  Eigen::MatrixXd R;
  Eigen::MatrixXd Q;
  
  Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  
  qr.compute(A);
  irank = lu_decomp.rank();
  
  if (irank == m + 1 || irank == n + 1 )
  {
    vQR.R = qr.matrixQR().template triangularView<Eigen::Upper>();
  } else {
    vQR.R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>(); 
  }

  if (bthin == false)
  {
    vQR.Q =  qr.householderQ();       // Full decomposition
  } else {
    
    vQR.Q = Eigen::MatrixXd::Identity(m,n);
    vQR.Q = qr.householderQ() * vQR.Q;    // Thin decomposition
  }
  
  return(vQR);
  
}




//' QR Decomposition 
//' 
//' This function compute QR decomposition (also called a QR factorization) 
//' of a matrix \code{A} into a product \code{A = QR} of an 
//' orthogonal matrix Q and an upper triangular matrix R.
//' 
//' @param X a real square matrix 
//' @param boolean thin, if thin = true returns Q thin  decomposition else returns Q full decomposition, default thin = false
//' @return List with orthogonal matrix \code{Q}  and upper triangular matrix \code{R}
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdQR( const Rcpp::RObject & X, Rcpp::Nullable<bool> thin = R_NilValue)
{
  
  auto dmtype = beachmat::find_sexp_type(X);
  Eigen::MatrixXd A;
  bool bthin;
  strQR decQR;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( X.isS4() == true){
      A = read_DelayedArray(X);
    }else {
      try{
        A = Rcpp::as<Eigen::MatrixXd >(X);
      }catch(std::exception &ex) {
        A = Rcpp::as<Eigen::VectorXd >(X);
      }
    }
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  
  if( thin.isNull()) {
    bthin = false;
  } else {
    bthin = Rcpp::as<bool> (thin);
  }
  
  decQR = rcpp_bdQR(A, bthin);
  
  return Rcpp::List::create(Rcpp::Named("Q") = decQR.Q,
                            Rcpp::Named("R") = decQR.R
  );
}




/*
//' QR Decomposition 
//' 
//' This function compute QR decomposition (also called a QR factorization) 
//' of a matrix \code{A} into a product \code{A = QR} of an 
//' orthogonal matrix Q and an upper triangular matrix R.
//' 
//' @param X a real square matrix 
//' @param boolean thin, if thin = true returns Q thin  decomposition else returns Q full decomposition, default thin = false
//' @return List with orthogonal matrix \code{Q}  and upper triangular matrix \code{R}
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdQR( const Rcpp::RObject & X, Rcpp::Nullable<bool> thin = R_NilValue)
{
  
  auto dmtype = beachmat::find_sexp_type(X);
  Eigen::MatrixXd A;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( X.isS4() == true){
      A = read_DelayedArray(X);
    }else {
      try{
        A = Rcpp::as<Eigen::MatrixXd >(X);
      }catch(std::exception &ex) {
        A = Rcpp::as<Eigen::VectorXd >(X);
      }
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  int m = A.rows(), n = A.cols();
  int irank;
  bool bthin;
  
  Eigen::MatrixXd R;
  Eigen::MatrixXd Q;
  
  Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  
  qr.compute(A);
  irank = lu_decomp.rank();
  
  if (irank == m + 1 || irank == n + 1 )
  {
    R = qr.matrixQR().template triangularView<Eigen::Upper>();
  } else {
    R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>(); 
  }
  
  if( thin.isNull()) {
    bthin = false;
  } else {
    bthin = Rcpp::as<bool> (thin);
  }
  
  if (bthin == false)
  {
    Q =  qr.householderQ();       // Full decomposition
  } else {
    
    Q = Eigen::MatrixXd::Identity(m,n);
    Q = qr.householderQ() * Q;    // Thin decomposition
  }
  
  return Rcpp::List::create(Rcpp::Named("Q") = Q,
                            Rcpp::Named("R") = R
  );
}

*/

// [[Rcpp::export]]
Rcpp::RObject review_decomposition(Eigen::MatrixXd R, int n)
{
  
  // if R is smaller than nxn
  while( R.cols() < n )
  {
    R.conservativeResize(R.rows()+1, R.cols()+1);
    R.col(R.cols()-1) = Eigen::VectorXd::Zero(R.rows());
    R.row(R.rows()-1) = Eigen::VectorXd::Zero(R.cols());
    R(R.rows()-1, R.cols()-1) = 0.0000000000000001;
  }
  
  
  if( (R.colwise().sum().array() == 0.0).any())
  {
    Eigen::VectorXd colsum = R.colwise().sum();
    for (int i = 0; i < n; i++){
      if( colsum(i) == 0 )   R(i, i) = 0.0000000000000001;
    }
  }
  
  return(Rcpp::wrap(R));
  
}





//' Solves matrix equations : A*X = B
//' 
//' 
//' 
//' @param R numerical or Delayed Array matrix. 
//' @param Z numerical or Delayed Array matrix.
//' @return X numerical matrix. 
//' @examples
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bddtrsm(Rcpp::RObject R, Rcpp::RObject Z, Rcpp::Nullable<int> threads = R_NilValue) 
{
  char Nchar='N';
  char Uchar='U';
  char Lchar='L';
  double done = 1.0;
  
  
  Eigen::MatrixXd A = Rcpp::as<Eigen::MatrixXd> (R);
  Eigen::FullPivLU<Eigen::MatrixXd>matinv(A);
  Eigen::MatrixXd B = Rcpp::as<Eigen::MatrixXd> (Z);
  
  int m = B.rows(), n = B.cols();
  int lda = std::max( 1, m );
  int ldb = std::max( 1, m );
  
  if(matinv.isInvertible()==true)
  {
    
    // dtrsm( char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double A, int LDA, double B, int LDB )
    dtrsm_(&Lchar, &Uchar, &Nchar, &Nchar, &m, &n, &done, A.data(), &lda, B.data(), &ldb );
    return(Rcpp::wrap(B));
    
  } else {

    int block_size = 128;
    if( block_size > std::min(A.rows(),A.cols()) || block_size > std::min(B.rows(),B.cols()) )  {
      block_size = std::min(  std::min(A.rows(),A.cols()), std::min(B.rows(),B.cols()));
    }
    // Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd>( bdpseudoinv(A) )* B;
    //.. MODIFICAT 06/2020 ..// Eigen::MatrixXd X = block_matrix_mul_parallel( rcpp_bdpseudoinv(A), B, block_size, threads );
    Eigen::MatrixXd X = Bblock_matrix_mul_parallel( rcpp_bdpseudoinv(A), B, block_size, threads );
    return(Rcpp::wrap(X));
  }
  
  return(0);
  
}






/***R

a <- matrix(c(1,-1,4,1,4,-2,1,4,2,1,-1,0), byrow = TRUE, ncol = 3)
bdQR(a)
ad <- DelayedArray(a)
bdQR(ad)

b <- matrix(c(-1,2,2), byrow = TRUE, ncol = 3)
solve(b)
bdpseudoinv(b)
bd <- DelayedArray(b)
bdpseudoinv(bd)


*/