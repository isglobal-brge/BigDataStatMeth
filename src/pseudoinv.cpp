// #include <RcppEigen.h>
#include "include/pseudoinv.h"

/*
using namespace std;

// dgemm_ is a symbol in the LAPACK-BLAS library files 
//    DGEMM  performs one of the matrix-matrix operations : C := alpha*op( A )*op( B ) + beta*C,
extern "C" {
  extern void dgemm_( char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int* );
}


// dgesvd_ is a symbol in the LAPACK-BLAS Level 3 
//    DGESVD computes the singular value decomposition (SVD) of a real M-by-N matrix A, 
//       optionally computing the left and/or right singular vectors
extern "C" {
  extern void dgesvd_( char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
}
*/


Eigen::MatrixXd rcpp_bdpseudoinv(Eigen::MatrixXd A)
{
  
  char Schar='S';
  char Cchar='C';
  //char Ochar='O';
  int ione = 1;
  double done = 1.0;
  double dzero = 0.0;
  //int izero = 0;
  
  int m = A.rows();
  int n = A.cols();
  int lda = m;
  int ldu = std::max(1,m);
  int ldvt =  std::max(1, std::min(m, n));
  int k = std::min(m,n);
  //int l = std::max(m,n);
  int lwork = std::max( 1, 4*std::min(m,n)* std::min(m,n) + 7*std::min(m, n) );
  int info = 0;
  
  Eigen::VectorXd s = Eigen::VectorXd::Zero(k);
  Eigen::VectorXd work = Eigen::VectorXd::Zero(lwork);
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(ldu,k);
  Eigen::MatrixXd vt = Eigen::MatrixXd::Zero(ldvt,n);
  
  
  // dgesvd_( char JOBU, char JOBVT, int M, int N, double* A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT, double WORK, int LWORK, int INFO  );
  dgesvd_( &Schar, &Schar, &m, &n, A.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &k, work.data(), &lwork, &info);
  
  double tempS;
  Eigen::MatrixXd pinv = Eigen::MatrixXd::Zero(m,n);
  
  for (int i = 0; i < k; i++){
    tempS = 1 / s[i];
    // zscal_ (int* N, double* DA, double* DX, int* INCX )
    dscal_( &m, &tempS, &(u(i*ldu)), &ione );
  }
  
  // dgemm_( char TRANSA, char TRANSB, int M, int N, int K, double ALPHA, double* A, int LDA, double* B, int LDB, double BETA, double* C, int LDC )	
  dgemm_( &Cchar, &Cchar, &n, &m, &k, &done, vt.data(), &k, u.data(), &m, &dzero, pinv.data(), &n );
  
  return(pinv);
  
}


/*
// 
//' Pseudo-Inverse
//' 
//' Compute the pseudo-inverse of a singular matrix
//' 
//' @param Singular matrix (m x n)
//' @return Pseudo-inverse matrix of A
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdpseudoinv(const Rcpp::RObject & X)
{
  
  auto dmtype = beachmat::find_sexp_type(X);
  char Schar='S';
  char Cchar='C';
  char Ochar='O';
  int ione = 1;
  double done = 1.0;
  double dzero = 0.0;
  int izero = 0;
  
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
  
  int m = A.rows();
  int n = A.cols();
  int lda = m;
  int ldu = std::max(1,m);
  int ldvt =  std::max(1, std::min(m, n));
  int k = std::min(m,n);
  int l = std::max(m,n);
  int lwork = std::max( 1, 4*std::min(m,n)* std::min(m,n) + 7*std::min(m, n) );
  int info = 0;
  
  Eigen::VectorXd s = Eigen::VectorXd::Zero(k);
  Eigen::VectorXd work = Eigen::VectorXd::Zero(lwork);
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(ldu,k);
  Eigen::MatrixXd vt = Eigen::MatrixXd::Zero(ldvt,n);
  
  
  // dgesvd_( char JOBU, char JOBVT, int M, int N, double* A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT, double WORK, int LWORK, int INFO  );
  dgesvd_( &Schar, &Schar, &m, &n, A.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &k, work.data(), &lwork, &info);
  
  double tempS;
  Eigen::MatrixXd pinv = Eigen::MatrixXd::Zero(m,n);
  
  for (int i = 0; i < k; i++){
    tempS = 1 / s[i];
    // zscal_ (int* N, double* DA, double* DX, int* INCX )
    dscal_( &m, &tempS, &(u(i*ldu)), &ione );
  }
  
  // dgemm_( char TRANSA, char TRANSB, int M, int N, int K, double ALPHA, double* A, int LDA, double* B, int LDB, double BETA, double* C, int LDC )	
  dgemm_( &Cchar, &Cchar, &n, &m, &k, &done, vt.data(), &k, u.data(), &m, &dzero, pinv.data(), &n );
  
  
  return(Rcpp::wrap(pinv));
  
}
*/


// 
//' Pseudo-Inverse
//' 
//' Compute the pseudo-inverse of a singular matrix
//' 
//' @param Singular matrix (m x n)
//' @return Pseudo-inverse matrix of A
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdpseudoinv(const Rcpp::RObject & X)
{
  
  auto dmtype = beachmat::find_sexp_type(X);
  char Schar='S';
  char Cchar='C';
  char Ochar='O';
  int ione = 1;
  double done = 1.0;
  double dzero = 0.0;
  int izero = 0;
  
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
  
  Eigen::MatrixXd pinv = rcpp_bdpseudoinv(A);
  
  return(Rcpp::wrap(pinv));
  
}

