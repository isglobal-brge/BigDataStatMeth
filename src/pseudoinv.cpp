// #include <RcppEigen.h>
#include "include/pseudoinv.h"


Eigen::MatrixXd rcpp_bdpseudoinv(Eigen::MatrixXd* A)
{
  
  char Schar='S';
  char Cchar='C';
  int ione = 1;
  double done = 1.0;
  double dzero = 0.0;

  int m = A->rows();
  int n = A->cols();
  int lda = m;
  int ldu = std::max(1,m);
  int ldvt =  std::max(1, std::min(m, n));
  int k = std::min(m,n);
  int lwork = std::max( 1, 4*std::min(m,n)* std::min(m,n) + 7*std::min(m, n) );
  int info = 0;
  
  Eigen::VectorXd s = Eigen::VectorXd::Zero(k);
  Eigen::VectorXd work = Eigen::VectorXd::Zero(lwork);
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(ldu,k);
  Eigen::MatrixXd vt = Eigen::MatrixXd::Zero(ldvt,n);
  
  
  // dgesvd_( char JOBU, char JOBVT, int M, int N, double* A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT, double WORK, int LWORK, int INFO  );
  dgesvd_( &Schar, &Schar, &m, &n, A->data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork, &info);
  
  Eigen::MatrixXd pinv = Eigen::MatrixXd::Zero(n,m);
  
  //.OpenMP.// omp_set_dynamic(1);
  
#pragma omp parallel for 
  for (int i = 0; i < k; i++){
    double tempS;
    if(s[i] > 1.0e-9)
      tempS = 1.0/s[i];
    else
      tempS = s[i];

    // zscal_ (int* N, double* DA, double* DX, int* INCX )
    dscal_( &m, &tempS, &(u(i*ldu)), &ione );
  }
  
  // dgemm_( char TRANSA, char TRANSB, int M, int N, int K, double ALPHA, double* A, int LDA, double* B, int LDB, double BETA, double* C, int LDC )	
  dgemm_( &Cchar, &Cchar, &n, &m, &k, &done, vt.data(), &k, u.data(), &m, &dzero, pinv.data(), &n );
  
  return(pinv);
  
}


// 
//' Pseudo-Inverse
//' 
//' Compute the pseudo-inverse of a singular matrix
//' 
//' @param X Singular matrix (m x n)
//' @return Pseudo-inverse matrix of A
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdpseudoinv(const Rcpp::RObject & X)
{
  
  auto dmtype = beachmat::find_sexp_type(X);

  try {
    
    Eigen::MatrixXd A;
    
    if ( dmtype == INTSXP || dmtype==REALSXP ) {
      if ( X.isS4() == true){
        A = read_DelayedArray(X);
      }else {
        try{
          A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
        }catch(std::exception &ex) {
          A = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(X);
        }
      }
      
    } else {
      throw std::runtime_error("unacceptable matrix type");
    }
    
    Eigen::MatrixXd pinv = rcpp_bdpseudoinv(&A);
    
    return(Rcpp::wrap(pinv));
    
  } catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }
  
  
}

