#include "include/qrdecomp.h"


using namespace std;

// 
// // Housholder's QR decomposition for OMP
// strQR rcpp_bdQR_parallel( Eigen::MatrixXd A, Rcpp::Nullable<int> threads  = R_NilValue )
// {
//   
//    strQR vQR;
//   int ithreads;
//    
//    try {
//      
//      A.transposeInPlace();
//      
//      int irows = A.rows(),
//          icols = A.cols();
//      
//     if(threads.isNotNull())
//     {
//       if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
//         ithreads = Rcpp::as<int> (threads);
//       else
//         ithreads = std::thread::hardware_concurrency()-1;
//     }
//     else    ithreads = std::thread::hardware_concurrency()-1; 
// 
//         omp_set_dynamic(0);  
//         omp_set_num_threads(ithreads);
// 
//     
//     bool verbose = true, definition = true;
//     int i,j;
//     double *vec;
//     
//     // matrixQ = create_matrix(irows,irows, verbose);
//     Eigen::MatrixXd matrixQ = Eigen::MatrixXd::Identity(irows,icols);
//     Eigen::MatrixXd tmpmat = Eigen::MatrixXd::Zero(irows,icols);
// 
//     int k,l,m;
//     for(i = 0; i < irows; i++){
//       
//       double x = 0;
//       vec = new double[irows-i];
//       for(j = i; j < irows; j++){
//         vec[j-i] = -A(j,i);
//         x += vec[j-i]*vec[j-i];
//       }
//       x = sqrt(x);
//       
//       if(vec[0] > 0) x = -x;
//       vec[0] = vec[0] + x;
//       x = 0;
//       for(j = 0; j < irows-i; j++){
//         x += vec[j]*vec[j];
//       }
//       x = sqrt(x);
//       
//       if(x > 0){
//         //normalizovat vec
//         for(j = 0; j < irows-i; j++){
//           vec[j] /= x;
//         }
//         
//         Eigen::MatrixXd pp = Eigen::MatrixXd::Zero(irows-i, irows-i);
//         for(k = 0; k < irows-i; k++){
// #pragma omp parallel for
//           for(l = 0; l < irows-i; l++){
//             if(k == l) pp(k,k) = 1 - 2*vec[k]*vec[l];
//             else pp(k,l) = -2*vec[k]*vec[l];
//           }
//         }
//         
// 
//         //R
//         double tm;
//         for(k = i; k < irows; k++){
// #pragma omp parallel for private(tm,m)
//           for(l = i; l < irows; l++){
//             tm = 0;
//             for(m = i; m < irows; m++){
//               tm += pp(k-i,m-i)*A(m,l);
//             }
//             tmpmat(k,l) = tm;
//           }
//         }
//         for(k = i; k < irows; k++){
// #pragma omp parallel for
//           for(l = i; l < irows; l++){
//             A(k,l) = tmpmat(k,l);
//             
//           }
//         }
//         
//         //Q
//         for(k = 0; k < irows; k++){
// #pragma omp parallel for private(tm,m)
//           for(l = i; l < irows; l++){
//             tm = 0;
//             for(m = i; m < irows; m++){
//               tm += matrixQ(k,m)*pp(m-i,l-i);
//             }
//             tmpmat(k,l) = tm;
//           }
//         }
//         for(k = 0; k < irows; k++){
// #pragma omp parallel for
//           for(l = i; l < irows; l++){
//             matrixQ(k,l) = tmpmat(k,l);
//           }
//         }
//       }
//     }
//       
//       vQR.Q = matrixQ;
//       vQR.R = A;
//     
//     
//   } catch(std::exception &ex) {
//     Rcpp::Rcout<< ex.what();
//     return vQR;
//   }
//   
//   return(vQR);
//   
//   }



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
//' @param thin boolean thin, if thin = true returns Q thin  decomposition else returns Q full decomposition, default thin = false
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



// //' QR Decomposition 
// //' 
// //' This function compute QR decomposition computes a QR
// //' factorization with column pivoting of a matrix A:  A*P = Q*R
// //' 
// //' This QR decomposition is like R native decomposition obtained with qr() function
// //' 
// //' @param X a real square matrix 
// //' @return List with  : a matrix with the same dimensions as x. The upper triangle contains the \(\bold{R}\) of the decomposition and the lower triangle contains information on the \(\bold{Q}\) of the decomposition (stored in compact form)
// //' @export
// // [[Rcpp::export]]
// Rcpp::RObject bdQR_compact( const Rcpp::RObject & A)
// {
//   
//   auto dmtype = beachmat::find_sexp_type(A);
//   Eigen::MatrixXd a;
//   int info = 0;
//   try {
//     
//     if ( dmtype == INTSXP || dmtype==REALSXP ) {
//       if ( A.isS4() == true){
//         a = read_DelayedArray(A);
//       }else {
//         try{
//           a = Rcpp::as<Eigen::Map<Eigen::MatrixXd > >(A);
//         }catch(std::exception &ex) {
//           a = Rcpp::as<Eigen::Map<Eigen::MatrixXd > >(A);
//         }
//       }
//     } else {
//       throw std::runtime_error("unacceptable matrix type");
//     }
//     
//     // Declare matrix variables
//     int m = a.rows();
//     int n = a.cols();
//     int minmn = std::min(m,n);
//     int lda = std::max( 1, m );
//     int lwork = 3*n+1;
//     double work[lwork];
//     Eigen::VectorXi jpvt = Eigen::VectorXi::Zero(n);
//     Eigen::VectorXd tau = Eigen::VectorXd::Zero(minmn);
//     
//     
//     // dgeqp3_( int* M, int* N, double* A, int* LDA, int* JPVT, double* TAU, double* WORK, int* LWORK, int* INFO);
//     dgeqp3_( &m, &n, a.data(), &lda, jpvt.data(), tau.data(), work, &lwork, &info);
//     
//     if(info<0){
//       throw std::runtime_error("the i-th argument had an illegal value.");
//     } else {
//       return Rcpp::List::create(Rcpp::Named("QR") = wrap(a),
//                                 Rcpp::Named("pivot") = wrap(jpvt),
//                                 Rcpp::Named("qraux") = wrap(tau)
//                                 );
//     }
//     
//   }catch(std::exception &ex) {
//     Rcpp::Rcout<< ex.what();
//     return wrap(-1);
//   }
//   
// }

// 
// //' QR Decomposition 
// //' 
// //' This function compute QR decomposition (also called a QR factorization) 
// //' of a matrix \code{A} into a product \code{A = QR} of an 
// //' orthogonal matrix Q and an upper triangular matrix R.
// //' 
// //' @param X a real square matrix 
// //' @param threads (optional) only if paral = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
// //' @return List with orthogonal matrix \code{Q}  and upper triangular matrix \code{R}
// //' @export
// // [[Rcpp::export]]
// Rcpp::RObject bdQR_parallel( const Rcpp::RObject & X, Rcpp::Nullable<int> threads = R_NilValue)
// {
//   
//   auto dmtype = beachmat::find_sexp_type(X);
//   Eigen::MatrixXd A;
//   strQR decQR;
//   
//   if ( dmtype == INTSXP || dmtype==REALSXP ) {
//     if ( X.isS4() == true){
//       A = read_DelayedArray(X);
//     }else {
//       try{
//         A = Rcpp::as<Eigen::MatrixXd >(X);
//       }catch(std::exception &ex) {
//         A = Rcpp::as<Eigen::VectorXd >(X);
//       }
//     }
//   } else {
//     throw std::runtime_error("unacceptable matrix type");
//   }
//   
//   decQR = rcpp_bdQR_parallel(A, threads);
//   
//   return Rcpp::List::create(Rcpp::Named("Q") = decQR.Q,
//                             Rcpp::Named("R") = decQR.R
//                             );
// }
// 




//' Solves matrix equations : A*X = B
//' 
//' 
//' 
//' @param R numerical or Delayed Array matrix. 
//' @param Z numerical or Delayed Array matrix.
//' @param threads integer with number of threads to use with parallelized execution
//' @return X numerical matrix. 
//' @examples
//'  a <- "Unused function"
//' @export
// [[Rcpp::export]]
Rcpp::RObject bddtrsm(Rcpp::RObject R, Rcpp::RObject Z, Rcpp::Nullable<int> threads = R_NilValue) 
{
  
  ////// IMPORTANT !!!! :  op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
  //                              where alpha is a scalar, X and B are m by n matrices
  
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
    //.. Modified 06/2020 ..// Eigen::MatrixXd X = block_matrix_mul_parallel( rcpp_bdpseudoinv(A), B, block_size, threads );
    Eigen::MatrixXd X = Bblock_matrix_mul_parallel( rcpp_bdpseudoinv(&A), B, block_size, threads );
    return(Rcpp::wrap(X));
  }
  
  return(0);
  
}






/***R

library(DelayedArray)
library(microbenchmark)
library(BigDataStatMeth)
library(devtools)

# reload(pkgload::inst("BigDataStatMeth"))

n <- 500
m <- 500

# R Object

set.seed(447)
a <- matrix(runif(n*m), nrow = n, ncol = m)
ad <- DelayedArray(a)

# a <- matrix(c(4,3,2,1,3,4,3,2,2,3,4,3,1,2,3,4), byrow = TRUE, nrow = 4)

resR <- qrR(a)
resP <- bdQR_parallel(a,2)

resR$Q
resP$Q

resR$R
resP$R



reload(pkgload::inst("BigDataStatMeth"))
sortida <- bdQR_parallel(a,2)

qrR <- function(A)
{
  res <- qr(A)
  Q <- qr.Q(res)
  R <- qr.R(res)
  return(list("Q" = Q,
              "R" = R))
}

reload(pkgload::inst("BigDataStatMeth"))
aix <- qrR(a)
aix$Q[1:5,1:5]
aix$R[1:5,1:5]

aix2 <- bdQR_parallel(a,2)
aix2$Q[1:5,1:5]
aix2$R[1:5,1:5]

reload(pkgload::inst("BigDataStatMeth"))
res <- microbenchmark( bdQR(a),
                       bdQR(ad),
                       qrR(a),
                       bdQR_compact(a),
                       bdQR_compact(ad),
                       bdQR_parallel(a,2),
                       bdQR_parallel(a,3),
                       bdQR_parallel(a,4),
                       times = 1)

res






bddtrsm(b)

n <- 700
m <- 500

# R Object

set.seed(222)
mat <- matrix(runif(n*m), nrow = n, ncol = m)
# matD<- DelayedArray(mat)

res <- microbenchmark::microbenchmark( solve(mat),
                                       bdpseudoinv(mat),
                                       # bdpseudoinv(matD),
                                       times = 1)

res

solve(a)


# TEST PSEUDO-INVERSE 

library(corpcor)

A <- matrix(c(5,-1,-1,-1,-1,5,-1,-1,-1,-1,5,-1,-1,-1,-1,5),byrow = TRUE, nrow = 3)

t <- matrix(c(7,2,3,4,5,3,4,5,3,2,6,3,4,7,8,9,3,5),byrow = TRUE, nrow = 3)

solve(A)
round(bdpseudoinv(t),2)

microbenchmark( bdpseudoinv(t),
                pseudoinverse(t) )

all.equal(zapsmall(bdpseudoinv(t)), zapsmall(pseudoinverse(t)))

zapsmall(bdpseudoinv(t))
zapsmall(pseudoinverse(t))



bddtrsm(A) 




# MATRIX EQUATION SOLVER

A <- matrix(c(7, -3, 6, 5, -2, 2, 2, -3, 8), byrow = TRUE, nrow = 3)
B <- matrix(c(5,11,10), nrow = 3)

A <- matrix(c(3, -6, 5, 2), byrow = TRUE, nrow = 2)
B <- matrix(c(7, 4, -5, 8), byrow = TRUE, nrow = 2)

A <- matrix(c(3, -1, 2, 1), byrow = TRUE, nrow = 2)
B <- matrix(c(5, 5), byrow = TRUE, nrow = 2)

A;B
X <- bddtrsm(A, B) ; X


solve(A,B)

bdd
  
  */