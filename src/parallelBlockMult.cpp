#include "include/parallelBlockMult.h"

// Documentació : http://www.netlib.org/lapack/lawnspdf/lawn129.pdf


Eigen::MatrixXd block_matrix_mul(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size)
{

  int M = A.rows();
  int K = A.cols();
  int N = B.cols();
  if( A.cols()==B.rows())
  {
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ; 
    
    int isize = block_size+1;
    int ksize = block_size+1;
    int jsize = block_size+1;
    
    for (int ii = 0; ii < M; ii += block_size)
    {
      if( ii + block_size > M ) isize = M - ii;
      for (int jj = 0; jj < N; jj += block_size)
      {
        if( jj + block_size > N) jsize = N - jj;
        for(int kk = 0; kk < K; kk += block_size)
        {
          if( kk + block_size > K ) ksize = K - kk;
          C.block(ii, jj, std::min(block_size,isize), std::min(block_size,jsize)) = 
            C.block(ii, jj, std::min(block_size,isize), std::min(block_size,jsize)) + 
            (A.block(ii, kk, std::min(block_size,isize), std::min(block_size,ksize)) * 
            B.block(kk, jj, std::min(block_size,ksize), std::min(block_size,jsize)));
          
          if( kk + block_size > K ) ksize = block_size+1;
        }
        if( jj + block_size > N ) jsize = block_size+1;
      }
      if( ii + block_size > M ) isize = block_size+1;
    }
    
    return(C);
    
  }else {
    throw std::range_error("non-conformable arguments");
  }
  
}


Eigen::MatrixXd block_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size)
{
  int ii=0, jj=0, kk=0;
  int chunk = 1;
  int tid;
  int M = A.rows();
  int K = A.cols();
  int N = B.cols();
  
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ;
  
  // Rcpp::Rcout << "Max threads: " << omp_get_max_threads() << "\n";
  
  omp_set_dynamic(1);   // omp_set_dynamic(0); omp_set_num_threads(4);
  
#pragma omp parallel shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
{

  tid = omp_get_thread_num();
  //només per fer proves dels threads i saber que està paralelitzant, sinó no cal tenir-ho descomentat
  /*
  if (tid == 0)   {
    Rcpp::Rcout << "Number of threads: " << omp_get_num_threads() << "\n";
  }
   */
  
#pragma omp for schedule (dynamic)
  
  
  for (int ii = 0; ii < M; ii += block_size)
  {
    // Rcpp::Rcout << "Number of threads: " << omp_get_num_threads() << "\n";
    for (int jj = 0; jj < N; jj += block_size)
    {
      for(int kk = 0; kk < K; kk += block_size)
      {
        C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) = 
          C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) + 
          (A.block(ii, kk, std::min(block_size,M - ii), std::min(block_size,K - kk)) * 
          B.block(kk, jj, std::min(block_size,K - kk), std::min(block_size,N - jj)));
        
      }
    }
  }
}
return(C);
}




//..// Eigen::MatrixXd blockmult(Rcpp::RObject a, Rcpp::RObject b, int block_size, bool paral)
// [[Rcpp::export]]
Eigen::MatrixXd blockmult(Rcpp::RObject a, Rcpp::RObject b, Rcpp::Nullable<int> block_size = R_NilValue, Rcpp::Nullable<bool> paral = R_NilValue)
{
  int iblock_size;
  bool bparal;// = Rcpp::as<double>;
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  Eigen::MatrixXd C;

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
  
  if ( b.isS4() == true)   
  {
    B = read_DelayedArray(b);
  }  else {
    B = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(b);
  } 
  
  if(block_size.isNotNull())
  {
    iblock_size = Rcpp::as<int> (block_size);

    if( iblock_size > std::min(A.rows(),A.cols()) || iblock_size > std::min(B.rows(),B.cols()) )
    {
      iblock_size = std::min(  std::min(A.rows(),A.cols()), std::min(B.rows(),B.cols()));
    }
    
  } else {
    // iblock_size = std::min(  std::min(A.rows(),A.cols()), std::min(B.rows(),B.cols()));
    iblock_size = 1;
  }
  
  if( paral.isNull()) {
    bparal = false;
  } else {
    bparal = Rcpp::as<bool> (paral);
  }
  
  if(bparal == true)
  {
    C = block_matrix_mul_parallel(A, B, iblock_size);
  }else if (bparal == false)
  {
    C = block_matrix_mul(A, B, iblock_size);
  }
  return(C);
}




/*** R

library(microbenchmark)
library(DelayedArray)
# A <- matrix(sample(1:10,100, replace = TRUE), ncol = 10);A

n <- 1000
p <- 1500
A <- matrix(rnorm(n*p), nrow=n, ncol=p)
  
  
n <- 1500
p <- 1000
B <- matrix(rnorm(n*p), nrow=n, ncol=p)

AD <- DelayedArray(A)
BD <- DelayedArray(B)


results <- microbenchmark( CP2 <- blockmult(A,B,128, FALSE),
                           CPP2 <- blockmult(A,B,128, TRUE),
                           CP2D <- blockmult(AD,BD,128, FALSE),
                           CPP2D <- blockmult(AD,BD,128, TRUE),
                           dgemm <- matv_gemm(A,B),
                           r <- A%*%B,
                           times = 2L)
  
print(summary(results)[, c(1:7)],digits=3)

stopifnot(all.equal(CP2,CP2D),
          all.equal(CPP2,CPP2D),
          all.equal(CP2,CPP2),
          all.equal(CP2,dgemm),
          all.equal(CP2,r))


  D <- block_matrix_mul_parallel(A,B,128)
  CP2 <- blockmult(A,B,128)
  
  results <- microbenchmark(CPP1 <- blockmult(A,B,64, TRUE),
                            CPP2 <- blockmult(A,B,128, TRUE),
                            CPP3 <- blockmult(A,B,256, TRUE),
                            CPP4 <- blockmult(A,B,512, TRUE),
                            CPP5 <- blockmult(A,B,1024, TRUE),
                            CP1 <- blockmult(A,B,64, FALSE),
                            CP2 <- blockmult(A,B,128, FALSE),
                            CP3 <- blockmult(A,B,256, FALSE),
                            CP4 <- blockmult(A,B,512, FALSE),
                            CP5 <- blockmult(A,B,1024, FALSE),
#C <- A%*%B,
                                                                  times = 2L)  # Proves multiplicacions x blocs
  
  print(summary(results)[, c(1:7)],digits=3)
  
  
  
  stopifnot(all.equal(block_matrix_mul_parallel(A,B,64),A%*%B), 
            all.equal(blockmult(A,B,128), A%*%B))
  
  CP <- particiona_matrius(A,B,2)
  CP
  C <- A%*%B
  stopifnot(all.equal(particiona_matrius(A,B,2),A%*%B ))
  
  
  
  
  
n <- 500
p <- 750
A <- matrix(rnorm(n*p), nrow=n, ncol=p)
B <- matrix(rnorm(p*p), nrow=p, ncol=p)

AD <- DelayedArray(A)
BD <- DelayedArray(B)

B <- t(A)

CP2 <- blockmult(AD,BD)

CP2[1:10,1:10]
CPD <- A%*%B
CPD[1:10,1:10]
all.equal(CPD,CP2)

CPP2 <- blockmult(A,B,128, TRUE)
stopifnot( all.equal(blockmult(A,B,128, FALSE),A%*%B),
           all.equal(blockmult(A,B,128, TRUE),A%*%B))
  
  E <- rcppparallel_blockmult(A,B,32)
  
  D <- block_matrix_mul_parallel(A,B,2)
  stopifnot(all.equal(block_matrix_mul_parallel(A,B,2),A%*%B))
  
  results <- microbenchmark(CP <- particiona_matrius(A,B),
                            C <- A%*%B, 
                               times = 5L)  # Prpves ,iñto`ñocacopms x blocs
  
  print(summary(results)[, c(1:7)],digits=3)
  
  
  
  
  A <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),byrow = TRUE, nrow = 4)
  B <- t(A)
  
  particiona_matrius(A,B)
  
  C <- A%*%B
  
# Partim matrius en 4 blocks :
#     1a Fila : A1  A2
#     2a Fila : A3  A4
  A00 <- A[1:2,1:2]
A01 <- A[1:2,3:4]
A10 <- A[3:4,1:2]
A11 <- A[3:4,3:4]

B00 <- B[1:2,1:2]
B01 <- B[1:2,3:4]
B10 <- B[3:4,1:2]
B11 <- B[3:4,3:4]

C00 <- A00%*%B00 + A01%*%B10;C00
  C01 <- A00%*%B01 + A01%*%B11;C01
    C10 <- A10%*%B00 + A11%*%B10;C10
      C11 <- A10%*%B01 + A11%*%B11;C11
        
        CT <- rbind(cbind(C00,C01),cbind(C10,C11))
        CT;CT
          
          */
