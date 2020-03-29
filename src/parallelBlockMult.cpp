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


Eigen::MatrixXd block_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, 
                                          int block_size, Rcpp::Nullable<int> threads  = R_NilValue)
{
  int ii=0, jj=0, kk=0;
  int chunk = 1, tid;
  unsigned int ithreads;
  int M = A.rows();
  int K = A.cols();
  int N = B.cols();
  
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ;
  if(block_size > std::min( N, std::min(M,K)) )
    block_size = std::min( N, std::min(M,K)); 
  
  if(threads.isNotNull())    ithreads = Rcpp::as<int> (threads);
  else    ithreads = std::thread::hardware_concurrency(); //omp_get_max_threads();

  omp_set_dynamic(1);   // omp_set_dynamic(0); omp_set_num_threads(4);
  omp_set_num_threads(ithreads);
  
#pragma omp parallel shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
{

  // tid = omp_get_thread_num();
  //només per fer proves dels threads i saber que està paralelitzant, sinó no cal tenir-ho descomentat
  // if (tid == 0)   {
  //   Rcpp::Rcout << "Number of threads: " << omp_get_num_threads() << "\n";
  // }
   
#pragma omp for schedule (static) 
  
  
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





//' Block matrix multiplication with Delayed Array Object
//' 
//' This function performs a block matrix-matrix multiplication with numeric matrix or Delayed Arrays
//' 
//' @param a a double matrix.
//' @param b a double matrix.
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @return numerical matrix
//' @examples
//' # with numeric matrix
//' m <- 500
//' k <- 1500
//' n <- 400
//' A <- matrix(rnorm(n*p), nrow=n, ncol=k)
//' B <- matrix(rnorm(n*p), nrow=k, ncol=n)
//' 
//' blockmult(A,B,128, TRUE)
//' 
//' # with Delaeyd Array
//' AD <- DelayedArray(A)
//' BD <- DelayedArray(B)
//' 
//' blockmult(AD,BD,128, TRUE)
//' 
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd blockmult(Rcpp::RObject a, Rcpp::RObject b, 
                          Rcpp::Nullable<int> block_size = R_NilValue, 
                          Rcpp::Nullable<bool> paral = R_NilValue,
                          Rcpp::Nullable<int> threads = R_NilValue)
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
      if ( TYPEOF(a) == INTSXP ) {
        A = Rcpp::as<Eigen::MatrixXi>(a).cast<double>()  ;
      } else{
        A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
      }
    }
    catch(std::exception &ex) { }
  }
  
  if ( b.isS4() == true)   
  {
    B = read_DelayedArray(b);
  }  else {
    
    if ( TYPEOF(b) == INTSXP ) {
      B = Rcpp::as<Eigen::MatrixXi>(b).cast<double>()  ;
    } else{
      B = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(b);
    }

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
    iblock_size = 128;
  }
  
  if( paral.isNull()) {
    bparal = false;
  } else {
    bparal = Rcpp::as<bool> (paral);
  }
  
  if(bparal == true)
  {
    C = block_matrix_mul_parallel(A, B, iblock_size, threads);
  }else if (bparal == false)
  {
    C = block_matrix_mul(A, B, iblock_size);
  }
  return(C);
}




/*** R

library(microbenchmark)
library(DelayedArray)
library(BigDataStatMeth)
# A <- matrix(sample(1:10,100, replace = TRUE), ncol = 10);A

n <- 1000
p <- 400
A <- matrix(rnorm(n*p), nrow=n, ncol=p)
  
  
n <- 400
p <- 1000
B <-  matrix(sample(1:10, 40, replace = TRUE), nrow = n, ncol = p)

CPP2 <- blockmult(A,B,128, TRUE)

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
  


*/
