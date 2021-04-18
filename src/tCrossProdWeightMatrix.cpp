#include "include/tCrossProdWeightMatrix.h"

using namespace std;
using namespace Rcpp;

// In-memory execution - Serial version
Eigen::MatrixXd Bblock_weighted_crossprod(const Eigen::MatrixXd& A, Eigen::MatrixXd& B, int block_size)
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
    
    B.resize(0,0);
    B.resize(C.rows(),C.cols());
    B = C;
    
    C.resize(B.rows(),A.rows());
    C = Eigen::MatrixXd::Zero(B.rows(),A.rows()) ;
    
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
            (B.block(ii, kk, std::min(block_size,isize), std::min(block_size,ksize)) * 
            A.block(jj, kk, std::min(block_size,jsize), std::min(block_size,ksize)).transpose() );
          
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




// In-memory execution - Parallel version
Eigen::MatrixXd Bblock_weighted_crossprod_parallel(const Eigen::MatrixXd& A, Eigen::MatrixXd& B, 
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
  
  if(threads.isNotNull()) 
  {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
      ithreads = Rcpp::as<int> (threads);
    else 
      ithreads = std::thread::hardware_concurrency()/2;
  }
  else    ithreads = std::thread::hardware_concurrency()/2; //omp_get_max_threads();

  omp_set_dynamic(0);   // omp_set_dynamic(0); omp_set_num_threads(4);
  omp_set_num_threads(ithreads);
  
#pragma omp parallel shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
{

  tid = omp_get_thread_num();
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


  B.resize(0,0);
  B.resize(C.rows(),C.cols());
  B = C;

  C.resize(B.rows(),A.rows());
  C = Eigen::MatrixXd::Zero(B.rows(),A.rows()) ;
  
  omp_set_dynamic(0);   
  omp_set_num_threads(ithreads);

#pragma omp parallel shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
{
  
  tid = omp_get_thread_num();
  
#pragma omp for schedule (static) 
  
  
  for (int ii = 0; ii < M; ii += block_size)
  {
    for (int jj = 0; jj < N; jj += block_size)
    {
      for(int kk = 0; kk < K; kk += block_size)
      {
        C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) = 
          C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) + 
          (B.block(ii, kk, std::min(block_size,M - ii), std::min(block_size,K - kk)) * 
          // A.block(kk, jj, std::min(block_size,K - kk), std::min(block_size,N - jj)));
          A.block(jj, kk, std::min(block_size,N - jj), std::min(block_size,K - kk)).transpose() );
      }
    }
  }
}


return(C);
}



//' Block matrix multiplication with Delayed Array Object
//' 
//' This function performs a Crossproduct with weigths matrix A%*%W%*%t(A) multiplication with numeric matrix or Delayed Arrays
//' 
//' @param a a double matrix.
//' @param w a Weighted matrix
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @return Matrix with crossproduct : 
//' @examples
//' 
//' library(DelayedArray)
//' 
//' # with numeric matrix
//' m <- 500
//' k <- 1500
//' n <- 400
//' A <- matrix(rnorm(n*k), nrow=n, ncol=k)
//' B <- matrix(rnorm(n*k), nrow=k, ncol=n)
//' 
//' blockmult(A,B,128, TRUE)
//' 
//' # with Delaeyd Array
//' AD <- DelayedArray(A)
//' BD <- DelayedArray(B)
//' 
//' blockmult(AD,BD,128, TRUE)
//' @export
// [[Rcpp::export]]
Rcpp::List tCrossprod_Weighted(Rcpp::RObject a, Rcpp::RObject w, 
                               Rcpp::Nullable<int> block_size = R_NilValue, 
                               Rcpp::Nullable<bool> paral = R_NilValue,
                               Rcpp::Nullable<int> threads = R_NilValue )
{
  
  int iblock_size;
  bool bparal; 

  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  Eigen::MatrixXd C;
  
  IntegerVector dsizeA, dsizeB;
  
  // Rcpp::Rcout<<"\n Tipus de dades :  "<<TYPEOF(a)<<"\n";
  // Rcpp::Rcout<<"\n Clase objecte :  "<<  as<std::string>(a.slot("class"))  <<"\n";
  
  
  try{
    
    // Get matrix sizes
    if ( a.isS4() == true )    
    {
      dsizeA = get_DelayedArray_size(a);
    } else { 
      try{  
        if ( TYPEOF(a) == INTSXP ) {
          dsizeA[0] = Rcpp::as<IntegerMatrix>(a).nrow();
          dsizeA[1] = Rcpp::as<IntegerMatrix>(a).ncol();
        }else{
          dsizeA[0] = Rcpp::as<NumericMatrix>(a).nrow();
          dsizeA[1] = Rcpp::as<NumericMatrix>(a).ncol();
        }
      }catch(std::exception &ex) { }
    }
    
    if ( w.isS4() == true)    
    {
      dsizeB = get_DelayedArray_size(w);
    } else { 
      try{  
        if ( TYPEOF(w) == INTSXP ) {
          dsizeB[0] = Rcpp::as<IntegerMatrix>(w).nrow();
          dsizeB[1] = Rcpp::as<IntegerMatrix>(w).ncol();
        }else{
          dsizeB[0] = Rcpp::as<NumericMatrix>(w).nrow();
          dsizeB[1] = Rcpp::as<NumericMatrix>(w).ncol();
        }
      }catch(std::exception &ex) { }
    }
    
    if(block_size.isNotNull())
    {
      iblock_size = Rcpp::as<int> (block_size);
  
    } else {
      iblock_size = std::min(  std::min(dsizeA[0],dsizeA[1]),  std::min(dsizeB[0],dsizeB[1]));
      if (iblock_size>128)
        iblock_size = 128;
    }
    
    if( paral.isNull()) {
      bparal = false;
    } else {
      bparal = Rcpp::as<bool> (paral);
    }
    
    
    
    // if number of elemenents < bigmat in all matrix work in memory else work with hdf5 files
    // // Good condition --> Forced to memory if( workmem == true || ( dsizeA[0]<bigmat && dsizeB[0]<bigmat && dsizeA[1]<bigmat && dsizeB[1]<bigmat))
    // {
      //..// Rcpp::Rcout<<"Working in memory...";
      

      
      /**********************************/
      /**** START IN-MEMORY PROCESSING **/
      /**********************************/
  
      // Read DelayedArray's a and b
      if ( a.isS4() == true)    
      {
        A = read_DelayedArray(a);
      } else {
        try{  
          if ( TYPEOF(a) == INTSXP ) {
            A = Rcpp::as<Eigen::MatrixXi>(a).cast<double>();
          } else{
            A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
          }
        }
        catch(std::exception &ex) { }
      }
      
      if ( w.isS4() == true) {
        B = read_DelayedArray(w);
      }  else {
        
        if ( TYPEOF(w) == INTSXP ) {
          B = Rcpp::as<Eigen::MatrixXi>(w).cast<double>();
        } else{
          B = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(w);
        }
        
      } 
      
      if(bparal == true) {
        C = Bblock_weighted_crossprod_parallel(A, B, iblock_size, threads);
      } else if (bparal == false)  {
        C = Bblock_weighted_crossprod(A, B, iblock_size);
      }
      
      
      return List::create(Named("matrix") = wrap(C));
      
      /********************************/
      /**** END IN-MEMORY PROCESSING **/
      /********************************/
  
    // } 
    // else 
    // {
      
      /********************************/
      /**** START ON-DISK PROCESSING **/
      /********************************/
      
     // TODO 
  
    // }
    
  } catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }
    
  
  // //..// return(C);
  // return List::create(Named("filename") = filename,
  //                     Named("dataset") = strsubgroupOut + "/C");
}



/***R

library(BigDataStatMeth)
library(microbenchmark)

n <- 1024

A <- matrix(runif(n*n), nrow = n, ncol = n)
B <- matrix(runif(n*n), nrow = n, ncol = n)


# D <- A%*%B%*%t(A)
# C <- Crossprod_Weighted(A,B,paral = FALSE)
# all.equal(D, C$matrix)



res <- microbenchmark(R <- A%*%B%*%t(A),
                      Serie<- tCrossprod_Weighted(A,B,paral = FALSE), 
                      Par_2cor<-tCrossprod_Weighted(A,B,paral = TRUE, block_size = 256, threads = 2),
                      Par_3cor<-tCrossprod_Weighted(A,B,paral = TRUE, block_size = 256, threads = 3),
                      Par_4cor<-tCrossprod_Weighted(A,B,paral = TRUE, block_size = 256, threads = 4),
               times = 3 )

res


microbenchmark(A%*%B%*%t(A), 
               Crossprod_Weighted(A,B,paral = FALSE))




C <- Crossprod_Weighted(A,B,paral = FALSE)


*/