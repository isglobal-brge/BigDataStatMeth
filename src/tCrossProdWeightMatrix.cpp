#include "include/tCrossProdWeightMatrix.h"

using namespace std;
using namespace Rcpp;

// In-memory execution - Serial version
Eigen::MatrixXd Bblock_weighted_tcrossprod(const Eigen::MatrixXd& A, Eigen::MatrixXd& B, int block_size)
{

  int M = A.rows();
  int K = A.cols();
  int N = B.cols();
  if( A.cols() == B.rows() && B.cols() == B.rows())
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
Eigen::MatrixXd Bblock_weighted_tcrossprod_parallel(const Eigen::MatrixXd& A, Eigen::MatrixXd& B, 
                                                   int block_size, Rcpp::Nullable<int> threads  = R_NilValue)
{
  int ii=0, jj=0, kk=0;
  int chunk = 1;
  int tid;
  unsigned int ithreads;
  int M = A.rows();
  int K = A.cols();
  int N = B.cols();
  
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ;
  if(block_size > std::min( N, std::min(M,K)) )
    block_size = std::min( N, std::min(M,K)); 
  
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

  //.OpenMP.//  omp_set_dynamic(0);   // omp_set_dynamic(0); omp_set_num_threads(4);
  //.OpenMP.//  omp_set_num_threads(ithreads);
  
//.OpenMP.// #pragma omp parallel shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
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
  
  //.OpenMP.//  omp_set_dynamic(0);   
  //.OpenMP.//  omp_set_num_threads(ithreads);

//.OpenMP.// #pragma omp parallel shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
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



//' Block matrix multiplication
//' 
//' This function performs a Crossproduct with weigths matrix A%*%W%*%t(A) multiplication with numeric matrix
//' 
//' @param A a double matrix.
//' @param W a Weighted matrix
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @return Matrix with A%*%W%*%t(A) product 
//' @examples
//' 
//' # with numeric matrix
//' m <- 500
//' k <- 1500
//' n <- 400
//' A <- matrix(rnorm(n*k), nrow=n, ncol=k)
//' B <- matrix(rnorm(n*k), nrow=k, ncol=n)
//' 
//' # Serial execution
//' Serie<- bdtCrossprod_Weighted(A, B, paral = FALSE)
//' 
//' # Parallel execution with 2 threads and blocks 256x256
//' Par_2cor <- bdtCrossprod_Weighted(A, B, paral = TRUE, block_size = 256, threads = 2)
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdtCrossprod_Weighted(Rcpp::RObject A, Rcpp::RObject W, 
                                    Rcpp::Nullable<int> block_size = R_NilValue, 
                                    Rcpp::Nullable<bool> paral = R_NilValue,
                                    Rcpp::Nullable<int> threads = R_NilValue )
{
  
  int iblock_size;
  bool bparal; 

  Eigen::MatrixXd mA;
  Eigen::MatrixXd B;
  Eigen::MatrixXd C;
  
  IntegerVector dsizeA = {0, 0}, 
                dsizeB = {0, 0};
  
  // Rcpp::Rcout<<"\n Tipus de dades :  "<<TYPEOF(A)<<"\n";
  // Rcpp::Rcout<<"\n Clase objecte :  "<<  as<std::string>(A.slot("class"))  <<"\n";
  
  
  try{
    
    // Get matrix sizes
    if ( A.isS4() == true )    
    {
      // dsizeA = get_DelayedArray_size(A);
      throw("Only numeric matrix allowd");
    } else { 
      try{  
        if ( TYPEOF(A) == INTSXP ) {
          dsizeA[0] = Rcpp::as<IntegerMatrix>(A).nrow();
          dsizeA[1] = Rcpp::as<IntegerMatrix>(A).ncol();
          
          if(dsizeA[0] == 0 || dsizeA[1] == 0) {
            dsizeA[0] = Rcpp::as<IntegerVector>(A).size();
            dsizeA[1] = 1;
          }
          
        }else{
          dsizeA[0] = Rcpp::as<NumericMatrix>(A).nrow();
          dsizeA[1] = Rcpp::as<NumericMatrix>(A).ncol();
          if(dsizeA[0] == 0 || dsizeA[1] == 0) {
            dsizeA[0] = Rcpp::as<NumericVector>(A).size();
            dsizeA[1] = 1;
          }
        }
      }catch(std::exception &ex) { }
    }
    
    if ( W.isS4() == true)    
    {
      // dsizeB = get_DelayedArray_size(W);
      throw("Only numeric matrix allowd");
    } else { 
      try{  
        if ( TYPEOF(W) == INTSXP ) {
          dsizeB[0] = Rcpp::as<IntegerMatrix>(W).nrow();
          dsizeB[1] = Rcpp::as<IntegerMatrix>(W).ncol();
          if(dsizeB[0] == 0 || dsizeB[1] == 0) {
            dsizeB[0] = Rcpp::as<IntegerVector>(W).size();
            dsizeB[1] = 1;
          }
          
        }else{
          dsizeB[0] = Rcpp::as<NumericMatrix>(W).nrow();
          dsizeB[1] = Rcpp::as<NumericMatrix>(W).ncol();
          if(dsizeB[0] == 0 || dsizeB[1] == 0) {
            dsizeB[0] = Rcpp::as<NumericVector>(W).size();
            dsizeB[1] = 1;
          }
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
      
    // All dimensions are ok before start process 
    if(dsizeB[0]==dsizeB[1] && dsizeA[1]==dsizeB[0])
    {
      
      /**********************************/
      /**** START IN-MEMORY PROCESSING **/
      /**********************************/
      
      // Read A and b
      if ( A.isS4() == true)    
      {
        // mA = read_DelayedArray(A);
        throw("Only numeric matrix allowd");
      } else {
        try{  
          if ( TYPEOF(A) == INTSXP ) {
            mA = Rcpp::as<Eigen::MatrixXi>(A).cast<double>();
          } else{
            mA = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(A);
          }
        }
        catch(std::exception &ex) { }
      }
      
      if ( W.isS4() == true) {
        // B = read_DelayedArray(W);
        throw("Only numeric matrix allowd");
      }  else {
        
        if ( TYPEOF(W) == INTSXP ) {
          B = Rcpp::as<Eigen::MatrixXi>(W).cast<double>();
        } else{
          B = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(W);
        }
        
      } 
      
      if(bparal == true) {
        C = Bblock_weighted_tcrossprod_parallel(mA, B, iblock_size, threads);
      } else if (bparal == false)  {
        C = Bblock_weighted_tcrossprod(mA, B, iblock_size);
      }
      
      
      return (wrap(C));
      
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
      
    } else {
      
      throw std::runtime_error("Dimension error");
      return wrap(-1);
    }
      
 
    
  } catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }
    
  
  // //..// return(C);
  // return List::create(Named("filename") = filename,
  //                     Named("dataset") = strsubgroupOut + "/C");
}



/***R

*/
