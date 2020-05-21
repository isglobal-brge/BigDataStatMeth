#include "include/parallelBlockMult_hdf5.h"

// Documentació : http://www.netlib.org/lapack/lawnspdf/lawn129.pdf


// Convert array with readed data in rowmajor to colmajor matrix : 
//    Eigen and R works with ColMajor and hdf5 in RowMajor
Eigen::MatrixXd RowMajorVector_to_ColMajorMatrix(double* datablock, int countx, int county)
{
  Eigen::MatrixXd mdata(countx, county);

  for (size_t i=0; i<countx; i++)
    for(size_t j=0;j<county;j++)
      mdata(i,j) = datablock[i*county+j];

  return(mdata);
}


// Only one instance reading
Eigen::MatrixXd GetCurrentBlock_hdf5( std::string filename, std::string dataset,
                                      int offsetx,int offsety, 
                                      int countx, int county)
{
  
  IntegerVector offset = IntegerVector::create(offsetx, offsety) ;
  IntegerVector count = IntegerVector::create(countx, county) ;
  IntegerVector stride = IntegerVector::create(1, 1) ;
  IntegerVector block = IntegerVector::create(1, 1) ;
  
  NumericMatrix data(countx, county);
  
  read_HDF5_matrix_subset(filename, dataset, offset, count, stride, block, REAL(data));
  
  Eigen::MatrixXd mat = RowMajorVector_to_ColMajorMatrix(REAL(data), countx, county);

  return( mat );
  
}


/***
// Parallel reading
Eigen::MatrixXd GetCurrentBlock_hdf5_parallel( std::string filename, std::string dataset,
                                      int offsetx,int offsety, 
                                      int countx, int county)
{
  
  IntegerVector offset = IntegerVector::create(offsetx, offsety) ;
  IntegerVector count = IntegerVector::create(countx, county) ;
  IntegerVector stride = IntegerVector::create(1, 1) ;
  IntegerVector block = IntegerVector::create(1, 1) ;
  
  NumericMatrix data(countx, county);
  
  read_HDF5_matrix_subset(filename, dataset, offset, count, stride, block, REAL(data));
  Eigen::MatrixXd mat = as<Eigen::MatrixXd>(data);
  return( mat.transpose() );
  
}
***/



//
// Working with C-matrix in memory
//
Eigen::MatrixXd hdf5_block_matrix_mul( IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                       std::string filename, std::string strsubgroup )
{
  int M = sizeA[0]; 
  int K = sizeA[1]; 
  int L = sizeB[0]; 
  int N = sizeB[1]; 
  

  if( K == L)
  {
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ; 
    
    int isize = hdf5_block+1;
    int ksize = hdf5_block+1;
    int jsize = hdf5_block+1;

    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1); 
    
    for (int ii = 0; ii < M; ii += hdf5_block)
    {

      if( ii + hdf5_block > M ) isize = M - ii;
      for (int jj = 0; jj < N; jj += hdf5_block)
      {
        
        if( jj + hdf5_block > N) jsize = N - jj;
        for(int kk = 0; kk < K; kk += hdf5_block)
        {
          if( kk + hdf5_block > K ) ksize = K - kk;
          
          // Get blocks from hdf5 file
          Eigen::MatrixXd A = GetCurrentBlock_hdf5( filename, strsubgroup + "/A" , ii, kk,
                                                    std::min(hdf5_block,isize),std::min(hdf5_block,ksize));
          Eigen::MatrixXd B = GetCurrentBlock_hdf5( filename, strsubgroup + "/B", kk, jj,
                                                    std::min(hdf5_block,ksize),std::min(hdf5_block,jsize));

          C.block(ii, jj, std::min(hdf5_block,isize), std::min(hdf5_block,jsize)) = 
              C.block(ii, jj, std::min(hdf5_block,isize), std::min(hdf5_block,jsize)) + A*B;

          if( kk + hdf5_block > K ) ksize = hdf5_block+1;
        }
        
        if( jj + hdf5_block > N ) jsize = hdf5_block+1;
      }
      
      if( ii + hdf5_block > M ) isize = hdf5_block+1;
    }
    
    return(C);
  }else {
      throw std::range_error("non-conformable arguments");
  }
}


//
// Working with C-matrix in hdf5 file
// 
int hdf5_block_matrix_mul_v2( IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                       std::string filename, std::string strsubgroup )
{
  int M = sizeA[0]; 
  int K = sizeA[1]; 
  int L = sizeB[0]; 
  int N = sizeB[1]; 
  
  IntegerVector stride = {1,1};
  IntegerVector block = {1,1};
  
  
  if( K == L)
  {
    //Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ; 
    
    {
      NumericMatrix c(M,N);
      write_HDF5_matrix(filename, strsubgroup + "/C", c);  
    }
    
    
    int isize = hdf5_block+1;
    int ksize = hdf5_block+1;
    int jsize = hdf5_block+1;
    
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1); 
    
    for (int ii = 0; ii < M; ii += hdf5_block)
    {
      
      if( ii + hdf5_block > M ) isize = M - ii;
      for (int jj = 0; jj < N; jj += hdf5_block)
      {
        
        if( jj + hdf5_block > N) jsize = N - jj;
        
        for(int kk = 0; kk < K; kk += hdf5_block)
        {
          if( kk + hdf5_block > K ) ksize = K - kk;
          
          // Get blocks from hdf5 file
          Eigen::MatrixXd A = GetCurrentBlock_hdf5( filename, strsubgroup + "/A" , ii, kk,
                                                    std::min(hdf5_block,isize),std::min(hdf5_block,ksize));
          Eigen::MatrixXd B = GetCurrentBlock_hdf5( filename, strsubgroup + "/B", kk, jj,
                                                    std::min(hdf5_block,ksize),std::min(hdf5_block,jsize));
          
        
          Eigen::MatrixXd C = GetCurrentBlock_hdf5( filename, strsubgroup + "/C", ii, jj,
                                                  std::min(hdf5_block,isize),std::min(hdf5_block,jsize));
        
          C = C + A*B;
          
          IntegerVector count = {std::min(hdf5_block,isize), std::min(hdf5_block,jsize)};
          IntegerVector offset = {ii,jj};
          
          write_HDF5_matrix_subset( filename, strsubgroup + "/C", offset, count, stride, block, Rcpp::wrap(C));
          
          if( kk + hdf5_block > K ) ksize = hdf5_block+1;
        }
        
        if( jj + hdf5_block > N ) jsize = hdf5_block+1;
      }
      
      if( ii + hdf5_block > M ) isize = hdf5_block+1;
    }
    
    return(0);
  }else {
    throw std::range_error("non-conformable arguments");
  }
}




Eigen::MatrixXd hdf5_block_matrix_mul_parallel( IntegerVector sizeA, IntegerVector sizeB, int block_size, 
                                               std::string filename, std::string strsubgroup, 
                                               Rcpp::Nullable<int> threads  = R_NilValue)
{
  int ii=0, jj=0, kk=0;
  int chunk = 1, tid;
  unsigned int ithreads;
  
  int M = sizeA[0]; 
  int K = sizeA[1]; 
  int L = sizeB[0]; 
  int N = sizeB[1]; 
  
  
  if( K == L)
  {
    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
  
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ;
    if(block_size > std::min( N, std::min(M,K)) )
      block_size = std::min( N, std::min(M,K)); 
    
    if(threads.isNotNull()) 
    {
      if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
        ithreads = Rcpp::as<int> (threads);
      else 
        ithreads = std::thread::hardware_concurrency();
    }
    else    ithreads = std::thread::hardware_concurrency() - 1; //omp_get_max_threads();
    
    omp_set_dynamic(1);   // omp_set_dynamic(0); omp_set_num_threads(4);
    omp_set_num_threads(ithreads);
    
    // H5File sharedfile = Open_hdf5_file(filename);
    
  #pragma omp parallel shared( A, B, C, chunk) private(ii, jj, kk, tid ) 
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
          // Get blocks from file
          A = GetCurrentBlock_hdf5( filename, strsubgroup + "/A" , ii, kk,
                                                    std::min(block_size,M - ii), std::min(block_size,K - kk));
          B = GetCurrentBlock_hdf5( filename, strsubgroup + "/B", kk, jj,
                                                    std::min(block_size,K - kk), std::min(block_size,N - jj));
          
          C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) = 
            C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) + A*B;
          /***
          C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) = 
            C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) + 
            (A.block(ii, kk, std::min(block_size,M - ii), std::min(block_size,K - kk)) * 
            B.block(kk, jj, std::min(block_size,K - kk), std::min(block_size,N - jj)));
          ***/
        }
      }
    }
  }
  
  return(C);
  
  }else {
    throw std::range_error("non-conformable arguments");
  }
}





Eigen::MatrixXd Bblock_matrix_mul(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size)
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


Eigen::MatrixXd Bblock_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, 
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
      ithreads = std::thread::hardware_concurrency();
  }
  else    ithreads = std::thread::hardware_concurrency() - 1; //omp_get_max_threads();

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
Eigen::MatrixXd Bblockmult_v2(Rcpp::RObject a, Rcpp::RObject b, 
                              Rcpp::Nullable<int> block_size = R_NilValue, 
                              Rcpp::Nullable<bool> paral = R_NilValue,
                              Rcpp::Nullable<int> threads = R_NilValue,
                              Rcpp::Nullable<double> bigmatrix = R_NilValue )
{
  int iblock_size, res, bigmat;
  bool bparal;// = Rcpp::as<double>;
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  Eigen::MatrixXd C;
  
  IntegerVector dsizeA, dsizeB;
  
  // hdf5 parameters
  std::string filename = "tmp_blockmult.hdf5";
  std::string strsubgroup = "Base.matrices/";
  
  // Remove old file
  RemoveFile(filename);
  
  if( bigmatrix.isNull()) {
    bigmat = 5000;
  } else {
    bigmat = Rcpp::as<double> (bigmatrix);
  }
  

  // Get matrix sizes
  if ( a.isS4() == true)    
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
  
  if ( b.isS4() == true)    
  {
    dsizeB = get_DelayedArray_size(b);
  } else { 
    try{  
      if ( TYPEOF(b) == INTSXP ) {
        dsizeB[0] = Rcpp::as<IntegerMatrix>(b).nrow();
        dsizeB[1] = Rcpp::as<IntegerMatrix>(b).ncol();
      }else{
        dsizeB[0] = Rcpp::as<NumericMatrix>(b).nrow();
        dsizeB[1] = Rcpp::as<NumericMatrix>(b).ncol();
      }
    }catch(std::exception &ex) { }
  }
  
  if(block_size.isNotNull())
  {
    iblock_size = Rcpp::as<int> (block_size);
    /*
    if( iblock_size > std::min(dsizeA[0],dsizeA[1]) || iblock_size > std::min(dsizeB[0],dsizeB[1]) )
    {
      iblock_size = std::min(  std::min(dsizeA[0],dsizeA[1]), std::min(dsizeB[0],dsizeB[1]));
    }
     */
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
  if( dsizeA[0]<bigmat && dsizeB[0]<bigmat && dsizeA[1]<bigmat && dsizeB[1]<bigmat)
  {
    Rcpp::Rcout<<"Working in memory...";
    
    /*********************************/
    /**** INICI FUNCIÓ EN MEMÒRIA ****/
    /*********************************/

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
    
    if(bparal == true)
    {
      C = Bblock_matrix_mul_parallel(A, B, iblock_size, threads);
    }else if (bparal == false)
    {
      C = Bblock_matrix_mul(A, B, iblock_size);
    }
    //..// return(C);
    
    
    /******************************/
    /**** FI FUNCIÓ EN MEMÒRIA ****/
    /******************************/

  } 
  else 
  {
    
    Rcpp::Rcout<<"\nWorking with hdf5...\n";
    
    // Read DelayedArray's a and b
    if ( a.isS4() == true)    
    {
      res = Create_hdf5_file(filename);
      res = create_HDF5_group(filename, strsubgroup );
      Eigen::MatrixXd tmp = read_DelayedArray(a);
      write_HDF5_matrix(filename, strsubgroup + "/A", Rcpp::wrap(tmp));
    } else {
      
      try{  
        
        if ( TYPEOF(a) == INTSXP ) {
          res = Create_hdf5_file(filename);
          res = create_HDF5_group(filename, strsubgroup );
          write_HDF5_matrix(filename, strsubgroup + "/A", Rcpp::as<IntegerMatrix>(a));
        } else{
          res = Create_hdf5_file(filename);
          res = create_HDF5_group(filename, strsubgroup );
          write_HDF5_matrix(filename, strsubgroup + "/A", Rcpp::as<NumericMatrix>(a));
        }
      }
      catch(std::exception &ex) { }
    }
    
    
    if ( b.isS4() == true)    
    {
      if(!ResFileExist(filename))  {
        res = Create_hdf5_file(filename);
        res = create_HDF5_group(filename, strsubgroup );
      }
      Eigen::MatrixXd tmp = read_DelayedArray(b);
      write_HDF5_matrix(filename, strsubgroup + "/B", Rcpp::wrap(tmp));
      
    } else {
      
      try{  
        if ( TYPEOF(b) == INTSXP ) {
          if(!ResFileExist(filename))  {
            res = Create_hdf5_file(filename);
            res = create_HDF5_group(filename, strsubgroup );
          }
          
          write_HDF5_matrix(filename, strsubgroup + "/B", Rcpp::as<IntegerMatrix>(b));
        } else{
          if(!ResFileExist(filename))  {
            res = Create_hdf5_file(filename);
            res = create_HDF5_group(filename, strsubgroup );
          }
          write_HDF5_matrix(filename, strsubgroup + "/B", Rcpp::as<NumericMatrix>(b));
        }
      }
      catch(std::exception &ex) { }
    } 
    

    if(bparal == true)
    {
      //.. hdf5 parallel not working only serial access ..// 
      //.. TODO : Work with parallel hdf5 access
      //..// C = hdf5_block_matrix_mul_parallel( dsizeA, dsizeB, iblock_size, filename, strsubgroup, threads );
      
      // Not parallel
      C = hdf5_block_matrix_mul( dsizeA, dsizeB, iblock_size, filename, strsubgroup );
      
    }else if (bparal == false)
    {
      // Not parallel
      int i = hdf5_block_matrix_mul_v2( dsizeA, dsizeB, iblock_size, filename, strsubgroup );
      C = Eigen::MatrixXd::Zero(2,2);
      
    }

  }
  
  return(C);
}





/*** R

library(DelayedArray)
library(BigDataStatMeth)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/BigDataStatMeth/tmp")

N <- 4
M <- 4
set.seed(123)
A <- matrix(rnorm(N*M,mean=0,sd=1), N, M)
Ad <- DelayedArray(A)
B <- matrix(sample.int(15, 20*100, TRUE), 20, 100)

A[1:2,1:2]*A[1:2,1:2]


dim(A)
dim(B)

A[1:5,1:5]
B[1:5,1:5]

C <- Bblockmult_v2(A,B, block_size = 2, bigmatrix = 100)


Bblockmult(Ad,B)

all.equal(C,A%*%B)

Bblockmult(B,B)





# library(microbenchmark)
# library(DelayedArray)
# library(BigDataStatMeth)
# # A <- matrix(sample(1:10,100, replace = TRUE), ncol = 10);A
# 
# n <- 1000
# p <- 400
# A <- matrix(rnorm(n*p), nrow=n, ncol=p)
#   
#   
# n <- 400
# p <- 1000
# B <-  matrix(sample(1:10, 40, replace = TRUE), nrow = n, ncol = p)
# 
# CPP2 <- blockmult(A,B,128, TRUE)
# 
# AD <- DelayedArray(A)
# BD <- DelayedArray(B)
# 
# results <- microbenchmark( CP2 <- blockmult(A,B,128, FALSE),
#                            CPP2 <- blockmult(A,B,128, TRUE),
#                            CP2D <- blockmult(AD,BD,128, FALSE),
#                            CPP2D <- blockmult(AD,BD,128, TRUE),
#                            dgemm <- matv_gemm(A,B),
#                            r <- A%*%B,
#                            times = 2L)
#   
# print(summary(results)[, c(1:7)],digits=3)
# 
# stopifnot(all.equal(CP2,CP2D),
#           all.equal(CPP2,CPP2D),
#           all.equal(CP2,CPP2),
#           all.equal(CP2,dgemm),
#           all.equal(CP2,r))
# 
# 
#   D <- block_matrix_mul_parallel(A,B,128)
#   CP2 <- blockmult(A,B,128)
#   
#   results <- microbenchmark(CPP1 <- blockmult(A,B,64, TRUE),
#                             CPP2 <- blockmult(A,B,128, TRUE),
#                             CPP3 <- blockmult(A,B,256, TRUE),
#                             CPP4 <- blockmult(A,B,512, TRUE),
#                             CPP5 <- blockmult(A,B,1024, TRUE),
#                             CP1 <- blockmult(A,B,64, FALSE),
#                             CP2 <- blockmult(A,B,128, FALSE),
#                             CP3 <- blockmult(A,B,256, FALSE),
#                             CP4 <- blockmult(A,B,512, FALSE),
#                             CP5 <- blockmult(A,B,1024, FALSE),
#                             #C <- A%*%B,
#                                                                   times = 2L)  # Proves multiplicacions x blocs
#   
#   print(summary(results)[, c(1:7)],digits=3)
#   
#   
#   
#   stopifnot(all.equal(block_matrix_mul_parallel(A,B,64),A%*%B), 
#             all.equal(blockmult(A,B,128), A%*%B))
#   


*/
