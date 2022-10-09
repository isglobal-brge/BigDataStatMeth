#include "include/parallelBlockMult_hdf5.h"

// Documentació : http://www.netlib.org/lapack/lawnspdf/lawn129.pdf


// Working with C-matrix in memory
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
    
    // Open file and get dataset
    H5File* file = new H5File( filename, H5F_ACC_RDWR );
    
    DataSet* datasetA = new DataSet(file->openDataSet(strsubgroup + "/A"));
    DataSet* datasetB = new DataSet(file->openDataSet(strsubgroup + "/B"));
    
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
          Eigen::MatrixXd A = GetCurrentBlock_hdf5( file, datasetA , ii, kk,
                                                    std::min(hdf5_block,isize),std::min(hdf5_block,ksize));
          Eigen::MatrixXd B = GetCurrentBlock_hdf5( file, datasetB, kk, jj,
                                                    std::min(hdf5_block,ksize),std::min(hdf5_block,jsize));

          C.block(ii, jj, std::min(hdf5_block,isize), std::min(hdf5_block,jsize)) = 
              C.block(ii, jj, std::min(hdf5_block,isize), std::min(hdf5_block,jsize)) + A*B;

          if( kk + hdf5_block > K ) ksize = hdf5_block+1;
        }
        
        if( jj + hdf5_block > N ) jsize = hdf5_block+1;
      }
      
      if( ii + hdf5_block > M ) isize = hdf5_block+1;
    }
    
    datasetA->close();
    datasetB->close();
    file->close();
    return(C);
  }else {
      throw std::range_error("non-conformable arguments");
  }
}


// IMPORTANT : R data stored as Col-major in hdf5 file (stored like original data from R)
// Working directly with C-matrix in hdf5 file
// If option paral·lel is enabled, this function loads hdf5 read blocks
// of medium size into memory and calculates the multiplication of blocks
// by applying the parallel algorithm Bblock_matrix_mul_parallel (in-memory process)
int hdf5_block_matrix_mul_hdf5( IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                std::string filename, std::string strsubgroupIN, std::string strsubgroupOUT,
                                int mem_block_size, bool bparal, Rcpp::Nullable<int> threads  = R_NilValue)
{
  int M = sizeA[0]; 
  int K = sizeA[1]; 
  int L = sizeB[0]; 
  int N = sizeB[1]; 
  
  IntegerVector stride = {1,1};
  IntegerVector block = {1,1};
  
  
  if( K == L)
  {
    
    int isize = hdf5_block+1;
    int ksize = hdf5_block+1;
    int jsize = hdf5_block+1;
    
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1); 
    
    // Create an empty dataset for C-matrix into hdf5 file
    //.commented 20201120 - warning check().// int res = create_HDF5_dataset( filename, strsubgroupOUT + "/C", N, M, "real");
    create_HDF5_dataset( filename, strsubgroupOUT + "/C", N, M, "real");
    
    // Open file and get dataset
    H5File* file = new H5File( filename, H5F_ACC_RDWR );
    
    DataSet* datasetA = new DataSet(file->openDataSet(strsubgroupIN + "/A"));
    DataSet* datasetB = new DataSet(file->openDataSet(strsubgroupIN + "/B"));
    DataSet* datasetC = new DataSet(file->openDataSet(strsubgroupOUT + "/C"));
    
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
          Eigen::MatrixXd A = GetCurrentBlock_hdf5( file, datasetA, ii, kk, 
                                                    std::min(hdf5_block,isize),std::min(hdf5_block,ksize));
          Eigen::MatrixXd B = GetCurrentBlock_hdf5( file, datasetB, kk, jj, 
                                                    std::min(hdf5_block,ksize),std::min(hdf5_block,jsize));
          
          Eigen::MatrixXd C = GetCurrentBlock_hdf5( file, datasetC, ii, jj, 
                                                    std::min(hdf5_block,isize),std::min(hdf5_block,jsize));
          
          if( bparal == false)
            C = C + A*B;
          else
            C = C + Bblock_matrix_mul_parallel(A, B, mem_block_size, threads);
          
          
          IntegerVector count = {std::min(hdf5_block,isize), std::min(hdf5_block,jsize)};
          IntegerVector offset = {ii,jj};
          
          write_HDF5_matrix_subset_v2( file, datasetC, offset, count, stride, block, Rcpp::wrap(C));
          
          if( kk + hdf5_block > K ) ksize = hdf5_block+1;
        }
        
        if( jj + hdf5_block > N ) jsize = hdf5_block+1;
      }
      
      if( ii + hdf5_block > M ) isize = hdf5_block+1;
    }
    
    datasetA->close();
    datasetB->close();
    datasetC->close();
    file->close();
    
    return(0);
  }else {
    throw std::range_error("non-conformable arguments");
  }
}







// IMPORTANT : R data stored as Row-major in hdf5 file (stored transposed data from R)
// Working directly with C-matrix in hdf5 file
// If option paral·lel is enabled, this function loads hdf5 read blocks
// of medium size into memory and calculates the multiplication of blocks
// by applying the parallel algorithm Bblock_matrix_mul_parallel (in-memory process)
// browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
//             else, indicates that R data is stored in hdf5 as column major
int hdf5_block_matrix_mul_hdf5_transposed( IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                    std::string filename, std::string strsubgroupIN, std::string strsubgroupOUT, 
                                    int mem_block_size, bool bparal, bool browmajor, 
                                    Rcpp::Nullable<int> threads  = R_NilValue)
{

  // Original data size without transpose 
  //  sizeA[0] :rows(A), sizeA[1] :cols(A), 
  //  sizeB[0] :rows(B), sizeB[1] :cols(B)
  
  int N = sizeA[0];
  int K = sizeA[1];
  
  int L = sizeB[0];
  int M = sizeB[1];
  
  IntegerVector stride = {1,1};
  IntegerVector block = {1,1};

  if( K == L)
  {
    
    int isize = hdf5_block+1;
    int ksize = hdf5_block+1;
    int jsize = hdf5_block+1;
    
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1); 
    
    // Create an empty dataset for C-matrix into hdf5 file
    //.commented 20201120 - warning check().// int res = create_HDF5_dataset( filename, strsubgroupOUT + "/C", M, N, "real");
    create_HDF5_dataset( filename, strsubgroupOUT + "/C", M, N, "real");
    
    // Open file and get dataset
    H5File* file = new H5File( filename, H5F_ACC_RDWR );
    
    DataSet* datasetA = new DataSet(file->openDataSet(strsubgroupIN + "/A"));
    DataSet* datasetB = new DataSet(file->openDataSet(strsubgroupIN + "/B"));
    DataSet* datasetC = new DataSet(file->openDataSet(strsubgroupOUT + "/C"));
    
    // Això haurien de ser les columnes de la matriu i no les files com ara....
    for (int ii = 0; ii < N; ii += hdf5_block)
    {
      
      if( ii + hdf5_block > N ) isize = N - ii;
      // Això haurien de ser files i no per columnes
      for (int jj = 0; jj < M; jj += hdf5_block)
      {
        
        if( jj + hdf5_block > M) jsize = M - jj;
        
        for(int kk = 0; kk < K; kk += hdf5_block)
        {
          if( kk + hdf5_block > K ) ksize = K - kk;
          
          // Get blocks from hdf5 file
          Eigen::MatrixXd A = GetCurrentBlock_hdf5( file, datasetA, kk, ii, 
                                                    std::min(hdf5_block,ksize),std::min(hdf5_block,isize));
          Eigen::MatrixXd B = GetCurrentBlock_hdf5( file, datasetB, jj, kk, 
                                                    std::min(hdf5_block,jsize),std::min(hdf5_block,ksize));
          Eigen::MatrixXd C = GetCurrentBlock_hdf5( file, datasetC, jj, ii, 
                                                    std::min(hdf5_block,jsize),std::min(hdf5_block,isize));
          
          if( bparal == false)
            C = C + B*A;
          else
            C = C + Bblock_matrix_mul_parallel(B, A, mem_block_size, threads);

          IntegerVector count = {std::min(hdf5_block,jsize), std::min(hdf5_block,isize)};
          IntegerVector offset = {jj,ii};
          
          write_HDF5_matrix_subset_v2( file, datasetC, offset, count, stride, block, Rcpp::wrap(C));
          
          if( kk + hdf5_block > K ) ksize = hdf5_block + 1;
        }
        
        if( jj + hdf5_block > M ) jsize = hdf5_block + 1;
      }
      
      if( ii + hdf5_block > N ) isize = hdf5_block + 1;
    }
    
    datasetA->close();
    datasetB->close();
    datasetC->close();
    file->close();
    
    return(0);
  }else {
    throw std::range_error("non-conformable arguments");
  }
}




// IMPORTANT : R data stored as Row-major in hdf5 file (stored transposed data from R)
// 
// This function differs from "hdf5_block_matrix_mul_hdf5_transposed" in parameters, in hdf5_block_matrix_mul_hdf5_transposed we have
// the matrix sizes and here we have the dataset name inside the hdf5 data file
// 
// 
// Working directly with C-matrix in hdf5 file
// If option paral·lel is enabled, this function loads hdf5 read blocks
// of medium size into memory and calculates the multiplication of blocks
// by applying the parallel algorithm Bblock_matrix_mul_parallel (in-memory process)
// browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
//             else, indicates that R data is stored in hdf5 as column major
int hdf5_block_matrix_mul_hdf5_indatasets_transposed( std::string matA, std::string matB, 
                                                      IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                                      std::string filename, std::string strsubgroupIN, 
                                                      std::string strsubgroupINB, std::string strsubgroupOUT, 
                                                      std::string strdatasetOUT, 
                                                      int mem_block_size, bool bparal, bool browmajor, 
                                                      Rcpp::Nullable<int> threads  = R_NilValue)
{
  
  // Original data size without transpose 
  //  sizeA[0] :rows(A), sizeA[1] :cols(A), 
  //  sizeB[0] :rows(B), sizeB[1] :cols(B)
  
  int K = sizeA[0];
  int N = sizeA[1];
  
  int M = sizeB[0];
  int L = sizeB[1];
  
  
  //..// Rcpp::Rcout<<"Estem tractant unes matriu de tamany : \n\tMatriu A : "<<sizeA<<"\n\tMatriu B : "<<sizeB<<"\n\tMatriu C : "<<M<<" "<<N;
  
  IntegerVector stride = {1,1};
  IntegerVector block = {1,1};
  
  if( K == L)
  {
    
    int isize = hdf5_block+1;
    int ksize = hdf5_block+1;
    int jsize = hdf5_block+1;
    
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1); 
    
    // Create an empty dataset for C-matrix into hdf5 file
    //.commented 20201120 - warning check().// int res = create_HDF5_dataset( filename, strsubgroupOUT + "/C", M, N, "real");
    
    //..// Rcpp::Rcout<<"\n Crearem dataset "<<  strsubgroupOUT + matA + "_x_" + matB<<"\n";
    //.2021/12/13.// create_HDF5_dataset( filename, strsubgroupOUT + matA + "_x_" + matB, M, N, "real");
    create_HDF5_dataset( filename, strsubgroupOUT + strdatasetOUT, M, N, "real");
    //..// Rcpp::Rcout<<"Dataset Creat\n";  
    // Open file and get dataset
    H5File* file = new H5File( filename, H5F_ACC_RDWR );
    
    DataSet* datasetA = new DataSet(file->openDataSet(strsubgroupIN + matA));
    DataSet* datasetB = new DataSet(file->openDataSet(strsubgroupINB + matB));
    //.2021/12/13.// DataSet* datasetC = new DataSet(file->openDataSet(strsubgroupOUT + matA + "_x_" + matB));
    DataSet* datasetC = new DataSet(file->openDataSet(strsubgroupOUT + strdatasetOUT));
    

    // Això haurien de ser les columnes de la matriu i no les files com ara....
    for (int ii = 0; ii < N; ii += hdf5_block)
    {
      
      if( ii + hdf5_block > N ) isize = N - ii;
      // Això haurien de ser files i no per columnes
      for (int jj = 0; jj < M; jj += hdf5_block)
      {
        
        if( jj + hdf5_block > M) jsize = M - jj;
        
        for(int kk = 0; kk < K; kk += hdf5_block)
        {
          if( kk + hdf5_block > K ) ksize = K - kk;
          
          // Get blocks from hdf5 file
          Eigen::MatrixXd A = GetCurrentBlock_hdf5( file, datasetA, kk, ii, 
                                                    std::min(hdf5_block,ksize),std::min(hdf5_block,isize));
          Eigen::MatrixXd B = GetCurrentBlock_hdf5( file, datasetB, jj, kk, 
                                                    std::min(hdf5_block,jsize),std::min(hdf5_block,ksize));
          Eigen::MatrixXd C = GetCurrentBlock_hdf5( file, datasetC, jj, ii, 
                                                    std::min(hdf5_block,jsize),std::min(hdf5_block,isize));
          
          if( bparal == false)
            C = C + B*A;
          else
            C = C + Bblock_matrix_mul_parallel(B, A, mem_block_size, threads);
          
          IntegerVector count = {std::min(hdf5_block,jsize), std::min(hdf5_block,isize)};
          IntegerVector offset = {jj,ii};
          
          write_HDF5_matrix_subset_v2( file, datasetC, offset, count, stride, block, Rcpp::wrap(C));
          
          if( kk + hdf5_block > K ) ksize = hdf5_block + 1;
        }
        
        if( jj + hdf5_block > M ) jsize = hdf5_block + 1;
      }
      
      if( ii + hdf5_block > N ) isize = hdf5_block + 1;
    }
    
    datasetA->close();
    datasetB->close();
    datasetC->close();
    file->close();
    
    return(0);
  }else {
    throw std::range_error("non-conformable arguments");
  }
}





int hdf5_block_matrix_mul_parallel( IntegerVector sizeA, IntegerVector sizeB, int block_size, 
                                               std::string filename, std::string strsubgroup, 
                                               Rcpp::Nullable<int> threads  = R_NilValue)
{
  //..// int ii=0, jj=0, kk=0;
  int chunk = 1, tid;
  unsigned int ithreads;
  
  int M = sizeA[0]; 
  int K = sizeA[1]; 
  int L = sizeB[0]; 
  int N = sizeB[1]; 
  
  
  if( K == L)
  {
    
    IntegerVector stride = {1,1};
    IntegerVector block = {1,1};
    
    if(block_size > std::min( N, std::min(M,K)) )
      block_size = std::min( N, std::min(M,K)); 
    
    // Open file for read
    H5File* file = new H5File( filename, H5F_ACC_RDWR );
    
    DataSet* datasetA = new DataSet(file->openDataSet(strsubgroup + "/A"));
    DataSet* datasetB = new DataSet(file->openDataSet(strsubgroup + "/B"));
    DataSet* datasetC = new DataSet(file->openDataSet(strsubgroup + "/C"));
    
    if(threads.isNotNull()) 
    {
      if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()) {
        ithreads = Rcpp::as<int> (threads);
      } else {
        ithreads = getDTthreads(0, true);
        //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;}
      }
    } else {
        ithreads = getDTthreads(0, true);
        //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;
    }
    
    //.OpenMP.// omp_set_dynamic(1);   // omp_set_dynamic(0); omp_set_num_threads(4);
    //.OpenMP.// omp_set_num_threads(ithreads);
    
    // H5File sharedfile = Open_hdf5_file(filename);
    
    //.OpenMP.// #pragma omp parallel shared( file, datasetA, datasetB, datasetC, chunk) private(ii, jj, kk, tid ) 
#pragma omp parallel num_threads(getDTthreads(ithreads, false)) shared( file, datasetA, datasetB, datasetC, chunk) private( tid ) 
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
          Eigen::MatrixXd A = GetCurrentBlock_hdf5( file, datasetA , ii, kk,
                                                    std::min(block_size,M - ii), std::min(block_size,K - kk));
          Eigen::MatrixXd B = GetCurrentBlock_hdf5( file, datasetB, kk, jj,
                                                    std::min(block_size,K - kk), std::min(block_size,N - jj));
          
          Eigen::MatrixXd C = GetCurrentBlock_hdf5( file, datasetC, ii, jj, 
                                                    std::min(block_size,M - ii),std::min(block_size,N - jj));
          C = C + A*B;
          
          IntegerVector count = {std::min(block_size,M - ii), std::min(block_size,N - jj)};
          IntegerVector offset = {ii,jj};
          
          write_HDF5_matrix_subset_v2( file, datasetC, offset, count, stride, block, Rcpp::wrap(C));

        }
      }
    }
  }
  
  return(0);
  
  }else {
    throw std::range_error("non-conformable arguments");
    return(-1);
  }
}




// In-memory execution - Serial version
Eigen::MatrixXd Bblock_matrix_mul(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, int block_size)
{
  
  
  // int M = A.rows();
  // int K = A.cols();
  // int N = B.cols();
  Eigen::MatrixXd C;
    
  try{
    
    int M = A.rows();
    int K = A.cols();
    int N = B.cols();
    
    if( A.cols()==B.rows())
    {
      C = Eigen::MatrixXd::Zero(M,N) ; 
      
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
      
      
    } else {
      throw std::range_error("non-conformable arguments");
    }    
    
  } catch(std::exception &ex) {
    Rcpp::Rcout<<"c++ error : Bblock_matrix_mul : " <<ex.what();
    return(C);
    
  } catch(...) { 
    ::Rf_error("c++ exception in Bblock_matrix_mul (unknown reason)"); 
  }
  

  return(C);
  
}




// In-memory execution - Parallel version
Eigen::MatrixXd Bblock_matrix_mul_parallel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, 
                                          int block_size, Rcpp::Nullable<int> threads  = R_NilValue)
{
  //..// int ii=0, jj=0, kk=0;
  int chunk = 1, tid;
  unsigned int ithreads;
  int M = A.rows();
  int K = A.cols();
  int N = B.cols();
  
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ;
  if(block_size > std::min( N, std::min(M,K)) )
    block_size = std::min( N, std::min(M,K)); 
  
  // if(threads.isNotNull()) 
  // {
  //   if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
  //     ithreads = (Rcpp::as<int> (threads)) * 1025;
  //   else 
  //     ithreads = (std::thread::hardware_concurrency()/2)*1025;
  // }
  // else    ithreads = (std::thread::hardware_concurrency()/2)*1025; //omp_get_max_threads();
  
  
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
  
  
  //.OpenMP.//omp_set_dynamic(1);   // omp_set_dynamic(0); omp_set_num_threads(4);
  //.OpenMP.//omp_set_num_threads(ithreads);
  
  //.OpenMP.//#pragma omp parallel shared(A, B, C, chunk) private(ii, jj, kk, tid ) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C, chunk) private(tid ) 
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
return(C);
}



//' Block matrix multiplication 
//' 
//' This function performs a block matrix-matrix multiplication with numeric matrix
//' 
//' @param a a double matrix.
//' @param b a double matrix.
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @param bigmatrix (optiona, default = 5000) maximum number of rows or columns to consider as big matrix and work with
//' hdf5 files, by default a matrix with more than 5000 rows or files is considered big matrix and computation is made in disk 
//' @param mixblock_size (optiona, default = 128), only if we are working with big matrix and parallel computation = true. 
//' Block size for mixed computation in big matrix parallel. Size of the block to be used to perform parallelized memory 
//' memory of the block read from the disk being processed.
//' @param outfile (optional) file name to work with hdf5 if we are working with big matrix in disk.
//' @param onmemory (optional) if onmemory = TRUE the multiplication is forced to execute in memory
//' @return A List with : 
//' \itemize{
//'   \item{"matrix"}{ Result matrix if execution has been performed in memory}
//'   \item{"filename"}{ HDF5 filename if execution has been performed in disk, HDF5 file contains : 
//'     \itemize{
//'       \item{"INPUT"}{hdf5 group with input matrix A and B}
//'       \item{"OUTPUT"}{hdf5 group with output matrix C}
//'     }with input and output matrix.
//'   }
//' }
//' @examples
//' 
//' # with numeric matrix
//' m <- 500
//' k <- 1500
//' n <- 400
//' A <- matrix(rnorm(n*k), nrow=n, ncol=k)
//' B <- matrix(rnorm(n*k), nrow=k, ncol=n)
//' 
//' C <- blockmult(A,B,128, TRUE)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List blockmult(Rcpp::RObject a, Rcpp::RObject b, 
                              Rcpp::Nullable<int> block_size = R_NilValue, 
                              Rcpp::Nullable<bool> paral = R_NilValue,
                              Rcpp::Nullable<int> threads = R_NilValue,
                              Rcpp::Nullable<double> bigmatrix = R_NilValue,
                              Rcpp::Nullable<double> mixblock_size = R_NilValue,
                              Rcpp::Nullable<std::string> outfile = R_NilValue,
                              Rcpp::Nullable<bool> onmemory = R_NilValue )
{
  
  int iblock_size, res, bigmat;
  bool bparal, workmem, transposed = false; 

  //..// std::string strsubgroup = "Base.matrices/";
  std::string strsubgroupIn = "INPUT/";
  std::string strsubgroupOut = "OUTPUT/";
  std::string filename;
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  Eigen::MatrixXd C;
  
  H5File* file = nullptr;
  
  IntegerVector dsizeA = {0, 0}, 
                dsizeB = {0, 0};
  
  
  // Rcpp::Rcout<<"\n Tipus de dades :  "<<TYPEOF(a)<<"\n";
  // Rcpp::Rcout<<"\n Clase objecte :  "<<  as<std::string>(a.slot("class"))  <<"\n";
  
  
  try{
    
    H5::Exception::dontPrint();
      
    // hdf5 parameters
    
    if( outfile.isNull()) {
      filename = "tmp_blockmult.hdf5";
    } else {
      filename = Rcpp::as<std::string> (outfile);
    }
    
    // Remove old file (if exists)
    RemoveFile(filename);
    
    if( bigmatrix.isNull()) {
      bigmat = 10000;
    } else {
      bigmat = Rcpp::as<double> (bigmatrix);
    }
    
    if( onmemory.isNull()) {
      workmem = false;
    } else {
      workmem = Rcpp::as<bool> (onmemory);
    }
    
  
    // Get matrix sizes
    if ( a.isS4() == true )    
    {
      // dsizeA = get_DelayedArray_size(a);
      throw("Only numeric matrix allowd");
    } else { 
      try{  
        if ( TYPEOF(a) == INTSXP ) {
          dsizeA[0] = Rcpp::as<IntegerMatrix>(a).nrow();
          dsizeA[1] = Rcpp::as<IntegerMatrix>(a).ncol();
          
          if(dsizeA[0] == 0 || dsizeA[1] == 0) {
            dsizeA[0] = Rcpp::as<IntegerVector>(a).size();
            dsizeA[1] = 1;
          }
          
        }else{
          dsizeA[0] = Rcpp::as<NumericMatrix>(a).nrow();
          dsizeA[1] = Rcpp::as<NumericMatrix>(a).ncol();
          if(dsizeA[0] == 0 || dsizeA[1] == 0) {
            dsizeA[0] = Rcpp::as<NumericVector>(a).size();
            dsizeA[1] = 1;
          }
        }
      }catch(std::exception &ex) { }
    }
    
    if ( b.isS4() == true)    
    {
      // dsizeB = get_DelayedArray_size(b);
      throw("Only numeric matrix allowd");
    } else { 
      try{  
        if ( TYPEOF(b) == INTSXP ) {
          dsizeB[0] = Rcpp::as<IntegerMatrix>(b).nrow();
          dsizeB[1] = Rcpp::as<IntegerMatrix>(b).ncol();
          if(dsizeB[0] == 0 || dsizeB[1] == 0) {
            dsizeB[0] = Rcpp::as<IntegerVector>(b).size();
            dsizeB[1] = 1;
          }
        }else{
          dsizeB[0] = Rcpp::as<NumericMatrix>(b).nrow();
          dsizeB[1] = Rcpp::as<NumericMatrix>(b).ncol();
          if(dsizeB[0] == 0 || dsizeB[1] == 0) {
            dsizeB[0] = Rcpp::as<NumericVector>(b).size();
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
    if( workmem == true || ( dsizeA[0]<bigmat && dsizeB[0]<bigmat && dsizeA[1]<bigmat && dsizeB[1]<bigmat))
    {
      //..// Rcpp::Rcout<<"Working in memory...";
      
      filename = "";
      
      /**********************************/
      /**** START IN-MEMORY PROCESSING **/
      /**********************************/
  
      // Read a and b
      if ( a.isS4() == true)    
      {
        // A = read_DelayedArray(a);
        throw("Only numeric matrix allowd");
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
      
      if ( b.isS4() == true) {
        // B = read_DelayedArray(b);
        throw("Only numeric matrix allowd");
      }  else {
        
        if ( TYPEOF(b) == INTSXP ) {
          B = Rcpp::as<Eigen::MatrixXi>(b).cast<double>();
        } else{
          B = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(b);
        }
        
      } 
      
      
      if(bparal == true) {
        C = Bblock_matrix_mul_parallel(A, B, iblock_size, threads);
      } else if (bparal == false)  {
        C = Bblock_matrix_mul(A, B, iblock_size);
      }
      
      
      return List::create(Named("matrix") = wrap(C), 
                          Named("filename") = filename);
      
      /********************************/
      /**** END IN-MEMORY PROCESSING **/
      /********************************/
  
    } 
    else 
    {
      
      // /** https://www.bioconductor.org/packages/release/bioc/vignettes/rhdf5/inst/doc/rhdf5.html **/
      
      /********************************/
      /**** START ON-DISK PROCESSING **/
      /********************************/
      
      //..// Rcpp::Rcout<<"Working on disk ...";
      if(!ResFileExist(filename))  {
        res = Create_hdf5_file(filename);
      }
      
      // Open file to work with
      file = new H5File( filename, H5F_ACC_RDWR );
      
      // Read  a and b
      if ( a.isS4() == true)    
      {
        //.Created before.// res = Create_hdf5_file(filename);
        //.Created before.// res = create_HDF5_group(filename, strsubgroupIn );
        
        //.. 20221009 ..// res = create_HDF5_group_ptr(file, strsubgroupIn );
        //.. 20221009 ..// write_DelayedArray_to_hdf5(filename, strsubgroupIn + "/A", a);
        throw("Only numeric matrix allowd");
  
      } else {
        
        try{  
          
          if ( TYPEOF(a) == INTSXP ) {
            //.Created before.// res = Create_hdf5_file(filename);
            //.Created before.// res = create_HDF5_group(filename, strsubgroupIn );
            res = create_HDF5_group_ptr(file, strsubgroupIn );
            
            //.Created before.// write_HDF5_matrix(filename, strsubgroupIn + "/A", Rcpp::as<IntegerMatrix>(a));
            write_HDF5_matrix_from_R_ptr(file, strsubgroupIn + "/A", Rcpp::as<IntegerMatrix>(a), transposed);
          } else{
            //.Created before.// res = Create_hdf5_file(filename);
            //.Created before.// res = create_HDF5_group(filename, strsubgroupIn );
            res = create_HDF5_group_ptr(file, strsubgroupIn );
            
            //.Created before.// write_HDF5_matrix(filename, strsubgroupIn + "/A", Rcpp::as<NumericMatrix>(a));
            write_HDF5_matrix_from_R_ptr(file, strsubgroupIn + "/A", Rcpp::as<NumericMatrix>(a), transposed );
            
          }
        }
        catch(std::exception &ex) { }
      }
      
      
      if ( b.isS4() == true)    
      {
        //.. 20221009 ..// if(!ResFileExist(filename))  {
        //.. 20221009 ..//   //.Created before.// res = Create_hdf5_file(filename);
        //.. 20221009 ..//   //.Created before.// res = create_HDF5_group(filename, strsubgroupIn );
        //.. 20221009 ..//   res = create_HDF5_group_ptr(file, strsubgroupIn );
        //.. 20221009 ..// }
        //.. 20221009 ..// write_DelayedArray_to_hdf5(filename, strsubgroupIn + "/B", b);
        throw("Only numeric matrix allowd");
        
      } else {
        
        try{  
          if ( TYPEOF(b) == INTSXP ) {
            if(!ResFileExist(filename))  {
              //.Created before.// res = Create_hdf5_file(filename);
              //.Created before.// res = create_HDF5_group(filename, strsubgroupIn );
              res = create_HDF5_group_ptr(file, strsubgroupIn );
            }
            
            //.Created before.// write_HDF5_matrix(filename, strsubgroupIn + "/B", Rcpp::as<IntegerMatrix>(b));
            write_HDF5_matrix_from_R_ptr(file, strsubgroupIn + "/B", Rcpp::as<IntegerMatrix>(b), transposed);
            
          } else{
            if(!ResFileExist(filename))  {
              //.Created before.// res = Create_hdf5_file(filename);
              //.Created before.// res = create_HDF5_group(filename, strsubgroupIn );
              res = create_HDF5_group_ptr(file, strsubgroupIn );
            }
            //.Created before.// write_HDF5_matrix(filename, strsubgroupIn + "/B", Rcpp::as<NumericMatrix>(b));
            write_HDF5_matrix_from_R_ptr(file, strsubgroupIn + "/B", Rcpp::as<NumericMatrix>(b), transposed);
          }
        }
        catch(std::exception &ex) { }
      } 
      
      //.Created before.// res = create_HDF5_group(filename, strsubgroupOut );
      res = create_HDF5_group_ptr(file, strsubgroupOut );
  
      if(bparal == true)
      {
        //.. TODO : Work with parallel hdf5 access
        //..// C = hdf5_block_matrix_mul_parallel( dsizeA, dsizeB, iblock_size, filename, strsubgroup, threads );
        
        int memory_block; // Block size to apply to read hdf5 data to paralelize calculus
        
        if(mixblock_size.isNotNull())
          memory_block = Rcpp::as<int> (mixblock_size);
        else 
          memory_block = 128;
        
        // Test mix versión read block from file and calculate multiplication in memory (with paral·lel algorithm)
        //.commented 20201120 - warning check().// int i = hdf5_block_matrix_mul_hdf5_transposed( dsizeA, dsizeB, iblock_size, filename, strsubgroupIn, strsubgroupOut, 
        //.commented 20201120 - warning check().//                                                memory_block, bparal,true, threads);
        hdf5_block_matrix_mul_hdf5_transposed( dsizeA, dsizeB, iblock_size, filename, strsubgroupIn, strsubgroupOut, 
                                               memory_block, bparal,true, threads);
        
        
        C = Eigen::MatrixXd::Zero(2,2);
        
      }else if (bparal == false)
      {
        
        // Not parallel
        //.commented 20201120 - warning check().// int i = hdf5_block_matrix_mul_hdf5_transposed( dsizeA, dsizeB, iblock_size, filename, strsubgroupIn, strsubgroupOut, 
        //.commented 20201120 - warning check().//                                                0, bparal,true, threads);
        hdf5_block_matrix_mul_hdf5_transposed( dsizeA, dsizeB, iblock_size, filename, strsubgroupIn, strsubgroupOut, 
                                               0, bparal,true, threads);
        C = Eigen::MatrixXd::Zero(2,2);
        
      }
      
      file->close();
  
    }
    
  }catch( FileIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception blockmult (File IException)" );
    return wrap(-1);
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    file->close();
    ::Rf_error( "c++ exception blockmult (DataSet IException)" );
    return wrap(-1);
  } catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }
    
  
  //..// return(C);
  return List::create(Named("filename") = filename,
                      Named("dataset") = strsubgroupOut + "/C");
}



