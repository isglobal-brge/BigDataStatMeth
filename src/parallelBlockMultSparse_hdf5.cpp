#include "include/parallelBlockMultSparse_hdf5.h"

// BASED ON : Gustavson’s algorithm (Colwise version)
//    References : 
//    * F. G. Gustavson, “Two fast algorithms for sparse matrices: Multiplication and permuted transposition” Trans. on Mathematical Software (TOMS),1978.
//    * MatRaptor: A Sparse-Sparse Matrix Multiplication Accelerator Based on Row-Wise Product
//    
//    Pros : Multiple rows from the output matrix can be computed in parallel.


// Parallel matrix multiplication with sparse matrices.¡
Eigen::SparseMatrix<double> matmult_sparse_parallel ( Eigen::Map<Eigen::SparseMatrix<double> > A, 
                                                      Eigen::Map<Eigen::SparseMatrix<double> > B,
                                                      Rcpp::Nullable<int> threads  = R_NilValue )
{
  
  int ii;
  int tid, iblocks; // chunk = 1,
  int ithreads;
  int M = A.rows(), N = B.cols();
  
  Eigen::MatrixXd C(M, N);
  
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
  
  // Columns x thread
  iblocks = N / ithreads;
  
  //.OpenMP.// omp_set_dynamic(0);   // omp_set_dynamic(0); omp_set_num_threads(4);
  //.OpenMP.// omp_set_num_threads(ithreads);
 
 
 //.OpenMP.//#pragma omp parallel shared( A, B, C) private(ii, tid) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared( A, B, C) private(ii, tid) 
{
  
  // tid = omp_get_thread_num();
  // //només per fer proves dels threads i saber que està paralelitzant, sinó no cal tenir-ho descomentat
  // if (tid == 0)   {
  //   Rcpp::Rcout << "Number of threads: " << omp_get_num_threads() << "\n";
  // }
  
#pragma omp for schedule(auto)
  
  for (int ii = 0; ii < ithreads; ii = ii+1)
  {
    // Rcpp::Rcout << "Number of threads: " << omp_get_num_threads() << "\n";
      //..// C_sp.col(ii)  = C_sp.col(ii) + A*B.col(ii);
      if( ii < ithreads-1){
        Eigen::SparseMatrix<double> C_sp = A * B.block( 0, ii*iblocks, M, iblocks );
        C.block( 0, ii*iblocks , M, iblocks ) =  Eigen::MatrixXd (C_sp);;
        
      } else {
        Eigen::SparseMatrix<double> C_sp = A * B.block(0, ii*iblocks, M, N-(ii*iblocks) );
        C.block( 0, ii*iblocks, M, N-(ii*iblocks) ) = Eigen::MatrixXd (C_sp);
      }
  }
  
}
return(C.sparseView());
  
  
  
}









// Simple matrix multiplication with sparse matrices.¡
Eigen::SparseMatrix<double> matmult_sparse( Eigen::Map<Eigen::SparseMatrix<double> > A, 
                                           Eigen::Map<Eigen::SparseMatrix<double> > B)
{
  
  Eigen::SparseMatrix<double> C_sp(A.rows(), B.cols());
  
  // C_sp = Bblock_matrix_mul_parallel(A, B, iblock_size, threads);
  C_sp = A*B; 
  
  return (C_sp);
  
}




//' Block matrix multiplication 
//' 
//' This function performs a block matrix-matrix multiplication with numeric matrix
//' 
//' @param A a sparse double matrix.
//' @param B a sparse double matrix.
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if paral = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @return Sparse matrix with results
//' @examples
//' 
//' library(Matrix)
//' library(BigDataStatMeth)
//' 
//' k <- 1e3
//' set.seed(1)
//' x_sparse <- sparseMatrix(
//'     i = sample(x = k, size = k),
//'     j = sample(x = k, size = k),
//'     x = rnorm(n = k)
//' )
//' set.seed(2)
//' y_sparse <- sparseMatrix(
//'     i = sample(x = k, size = k),
//'     j = sample(x = k, size = k),
//'     x = rnorm(n = k)
//' )
//' 
//' d <- bdblockmult_sparse(x_sparse, y_sparse)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdblockmult_sparse(Rcpp::RObject A, Rcpp::RObject B, Rcpp::Nullable<bool> paral = R_NilValue, Rcpp::Nullable<int> threads  = R_NilValue )
{
  
  std::string spMatrix("dgCMatrix");
  bool bSparseA = true, bSparseB = true;
  bool bparal = false;

  try{
    
    if( spMatrix.compare( as<std::string>(A.slot("class")))!=0 ){
      bSparseA = false;
      Rf_error("Matrix A isn't Sparsed");
    }
    
    if( spMatrix.compare( as<std::string>(B.slot("class")))!=0 ){
      bSparseB = false;
      Rf_error("Matrix B isn't Sparsed");
    }
    
    if( paral.isNull()) {
      bparal = false;
    } else {
      bparal = Rcpp::as<bool> (paral);
    }
    
    
    if( bSparseA == true || bSparseB == true) 
    {
      Eigen::Map<Eigen::SparseMatrix<double> > A_sp = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double> > >(A);
      Eigen::Map<Eigen::SparseMatrix<double> > B_sp = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double> > >(B);
      Eigen::SparseMatrix<double> C_sp(A_sp.rows(), B_sp.cols());
      
      if(bparal == false){
        // C_sp = Bblock_matrix_mul_parallel(A, B, iblock_size, threads);
        C_sp = A_sp*B_sp; 
      
      } else {
        
        C_sp = matmult_sparse_parallel(A_sp, B_sp, threads);
        
      }
      
      return (wrap(C_sp));
      
    } else {
      
      Rf_error("No sparse matrix found");
      return(wrap(-1));
      
    }

  } catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }

   
}



