#include "include/mlr_mr.h"

using namespace RcppEigen;
using namespace Rcpp;


// Apply mlr_mr algorith to perform lm regression with big data
// Algorithm from : 
//    https://isglobal-brge.github.io/Master_Modelling/dealing-with-big-data-in-r.html#linear-regression-for-big-data
Eigen::MatrixXd Rcpp_mlr_mr(Eigen::MatrixXd x, Eigen::MatrixXd y, int iblocks, Rcpp::Nullable<int> threads  = R_NilValue )
{
  
  int chunk = 1; // , tid;
  unsigned int ithreads;
  int irows =  x.rows(), 
      icols = x.cols();
  
  // Get number of threads
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
  

  //.OpenMP.// omp_set_num_threads(ithreads); 
  
  // Define variables
  Eigen::MatrixXd R1(iblocks*x.cols(), x.cols());
  Eigen::MatrixXd Q1( irows, icols );
  Eigen::VectorXd indexQ1( iblocks+1 );
  
//.OpenMP.// #pragma omp parallel shared( R1, Q1, indexQ1, chunk) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared( R1, Q1, indexQ1, chunk) 
  {
    // First steps --> (1) Read block
    //                 (2) Make decomposition by block and reduce result matrix
    
    double block_size = std::floor((double)irows/(double)iblocks);
    int maxsizetoread = block_size;
    int iSizetoRead; 

    
#pragma omp parallel for schedule(static)
    for( int i = 0; i<iblocks ; i++) 
    {
      strQR decQR;
      
      if( i+1 == iblocks && irows - maxsizetoread!=0)
        iSizetoRead = irows - (i*block_size);
      else
        iSizetoRead = maxsizetoread;
      
      // Read current block
      Eigen::MatrixXd obsBlock = GetCurrentBlock( x, i*block_size, 0, iSizetoRead, icols);
  
      // QR Decomposition
      decQR = rcpp_bdQR(obsBlock, true);
      
      // Complete R1
      R1.block( i*x.cols(), 0, decQR.R.rows(), decQR.R.cols() ) = decQR.R;
  
      if (i == 0) { indexQ1(i) = 0; } 
      indexQ1(i+1) = indexQ1(i) + decQR.Q.rows();
  
      // Complete Q1
      Q1.block(indexQ1(i), 0, decQR.Q.rows(), decQR.Q.cols() ) = decQR.Q;
  
    }
  } 

   
  // QR decomposition from R1
  strQR decR1;
  decR1 = rcpp_bdQR(R1, true);
  
  
  Eigen::MatrixXd Q3;
  Eigen::MatrixXd V(icols, iblocks);
  
  
  double block_size = std::floor((double)decR1.Q.rows()/(double)iblocks);
  int maxsizetoread = block_size;
  int iSizetoRead; 
  int startnext = 0;

  
#pragma omp parallel for ordered schedule(dynamic) 
  for( int i = 0; i< iblocks ; i++) 
  {
    
    Eigen::MatrixXd Q1BlockDiv;
    Eigen::MatrixXd Q2BlockDiv;
    
    if( i+1 == iblocks && decR1.Q.rows() - maxsizetoread!=0){
      iSizetoRead = decR1.Q.rows() - (i*block_size);
    }else{
      iSizetoRead = maxsizetoread;
    }

    // Get block from Q2    
    Q2BlockDiv = GetCurrentBlock( decR1.Q, i*block_size, 0, iSizetoRead, decR1.Q.cols());

    // Get block from Q1
    Q1BlockDiv = GetCurrentBlock( Q1, indexQ1(i), 0, ((indexQ1(i+1)) - indexQ1(i)), Q1.cols());

    // Get Q3 and V
    Q3 = Bblock_matrix_mul_parallel(Q1BlockDiv, Q2BlockDiv, 256, threads);
    
    Eigen::MatrixXd YBlock  = GetCurrentBlock( y, indexQ1(i), 0, Q3.rows(), 1);
    
    V.block(0, i, V.rows(), 1) = Bblock_matrix_mul_parallel( Q3.adjoint(), YBlock , 256, threads );
    
    startnext = startnext + Q3.adjoint().cols();

  }


  // Get Betas
  Eigen::MatrixXd beta = Bblock_matrix_mul_parallel( decR1.R.inverse(), V.rowwise().sum(), 256, threads);

  return(beta);
}



// // ' Linear regression using MLR-MR algorithm
// // '
// // ' Linear regression for Big Data using MLR-MR algorithm
// // '
// // ' @param X, numerical matrix with paired observations of the predictor variable X
// // ' @param Y, numerical matrix column with response variable
// // ' @param blocks, integer with number of blocks we want to split matrix if null matrix is splited in blocks as maximum of 1000 variables per block
// // ' @param threads, threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
// // ' @return Lineal regression coefficients
// //' @export
// [[Rcpp::export(.bdMLR_MR)]]
Rcpp::RObject bdMLR_MR(Rcpp::RObject X, Rcpp::RObject y, int blocks, Rcpp::Nullable<int> threads  = R_NilValue) 
{
  
  try
  {
    Eigen::MatrixXd eX;
    Eigen::MatrixXd eY;
    
    // Read DelayedArray's X   
    if ( X.isS4() == true)    
    {
      eX = read_DelayedArray(X); // Modify this function
    } else {
      try{
        eX = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
      }
      catch(std::exception &ex) { }
    }
    
    // Read DelayedArray's Y  
    if ( y.isS4() == true)    
    {
      eY = read_DelayedArray(y); // Modify this function
    } else {
      try{
        eY = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(y);
      }
      catch(std::exception &ex) { }
    }

    Eigen::MatrixXd beta = Rcpp_mlr_mr(eX, eY, blocks, threads);

    return (wrap(beta));
    
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  
  // other cases
  return (wrap(-1));
  
}


/*** R

*/
