#include "include/parallelCholesky.h"


Eigen::MatrixXd Cholesky_decomposition_parallel( Eigen::MatrixXd& A, Rcpp::Nullable<int> threads = R_NilValue )
{
  
  int dimensionSize = A.rows();
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(dimensionSize,dimensionSize);
  L.triangularView<Eigen::Lower>() = A.triangularView<Eigen::Lower>();
  int i,k,chunk = 1; 
  double sum = 0;
  unsigned int ithreads;
  
  // int dimensionSize = A.rows();
  if( A.rows()== A.cols() )
  {
    for (int j = 0; j < dimensionSize; j++) 
    {
      if(j==0)
        L(j,j) = std::sqrt(A(j,j));
      else 
        L(j,j) = std::sqrt(A(j,j) - (L.row(j).head(j).array().pow(2).sum() ));
      
      if(threads.isNotNull()) 
      {
        if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
          ithreads = Rcpp::as<int> (threads);
        else 
          ithreads = std::thread::hardware_concurrency() /2;
      }
      else    ithreads = std::thread::hardware_concurrency()/2; //omp_get_max_threads();
      
      omp_set_num_threads(ithreads);
      
#pragma omp parallel for private(i,k,sum) shared (A,L,j) schedule(static) if (j < dimensionSize - chunk)
      for (int i = j + 1; i < dimensionSize; i++) 
      {
        sum = 0;
        for (int k = 0; k < j; k++) {
          sum = sum + L(i,k) * L(j,k);
        }
        L(i,j) = ( 1/L(j,j)*(A(i,j)- sum));
      }
    }
    
    L.triangularView<Eigen::StrictlyUpper>()=L.adjoint().triangularView<Eigen::StrictlyUpper>();
    
    return(L);
  } else {
    throw std::range_error("non-conformable arguments");
  }
  
}



Eigen::VectorXd Forward_Substituion_parallel(Eigen::MatrixXd L, Eigen::VectorXd y, Rcpp::Nullable<int> threads = R_NilValue)
{
  
  int i, j;
  int n = L.cols();
  unsigned int ithreads;
  
  if(threads.isNotNull()) 
  {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
      ithreads = Rcpp::as<int> (threads);
    else 
      ithreads = std::thread::hardware_concurrency()/2;
  }
  else    ithreads = std::thread::hardware_concurrency() /2; //omp_get_max_threads();
  
  omp_set_num_threads(ithreads);
  
  if( L.rows()== L.cols() )
  {
    for( i=1; i<n; i++)
#pragma omp parallel shared(L,y) private(j)
{
  for(j=i; j<n; j++)
  {
    y[j] = y[j] - L(j,i-1)*y[i-1];
  }
}
    return(y);
  } else {
    throw std::range_error("non-conformable arguments");
  }
}

/*
Eigen::MatrixXd Backward_Substitution()
{

}
*/


Eigen::MatrixXd Inverse_of_Cholesky_decomposition_parallel( Eigen::MatrixXd& A, Eigen::MatrixXd& L, Rcpp::Nullable<int> threads = R_NilValue )
{
  
  int dimensionSize = A.rows();
  Eigen::MatrixXd R = L.triangularView<Eigen::StrictlyLower>();
  R.triangularView<Eigen::StrictlyUpper>() =  Eigen::MatrixXd::Constant(dimensionSize,dimensionSize,-1).triangularView<Eigen::StrictlyUpper>();
  R.diagonal() = A.diagonal();
  
  int j, k; 
  double sum = 0;
  unsigned int ithreads;
  
  if(threads.isNotNull()) 
  {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
      ithreads = Rcpp::as<int> (threads);
    else 
      ithreads = std::thread::hardware_concurrency()/2;
  }
  else    ithreads = std::thread::hardware_concurrency() /2; //omp_get_max_threads();
  
  omp_set_num_threads(ithreads);
  
  if( L.rows()== L.cols() )
  {
    for (int i = 0; i < dimensionSize; i++) 
    {
      R(i,i) = 1/L(i,i);
#pragma omp parallel for private(j,k,sum) shared (L,R,A,i,dimensionSize) schedule(static) if (j < dimensionSize)
      for (int j = i + 1; j < dimensionSize; j++) 
      {
        sum = 0;
        for (int k = i; k < j; k++) {
          sum = sum - R(j,k) * R(k,i);
        }
        R(j,i) = sum / L.diagonal()[j];
      }
    }
    return(R);
    
  } else {
    throw std::range_error("non-conformable arguments");
  }
  
}




Eigen::MatrixXd Inverse_Matrix_Cholesky_parallel( Eigen::MatrixXd L, Rcpp::Nullable<int> threads = R_NilValue  )
{
  int i, j, k;
  int dimensionSize = L.rows();
  unsigned int ithreads;
  
  Eigen::MatrixXd InvCh = Eigen::MatrixXd::Zero(dimensionSize,dimensionSize);
  InvCh = L.triangularView<Eigen::Lower>();
  
  if(threads.isNotNull()) 
  {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
      ithreads = Rcpp::as<int> (threads);
    else 
      ithreads = std::thread::hardware_concurrency()/2;
  }
  else    ithreads = std::thread::hardware_concurrency() /2; //omp_get_max_threads();
  
  omp_set_num_threads(ithreads);
  
  if( L.rows()== L.cols() )
  {
    
    for (i = 0; i < dimensionSize; i++) 
    {
      InvCh(i,i) = std::pow(InvCh(i,i),2);
#pragma omp parallel for private(j,k) shared (InvCh,i,L,dimensionSize) schedule(static) if (j < dimensionSize)
      for (k = i + 1; k < dimensionSize; k++) {
        InvCh(i,i) += InvCh(k,i) * InvCh(k,i);
      }
      for (j = i + 1; j < dimensionSize; j++) 
      {
        for (k = j; k < dimensionSize; k++) 
        {
          InvCh(i,j) += InvCh(k,i) * InvCh(k,j);
        }
      }
    }
    
    InvCh.triangularView<Eigen::StrictlyLower>()=InvCh.adjoint();
    return(InvCh);
    
  } else {
    throw std::range_error("non-conformable arguments");
  }
  
}




// [[Rcpp::export]]
Eigen::MatrixXd CholFactor(Rcpp::RObject a)
{
  Eigen::MatrixXd A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
  
  if( A.rows()== A.cols() ) {
    
    return(Cholesky_decomposition_parallel(A));
    
  } else {
    
    throw std::range_error("non-conformable arguments");
  }
}


// [[Rcpp::export]]
Eigen::MatrixXd CholSolve(Rcpp::RObject a, Rcpp::RObject b)
{
  
  Eigen::MatrixXd A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
  Eigen::VectorXd B = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(b);
  
  if( A.rows()== A.cols() )   {
    
    Eigen::MatrixXd L = Cholesky_decomposition_parallel(A);
    return(Forward_Substituion_parallel( L.triangularView<Eigen::Lower>(), B));
    
  } else   {
    
    throw std::range_error("non-conformable arguments");
    
  }
  
}


// [[Rcpp::export]]
Eigen::MatrixXd inversechol_par(Rcpp::RObject a, Rcpp::Nullable<int> threads = R_NilValue)
{
  
  Eigen::MatrixXd A; // = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
  
  if ( a.isS4() == true)    
  {
    A = read_DelayedArray(a);
  } else {
    try{  
      A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
    }
    catch(std::exception &ex) { }
  }
  
  
  if( A.rows()== A.cols() )   {
    Eigen::LLT<Eigen::MatrixXd> lltOfA(A);
    if(lltOfA.info() == Eigen::NumericalIssue)
    {
      throw std::runtime_error("Possibly non semi-positive definitie matrix!");
    } else
    {
      Eigen::MatrixXd L = Cholesky_decomposition_parallel(A);
      Eigen::MatrixXd Z = Inverse_of_Cholesky_decomposition_parallel(A,L);
      return(  Inverse_Matrix_Cholesky_parallel(Inverse_of_Cholesky_decomposition_parallel(A,L, threads)));  
    }

  } else   {
    
    throw std::range_error("non-conformable arguments");
    
  }
  
}





/*** R

*/