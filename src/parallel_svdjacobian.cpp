#include "include/parallel_svdjacobian.h"
//
// https://github.com/lixueclaire/Parallel-SVD/blob/master/OMP_SVD.cpp



void svdjacob (Eigen::MatrixXd U_t, int M, int N, Eigen::MatrixXd& U, Eigen::MatrixXd& V, 
               Eigen::VectorXd& S, double &error, int &iter, Rcpp::Nullable<int> threads = R_NilValue)
{
  
  Eigen::MatrixXd V_t = Eigen::MatrixXd::Identity(N,M);
  double t, converge;
  Eigen::VectorXd C(100);
  int *I1, *I2;
  unsigned int ithreads;
  
  // int iter = 0;
  converge = 1.0;
  
  I1 = new int [N];
  I2 = new int [N];
  
  while(converge > epsilon)
  { 		
    converge = 0.0;	
    iter++;				//counter of loops
    
    for (int l = 1; l < M; l ++) 
    {
      int r1 = 0, r2 = 0;
      for (int i = 0; i + l < M; i++) 
      {
        if (i % (2 * l) < l)
        {
          r1 = r1+1;
          I1[r1] = i;
        } else {
          r2 = r2+1;
          I2[r2] = i;
        }
      }
      
      C = Eigen::VectorXd::Zero(C.size());
      
      if(threads.isNotNull()) 
      {
        if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
          ithreads = Rcpp::as<int> (threads);
        else 
          ithreads = std::thread::hardware_concurrency();
      }
      else    ithreads = std::thread::hardware_concurrency(); //omp_get_max_threads(); 
      
      omp_set_num_threads(ithreads);
      
// #pragma omp parallel for num_threads(omp_get_max_threads())
#pragma omp parallel for num_threads(ithreads)
      for (int p = 1; p <= r1; p++)
      {
        // int k = omp_get_thread_num();
        int k = omp_get_max_threads();
        int i = I1[p], j = i + l;
        double alpha = 0, beta = 0, gamma = 0;
        double zeta, t, c, s;
        for (int kl = 0; kl < N; kl++) 
        {
          alpha = alpha + U_t(i,kl)*U_t(i,kl);
          beta = beta + U_t(j,kl)*U_t(j,kl);
          gamma = gamma + (U_t(i,kl) * U_t(j,kl));
        }
        
        C[k] = std::max(C[k], std::abs(gamma)/sqrt(alpha*beta));
        
        zeta = (beta - alpha) / (2.0 * gamma);
        t = sgn(zeta) / (std::abs(zeta) + sqrt(1.0 + (zeta*zeta))); // tg of angle
        c = 1.0 / (sqrt (1.0 + (t*t)));	// cos
        s = c*t;							          // sin
        for(int kl=0; kl<N; kl++){
          t = U_t(i,kl);
          U_t(i,kl) = c*t - s*U_t(j,kl);
          U_t(j,kl) = s*t + c*U_t(j,kl);
          
          t = V_t(i,kl);
          V_t(i,kl) = c*t - s*V_t(j,kl);
          V_t(j,kl) = s*t + c*V_t(j,kl);
          
        }
      }
      
#pragma omp parallel for num_threads(ithreads)
      for (int p = 1; p <= r2; p++){
        int k = omp_get_thread_num();
        int i = I2[p], j = i + l;
        double alpha = 0, beta = 0, gamma = 0;
        double zeta, t, c, s;
        for (int kl = 0; kl < N; kl++) {
          alpha = alpha + (U_t(i,kl) * U_t(i,kl));
          beta = beta + (U_t(j,kl) * U_t(j,kl));
          gamma = gamma + (U_t(i,kl) * U_t(j,kl));
        }
        
        C[k] = std::max(C[k], std::abs(gamma)/sqrt(alpha*beta));

        zeta = (beta - alpha) / (2.0 * gamma);
        t = sgn(zeta) / (std::abs(zeta) + sqrt(1.0 + (zeta*zeta))); //tan
        c = 1.0 / (sqrt (1.0 + (t*t)));	//cos
        s = c*t;							          //sin
        for(int kl=0; kl<N; kl++){
          t = U_t(i,kl);
          U_t(i,kl) = c*t - s*U_t(j,kl);
          U_t(j,kl) = s*t + c*U_t(j,kl);
          
          t = V_t(i,kl);
          V_t(i,kl) = c*t - s*V_t(j,kl);
          V_t(j,kl) = s*t + c*V_t(j,kl);
          
        }
      }
      
      converge = std::max(converge,C.maxCoeff());
      
    }
  }
  
  //
  for(int i =0; i<M; i++)
  {
    t = sqrt(U_t.row(i).array().pow(2).sum());
    for(int j=0; j<N;j++)
    {
      U_t(i,j) = U_t(i,j) / t;
      if(i == j) {
        S[i] = t;
      }
    }
  }
  
  U = U_t.adjoint();
  V = V_t.adjoint();
  
  error = converge;
  // iter = acum;
  
}




// [[Rcpp::export]]
Rcpp::RObject JacobianSVD(Rcpp::RObject X)
{
  try
  {
    Eigen::MatrixXd A;
    
    // Read DelayedArray's X   
    if ( X.isS4() == true)    
    {
      A = read_DelayedArray(X);
    } else {
      try{  
        A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
      }
      catch(std::exception &ex) { forward_exception_to_r(ex); }
    }
    
    A.transposeInPlace();
  
    int n = A.rows(),  m = A.cols(), it = 0;
    double error; 
    
    Eigen::MatrixXd U(n,m), V(n,m);
    Eigen::VectorXd S(n);
    
    svdjacob(A, n, m, U, V, S, error, it);
    
    return Rcpp::List::create(Rcpp::Named("d") = S,
                              Rcpp::Named("u") = U,
                              Rcpp::Named("v") = V,
                              Rcpp::Named("epsilon") = error,
                              Rcpp::Named("iter") = it);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  
}



// [[Rcpp::export]]
Rcpp::RObject bdtsvd(Rcpp::RObject X, Rcpp::Nullable<int> k = R_NilValue) 
{
  
  try  
  {
    Eigen::MatrixXd eX;
    int estimateRank=20;
    
    if ( X.isS4() == true)    
    {
      eX = read_DelayedArray(X);
    } else {
      try{  
        eX = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
      }
      catch(std::exception &ex) { forward_exception_to_r(ex); }
    }
    
    // k (default): 20
    if(k.isNotNull())  estimateRank = Rcpp::as<int> (k);
    else {
      if( std::min(eX.cols(),eX.rows())>20 ) {
        estimateRank = 20;
      } else {
        estimateRank = std::min(eX.cols(),eX.rows());
      } 
    }
    
    RedSVD::RedSVD<Eigen::MatrixXd>svd(eX, estimateRank);

    return Rcpp::List::create(Rcpp::Named("d") = svd.singularValues() ,
                              Rcpp::Named("u") = svd.matrixU(),
                              Rcpp::Named("v") = svd.matrixV());

  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  
}



// [[Rcpp::export]]
Rcpp::RObject bdsvd(Rcpp::RObject X) 
{
  try
  {
    Eigen::MatrixXd eX;
    
    // Read DelayedArray's X   
    if ( X.isS4() == true)    
    {
      eX = read_DelayedArray(X);
    } else {
      try{  
        eX = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
      }
      catch(std::exception &ex) { }
    }

    if( eX.rows()<16)
    {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(eX, Eigen::ComputeThinU | Eigen::ComputeThinV);
      return Rcpp::List::create(Rcpp::Named("d") = svd.singularValues() ,
                                Rcpp::Named("u") = svd.matrixU(),
                                Rcpp::Named("v") = svd.matrixV());
      
    } else  {
      Eigen::BDCSVD<Eigen::MatrixXd> svd(eX, Eigen::ComputeThinU | Eigen::ComputeThinV);
      return Rcpp::List::create(Rcpp::Named("d") = svd.singularValues() ,
                                Rcpp::Named("u") = svd.matrixU(),
                                Rcpp::Named("v") = svd.matrixV());
    }
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  
}



/***R

library(microbenchmark)
library(DelayedArray)
A <- matrix(seq(1:16), byrow = TRUE, nrow = 4) 

U <- matrix(c(.6, .8, .8, -.6),byrow = TRUE,nrow = 2)
V <- sqrt(2)/2*(matrix(c(1,1,1,-1),byrow = TRUE,nrow = 2))
S <- diag(c(5,4))

A <- U%*%S%*%V
JacobianSVD(A)
jsvd <- JacobianSVD(A)


a <- c(1,3,2,5,6,4,7,8,9)


n <- 1500
p <- 1500
A <- matrix(rnorm(n*p), nrow=n, ncol=p)
AD <- DelayedArray(A);

results <- microbenchmark(# jsvd <- JacobianSVD(A),
                          svd <- svd(A),
                          cppsvd <- bdsvd(A),
                          delcppsvd <- bdsvd(AD),
                          cpptsvd <- bdtsvd(A, 50),
                          delcpptsvd <- bdtsvd(AD, 50),
                          times = 5L, unit='s')  # Inversa amb Cholesky
  
print(summary(results)[, c(1:7)],digits=3)

delcppsvd$d[1:10]
cppsvd$d[1:10]
svd$d[1:10]
delcpptsvd$d[1:10]
cpptsvd$d[1:10]



svd(A)



*/