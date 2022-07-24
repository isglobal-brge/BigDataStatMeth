#include "include/parallelCholesky.h"




Eigen::MatrixXd Cholesky_decomposition_parallel( Eigen::MatrixXd& A, Rcpp::Nullable<int> threads = R_NilValue )
{
    
    int dimensionSize = A.rows();
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(dimensionSize,dimensionSize);
    L.triangularView<Eigen::Lower>() = A.triangularView<Eigen::Lower>();
    int i,k,chunk = 1; 
    // double sum = 0;
    unsigned int ithreads;
    
    if(threads.isNotNull()) {
        if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
            ithreads = Rcpp::as<int> (threads);
        } else {
            ithreads = getDTthreads(0, true);
        }
    } else {
        ithreads = getDTthreads(0, true);
    }
    
    
    if( A.rows()== A.cols() )
    {
        for (int j = 0; j < dimensionSize; j++)  {
            
            if(j==0) {
                L(j,j) = std::sqrt(A(j,j));
            } else {
                L(j,j) = std::sqrt(A(j,j) - (L.row(j).head(j).array().pow(2).sum() ));
            }
            
#pragma omp parallel for num_threads(getDTthreads(ithreads, true)) shared (A,L,j) schedule(static)
            for (int i = j + 1; i < dimensionSize; i++) 
            {
                double sum = 0;
                if( j > 0) {
                    sum = (L.block(i, 0, 1, j).array() * L.block(j, 0, 1, j).array()).array().sum();
                } else {
                    sum = 0;
                }
                L(i,j) =  (1/L(j,j)*(A(i,j) - sum));
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
  
  if(threads.isNotNull()) {
      if (Rcpp::as<unsigned int> (threads) <= std::thread::hardware_concurrency()) {
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
  
  if( L.rows()== L.cols() )
  {
    for( i=1; i<n; i++)
//.OpenMP.// #pragma omp parallel shared(L,y) private(j)
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(L,y) private(j)
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






Eigen::MatrixXd Inverse_of_Cholesky_decomposition_parallel(  int dimensionSize, Eigen::MatrixXd& L, Rcpp::Nullable<int> threads = R_NilValue )
{
    
    Eigen::MatrixXd TT = L.triangularView<Eigen::StrictlyLower>();
    TT.triangularView<Eigen::StrictlyUpper>() =  Eigen::MatrixXd::Constant(dimensionSize,dimensionSize,0).triangularView<Eigen::StrictlyUpper>();
    TT.diagonal() = L.diagonal().cwiseInverse();
    
    // Rcpp::Rcout<<" Matriu TT 1: \n"<<TT<<"\n";
    // Rcpp::Rcout<<" Matriu TT 2: \n"<<TT<<"\n";
    // Rcpp::Rcout<<" Diagonal TT: \n"<<TT.diagonal()<<"\n";
    
    int j;// , k; 
    unsigned int ithreads;
    
    if(threads.isNotNull()) {
        if (Rcpp::as<unsigned int> (threads) <= std::thread::hardware_concurrency()){
            ithreads = Rcpp::as<int> (threads);
        } else {
            ithreads = getDTthreads(0, true);
        }
    } else {
        ithreads = getDTthreads(0, true);
    }
    
    if( L.rows() == L.cols() ) {

        for (int j = 1; j < dimensionSize; j++) 
        {
            Eigen::VectorXd vR(j);
            Eigen::ArrayXd ar_j = TT.block(j, 0, 1, j).transpose().array();
            // Rcpp::Rcout<<"\nar_j:\n"<<ar_j<<"\n";
#pragma omp parallel for num_threads(getDTthreads(ithreads, true)) shared (j, TT, vR, ar_j) schedule(static)
            for (int i = 0; i < j ; i++) {
                Eigen::ArrayXd ar_i = TT.block(0, i, j, 1).array();
                // Rcpp::Rcout<<"\nar_i:\n"<<ar_i<<"\n";
                vR(i) = ((ar_j.transpose() * ar_i.transpose()) * (-1)).sum() / L.diagonal()[j];
                // Rcpp::Rcout<<"\nvR(i):\n"<<vR<<"\n";
            }
            
            TT.block(j, 0, 1, j) = vR.transpose();
            // Rcpp::Rcout<<"\nTT:\n"<<TT<<"\n";
            
        }
        
        return(TT);
        
    } else {
        throw std::range_error("non-conformable arguments");
    }
    
}





Eigen::MatrixXd Inverse_Matrix_Cholesky_parallel( Eigen::MatrixXd InvCh, Rcpp::Nullable<int> threads = R_NilValue  )
{
    int i, j;
    int dimensionSize = InvCh.rows();
    unsigned int ithreads;
    Eigen::VectorXd newDiag(dimensionSize);
    
    //..// Eigen::MatrixXd InvCh = Eigen::MatrixXd::Zero(dimensionSize,dimensionSize);
    InvCh = InvCh.triangularView<Eigen::Lower>();
    
    if(threads.isNotNull()) {
        if (Rcpp::as<unsigned int> (threads) <= std::thread::hardware_concurrency()){
            ithreads = Rcpp::as<int> (threads);
        } else {
            ithreads = getDTthreads(0, true);
        }
    } else {
        ithreads = getDTthreads(0, true);
    }
    
    if( InvCh.rows()== InvCh.cols() )
    {
        
        for (i = 0; i < dimensionSize; i++) 
        {
            newDiag(i) = InvCh.block(i, i, dimensionSize-i, 1).array().pow(2).sum();

#pragma omp parallel for num_threads(getDTthreads(ithreads, true)) private(j) shared (InvCh,i,dimensionSize) schedule(static)
            for (j = i + 1; j < dimensionSize; j++)  {
                InvCh(i,j) = (InvCh.block(j, i, dimensionSize-j, 1).array() * InvCh.block(j, j, dimensionSize-j, 1).array()).sum();
            }
        }
        
        InvCh.diagonal() = newDiag;
        
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
    
    try {  
        if ( TYPEOF(a) == INTSXP ) {
            A = Rcpp::as<Eigen::MatrixXi>(a).cast<double>();
        } else {
            A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(a);
        }
    }
    catch(std::exception &ex) { }
    
    
    if( A.rows() == A.cols() )   {
        Eigen::LLT<Eigen::MatrixXd> lltOfA(A);
        if(lltOfA.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("Possibly non semi-positive definitie matrix!");
            
        } else {
            Eigen::MatrixXd L = Cholesky_decomposition_parallel(A, threads);
            // Rcpp::Rcout<<"\n La matriu L val : \n"<<L<<"\n";
            
            Eigen::MatrixXd Z = Inverse_of_Cholesky_decomposition_parallel( A.rows(), L, threads);
            
            // Rcpp::Rcout<<"\n La matriu Z val : \n"<<Z<<"\n";
            
            return(  Inverse_Matrix_Cholesky_parallel(Z));  
        }
        
    } else   {
        
        throw std::range_error("non-conformable arguments");
        
    }
    
}




/*** R

library(BigDataStatMeth)
library(microbenchmark)



A  <- matrix(sample.int(10, 25, replace = TRUE), ncol = 5)
A <- crossprod(A)
# b <- diag(100)

F <- BigDataStatMeth:::inversechol_par(A)


E <- BigDataStatMeth:::inversechol_par(A)

microbenchmark::microbenchmark( T <- solve(A),
                                E <- BigDataStatMeth:::inversechol_par(A),
                                F <- BigDataStatMeth:::inversechol_par(A)
                              )







A = as.matrix(data.frame(c(3,4,3),c(4,8,6),c(3,6,9)))
colnames(A) <- NULL


A  <- matrix(sample.int(1000, 1000, replace = TRUE), ncol = 100)
A <- crossprod(A)
b <- diag(100)

microbenchmark::microbenchmark( T <- solve(A),
                                C <- BigDataStatMeth::bdInvCholesky(A),
                                E <- BigDataStatMeth:::inversechol_par(A,threads = 2),
                                F <- BigDataStatMeth:::inversechol_par(A,threads = 3),
                                G <- BigDataStatMeth:::inversechol_par(A,threads = 4),
                                ST <- BigDataStatMeth:::CholSolve(A,b),
                                bdSolve(A, b),
                                solve(A, b)
                                )


chol(A)
solve(A)[1:5,1:5]
solve(A, b)[1:5,1:5]

bdSolve(A, A)[1:5,1:5]
solve(A, A)[1:5,1:5]


solve(A)

# Equivalents utilitzant BigDataStatMeth 
BigDataStatMeth:::inversechol_par(A)
BigDataStatMeth::bdInvCholesky(A)

# BigDataStatMeth::bdInvCholesky(A)

BigDataStatMeth:::inversechol_par(A)









library(BigDataStatMeth)
library(microbenchmark)


set.seed(1234)
A  <- matrix(sample.int(10, 10000000, replace = TRUE), ncol = 10000)
A <- crossprod(A)
# A
# b <- diag(100)

microbenchmark::microbenchmark( T <- solve(A),
                                F <- BigDataStatMeth:::inversechol_par(A),
                                G2 <- BigDataStatMeth:::inversechol_par(A,threads = 4),
                                # Z <- bdInvCholesky(A), 
                                times = 1
    )
all.equal(T,G2)
all.equal(T,F)





*/
