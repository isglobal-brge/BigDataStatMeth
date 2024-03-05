#ifndef BIGDATASTATMETH_ALGEBRA_MEM_MULTIPLICATION_HPP
#define BIGDATASTATMETH_ALGEBRA_MEM_MULTIPLICATION_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
// #include <thread>

namespace BigDataStatMeth {

extern inline Eigen::MatrixXd block_matrix_mul( Eigen::MatrixXd& A, Eigen::MatrixXd& B, int block_size);
// extern inline Eigen::MatrixXd block_matrix_mul_parallel( Eigen::MatrixXd& A, Eigen::MatrixXd& B, int block_size, Rcpp::Nullable<int> threads);

// In-memory execution - Serial version
extern inline Eigen::MatrixXd block_matrix_mul( Eigen::MatrixXd& A, Eigen::MatrixXd& B, int block_size)
{
    
    Eigen::MatrixXd C;
    
    try{
        
        int M = A.rows();
        int K = A.cols();
        int N = B.cols();
        
        if( K == B.rows())
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
    template<typename T>
    extern inline Eigen::MatrixXd block_matrix_mul_parallel( T X, T Y, 
                                               int block_size, Rcpp::Nullable<int> threads  = R_NilValue)
    {
        
        try{
            static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
                          std::is_same<T, Eigen::Map< Eigen::MatrixXd >>::value || 
                          std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                          std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
                          "Error - type not allowed");
            
            Eigen::MatrixXd A = X,
                            B = Y;
            
            //..// int ii=0, jj=0, kk=0;
            int chunk = 1, tid;
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
                }
            } else {
                ithreads = getDTthreads(0, true);
            }
            
            #pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C, chunk) private(tid ) 
            {
                
                tid = omp_get_thread_num();
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
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception multiplication: "<<ex.what()<< " \n";
        }
    
    }

}

#endif // BIGDATASTATMETH_ALGEBRA_MEM_MULTIPLICATION_HPP