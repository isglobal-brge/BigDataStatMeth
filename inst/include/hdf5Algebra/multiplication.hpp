#ifndef BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_HPP
#define BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
// #include <thread>

namespace BigDataStatMeth {


    extern inline void multiplication( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
                                       int hdf5_block, int mem_block_size, bool bparal, Rcpp::Nullable<int> threads);




    // In-memory execution - Parallel version
    extern inline Eigen::MatrixXd Bblock_matrix_mul_parallel( Eigen::MatrixXd A, Eigen::MatrixXd B, 
                                                             int block_size, Rcpp::Nullable<int> threads  = R_NilValue)
    {
        
        unsigned int ithreads;
        Eigen::MatrixXd C;
        
        try {

            int M = A.rows();
            int K = A.cols();
            int N = B.cols();
            
            C = Eigen::MatrixXd::Zero(M,N) ;
            
            if( A.rows() == B.cols())
            {
                
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
                
    #pragma omp parallel num_threads(ithreads) shared(A, B, C) //..// , chunk) private(tid ) 
                {
                    //..// tid = omp_get_thread_num();
                    
    #pragma omp for schedule (dynamic) 
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
                
            } else {
                throw std::range_error("multiplication error: non-conformable arguments");
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception multiplication: "<<ex.what()<< " \n";
        }
        
        return(C);
        
    }
    
    
    
    
    extern inline void multiplication( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
                                       int hdf5_block, int mem_block_size, bool bparal, Rcpp::Nullable<int> threads  = R_NilValue) 
    {
        
        try {
            
            hsize_t K = dsA->nrows();
            hsize_t N = dsA->ncols();
            
            hsize_t M = dsB->nrows();
            // hsize_t L = dsB->ncols();
            
            if( dsA->nrows() == dsB->ncols())
            {
                
                hsize_t isize = hdf5_block + 1,
                    ksize = hdf5_block + 1,
                    jsize = hdf5_block + 1;
                
                std::vector<hsize_t> stride = {1, 1};
                std::vector<hsize_t> block = {1, 1};
                
                dsC->createDataset( M, N, "real"); 
                
                for (hsize_t ii = 0; ii < N; ii += hdf5_block)
                {
                    if( ii + hdf5_block > N ) isize = N - ii;
                    
                    for (hsize_t jj = 0; jj < M; jj += hdf5_block)
                    {
                        if( jj + hdf5_block > M) jsize = M - jj;
                        
                        for(hsize_t kk = 0; kk < K; kk += hdf5_block)
                        {
                            if( kk + hdf5_block > K ) ksize = K - kk;
                            
                            hsize_t iRowsA = std::min<hsize_t>(hdf5_block,ksize),
                                iColsA = std::min<hsize_t>(hdf5_block,isize),
                                iRowsB = std::min<hsize_t>(hdf5_block,jsize),
                                iColsB = std::min<hsize_t>(hdf5_block,ksize);
                            
                            std::vector<double> vdA( iRowsA * iColsA ); 
                            dsA->readDatasetBlock( {kk, ii}, {iRowsA, iColsA}, stride, block, vdA.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), iRowsA, iColsA );
                            
                            std::vector<double> vdB( iRowsB * iColsB ); 
                            dsB->readDatasetBlock( {jj, kk}, {iRowsB, iColsB}, stride, block, vdB.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), iRowsB, iColsB );
                            
                            std::vector<double> vdC( iRowsB * iColsA ); 
                            dsC->readDatasetBlock( {jj, ii}, {iRowsB, iColsA}, stride, block, vdC.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> C (vdC.data(), iRowsB, iColsA );
                            
                            if( bparal == false) {
                                C = C + B * A;
                            } else {
                                C = C + Bblock_matrix_mul_parallel(B, A, mem_block_size, threads);
                            }
                            
                            std::vector<hsize_t> offset = {jj,ii};
                            std::vector<hsize_t> count = {iRowsB, iColsA};
                            dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, true);
                            
                            if( kk + hdf5_block > K ) ksize = hdf5_block + 1;
                        }
                        
                        if( jj + hdf5_block > M ) jsize = hdf5_block + 1;
                    }
                    
                    if( ii + hdf5_block > N ) isize = hdf5_block + 1;
                }
                
            }else {
                throw std::range_error("multiplication error: non-conformable arguments");
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception multiplication: "<<ex.what()<< " \n";
            // return(dsC);
            return void();
            
        }
        
        return void();
        // return(dsC);
    }
    
    

}

#endif // BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_HPP