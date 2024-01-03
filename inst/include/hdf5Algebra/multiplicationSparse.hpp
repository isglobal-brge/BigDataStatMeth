#ifndef BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_SPARSE_HPP
#define BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_SPARSE_HPP

#include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include <thread>

namespace BigDataStatMeth {


// Eigen not allow blocks 
// extern inline Eigen::SparseMatrix<double> Bblock_matrix_mul_parallel_sparse(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B,
//                                                                             int block_size, Rcpp::Nullable<int> threads);

extern inline BigDataStatMeth::hdf5Dataset* multiplication( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
                                                            hsize_t hdf5_block, hsize_t mem_block_size, bool bparal, bool browmajor, Rcpp::Nullable<int> threads);

extern inline BigDataStatMeth::hdf5Dataset* multiplicationSparse( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
                                                                  hsize_t hdf5_block, hsize_t mem_block_size, bool bparal, bool browmajor, Rcpp::Nullable<int> threads  = R_NilValue) 
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
                    // Això haurien de ser files i no per columnes
                    for (hsize_t jj = 0; jj < M; jj += hdf5_block)
                    {
                        
                        if( jj + hdf5_block > M) jsize = M - jj;
                        
                        for(hsize_t kk = 0; kk < K; kk += hdf5_block)
                        {
                            
                            if( kk + hdf5_block > K ) ksize = K - kk;
                            
                            
                            hsize_t iRowsA = std::min(hdf5_block,ksize),
                                    iColsA = std::min(hdf5_block,isize),
                                    iRowsB = std::min(hdf5_block,jsize),
                                    iColsB = std::min(hdf5_block,ksize);
                            
                            std::vector<double> vdA( iRowsA * iColsA ); 
                            dsA->readDatasetBlock( {kk, ii}, {iRowsA, iColsA}, stride, block, vdA.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), iRowsA, iColsA );
                            
                            std::vector<double> vdB( iRowsB * iColsB ); 
                            dsB->readDatasetBlock( {jj, kk}, {iRowsB, iColsB}, stride, block, vdB.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), iRowsB, iColsB );
                            
                            std::vector<double> vdC( iRowsB * iColsA ); 
                            dsC->readDatasetBlock( {jj, ii}, {iRowsB, iColsA}, stride, block, vdC.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> C (vdC.data(), iRowsB, iColsA );
                            
                            // if( bparal == false) {
                            //     C = C + B * A;
                            // } else {
                            //     C = C + Bblock_matrix_mul_parallel(B, A, mem_block_size, threads);
                            // }
                            
                            C = C + (B.sparseView() * A.sparseView()).toDense() ;
                            
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
                throw std::range_error("multiplicationSparse error: non-conformable arguments");
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception multiplicationSparse: "<<ex.what()<< " \n";
            return(dsC);
        }
        
        return(dsC);

    }




//  ----------------------------
//      FALTA REPASSAR !!!!
//  ----------------------------
//
// // In-memory execution - Parallel version 
// 
// 
// extern inline Eigen::SparseMatrix<double> Bblock_matrix_mul_parallel_sparse(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B,
//                                                          int block_size, Rcpp::Nullable<int> threads  = R_NilValue)
// {
//     //..// int ii=0, jj=0, kk=0;
//     //..// int chunk = 1, tid;
//     unsigned int ithreads;
//     int M = A.rows();
//     int K = A.cols();
//     int N = B.cols();
// 
//     //..// Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N) ;
// 
//     // Eigen::SparseMatrix<double> C(M, N);
//     Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,N);
// 
//     if(block_size > std::min( N, std::min(M,K)) )
//         block_size = std::min( N, std::min(M,K));
// 
//     if(threads.isNotNull()) {
//         if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
//             ithreads = Rcpp::as<int> (threads);
//         } else {
//             ithreads = getDTthreads(0, true);
//         }
//     } else {
//         ithreads = getDTthreads(0, true);
//     }
// 
//     #pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C) //..// , chunk) private(tid )
//     {
//             //..// tid = omp_get_thread_num();
// 
//     #pragma omp for schedule (dynamic)
//             // Sparse matrix only works with continuous and complete rows or cols (rowmajor/colmajor) 
//             // not by (x,y, sizex, sizey)
//             for (int ii = 0; ii < N; ii += block_size)
//             {
//                  int cur_blockSize = std::min(block_size, N - ii);
//                 
//                 // C.block( ii, 0, cur_blockSize, M) = (C.block( ii, 0, cur_blockSize, M).sparseView() +
//                 //     ( A.block(ii, 0, cur_blockSize, M) * B.block(ii, 0, cur_blockSize, M) ) ).toDense();
//                     
//                 Rcpp::Rcout<<"Matriu C - Sparse View:\n"<<C.sparseView().middleRows( ii, cur_blockSize)<<"\n";
//                 Rcpp::Rcout<<"Matriu C - Dense View:\n"<<C.block( ii, 0, cur_blockSize, M)<<"\n";
//                 Rcpp::Rcout<<"Matriu A:\n"<<A.middleRows( ii, cur_blockSize)<<"\n";
//                 Rcpp::Rcout<<"Matriu B:\n"<<B.middleRows( ii, cur_blockSize)<<"\n";
//                 
//                 C.block( ii, 0, cur_blockSize, M) = (C.block( ii, 0, cur_blockSize, M).sparseView() +
//                 ( A.middleRows( ii, cur_blockSize) * B.middleCols( ii, cur_blockSize) ) ).toDense();
//                 
//             }
//     }
//     return( C.sparseView() );
// }

}

#endif // BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_SPARSE_HPP