#ifndef BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_SPARSE_HPP
#define BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_SPARSE_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
// #include <thread>

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

}

#endif // BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_SPARSE_HPP