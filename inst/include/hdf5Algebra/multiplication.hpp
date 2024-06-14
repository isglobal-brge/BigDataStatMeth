#ifndef BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_HPP
#define BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include "memAlgebra/memMultiplication.hpp"
// #include <thread>

namespace BigDataStatMeth {


    extern inline void multiplication( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
                                       Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> hdf5_block, Rcpp::Nullable<int> threads);




extern inline void getBlockPositionsSizes_hdf5( hsize_t maxPosition, hsize_t blockSize, std::vector<hsize_t>& starts, std::vector<hsize_t>& sizes ){

        hsize_t isize = blockSize + 1;

        for (hsize_t ii = 0; ii < maxPosition; ii += blockSize)
        {
            if( ii + blockSize > maxPosition ) {
                isize = maxPosition - ii; }

            hsize_t sizetoRead = getOptimBlockSize( maxPosition, blockSize, ii, isize);

            starts.push_back(ii);
            sizes.push_back(sizetoRead);

            if( ii + blockSize > maxPosition ) isize = blockSize + 1;
            if( sizetoRead > blockSize ) {
                ii = ii - blockSize + sizetoRead; }
        }

    }


    // In-memory execution - Parallel version
    // 
    //  IMPORTANT : FUNCIÓ MODIFICADA EL 2024/04/06  I NO TESTEJADA !!!!
    // 
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
                
                std::vector<hsize_t> vsizetoRead, vstart;
                
                // ithreads = get_number_threads(threads, R_NilValue);
                
                getBlockPositionsSizes_hdf5( N, block_size, vstart, vsizetoRead );
                
                int ithreads = get_number_threads(threads, R_NilValue);
                int chunks = vstart.size()/ithreads;
                
                #pragma omp parallel num_threads(ithreads) shared(A, B, C) //..// , chunk) private(tid ) 
                {
                    //..// tid = omp_get_thread_num();
                    
                    // #pragma omp for schedule (dynamic) 
                    // for (int ii = 0; ii < M; ii += block_size)
                    #pragma omp for schedule (dynamic, chunks) // collapse(3)
                    for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                    {
                        // Rcpp::Rcout << "Number of threads: " << omp_get_num_threads() << "\n";
                        for (int jj = 0; jj < N; jj += block_size)
                        {
                            for(int kk = 0; kk < K; kk += block_size)
                            {
                                // C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) = 
                                //     C.block(ii, jj, std::min(block_size,M - ii), std::min(block_size,N - jj)) + 
                                //     (A.block(ii, kk, std::min(block_size,M - ii), std::min(block_size,K - kk)) * 
                                //     B.block(kk, jj, std::min(block_size,K - kk), std::min(block_size,N - jj)));
                                
                                C.block(vstart[ii], jj, vsizetoRead[ii], std::min(block_size,N - jj)) = 
                                    C.block(vstart[ii], jj, vsizetoRead[ii], std::min(block_size,N - jj)) + 
                                    (A.block(vstart[ii], kk, vsizetoRead[ii], std::min(block_size,K - kk)) * 
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
                                       Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> hdf5_block, Rcpp::Nullable<int> threads  = R_NilValue) 
    {
        
        try {
        
            int ihdf5_block;
            hsize_t K = dsA->nrows();
            hsize_t N = dsA->ncols();

            hsize_t M = dsB->nrows();
            // hsize_t L = dsB->ncols();
             
             if( hdf5_block.isNotNull()) {
                 ihdf5_block =  Rcpp::as<int>(hdf5_block);
             } else {
                 ihdf5_block =  MAXBLOCKSIZE/3;
             }

            if( dsA->nrows() == dsB->ncols())
            {

                hsize_t ksize = ihdf5_block + 1,
                        jsize = ihdf5_block + 1;

                std::vector<hsize_t> stride = {1, 1},
                                     block = {1, 1},
                                     vsizetoRead, vstart;
                
                dsC->createDataset( M, N, "real");
                
                getBlockPositionsSizes_hdf5( N, ihdf5_block, vstart, vsizetoRead );
                
                int ithreads = get_number_threads(threads, R_NilValue);
                int chunks = vstart.size()/ithreads;

                #pragma omp parallel num_threads(ithreads) shared(dsA, dsB, dsC, chunks, vstart, vsizetoRead)
                {
                    
                    #pragma omp for schedule (dynamic, chunks)
                    for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                        // for (hsize_t ii = 0; ii < N; ii += ihdf5_block)
                    {
                        // if( ii + ihdf5_block > N ) isize = N - ii;
                        
                        for (hsize_t jj = 0; jj < M; jj += ihdf5_block)
                        {
                            if( jj + ihdf5_block > M) jsize = M - jj;
                            
                            for(hsize_t kk = 0; kk < K; kk += ihdf5_block)
                            {
                                if( kk + ihdf5_block > K ) ksize = K - kk;
                                
                                hsize_t iRowsA = std::min<hsize_t>(ihdf5_block, ksize),
                                    iColsA = std::min<hsize_t>(ihdf5_block, vsizetoRead[ii]),
                                    iRowsB = std::min<hsize_t>(ihdf5_block, jsize),
                                    iColsB = std::min<hsize_t>(ihdf5_block, ksize);
                                
                                std::vector<double> vdA( iRowsA * iColsA );
                                #pragma omp critical(accessFile) 
                                {
                                    Rcpp::Rcout<<"\nLlegint dsA (inici - Fi) + (files - columnes): ( "<< kk << " - "<<vstart[ii]<<" ) + ("<<iRowsA<<" - "<<iColsA<<" )";
                                    dsA->readDatasetBlock( {kk, vstart[ii]}, {iRowsA, iColsA}, stride, block, vdA.data() );
                                }
                                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), iRowsA, iColsA );
                                
                                std::vector<double> vdB( iRowsB * iColsB );
                                #pragma omp critical(accessFile) 
                                {
                                    Rcpp::Rcout<<"\nLlegint dsB (inici - Fi) + (files - columnes): ( "<< jj << " - "<<kk<<" ) + ("<<iRowsB<<" - "<<iColsB<<" )";
                                    dsB->readDatasetBlock( {jj, kk}, {iRowsB, iColsB}, stride, block, vdB.data() );
                                }
                                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), iRowsB, iColsB );
                                
                                std::vector<double> vdC( iRowsB * iColsA );
                                #pragma omp critical(accessFile) 
                                {
                                    Rcpp::Rcout<<"\nLlegint dsC (inici - Fi) + (files - columnes): ( "<< jj << " - "<<vstart[ii]<<" ) + ("<<iRowsB<<" - "<<iColsA<<" )";
                                    dsC->readDatasetBlock( {jj, vstart[ii]}, {iRowsB, iColsA}, stride, block, vdC.data() );
                                }
                                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> C (vdC.data(), iRowsB, iColsA );
                                
                                Rcpp::Rcout<<"\nMatriu A:\n"<<A;
                                Rcpp::Rcout<<"\nMatriu B:\n"<<B;
                                Rcpp::Rcout<<"\nMatriu C:\n"<<C;
                                
                                C = C + B * A;
                                
                                Rcpp::Rcout<<"\nMatriu C-Resultant:"<<C;
                                
                                std::vector<hsize_t> offset = {jj,vstart[ii]};
                                std::vector<hsize_t> count = {iRowsB, iColsA};
                                
                                #pragma omp critical(accessFile) 
                                {
                                    dsC->writeDatasetBlock(vdC, offset, count, stride, block);
                                }
                                
                                if( kk + ihdf5_block > K ) ksize = ihdf5_block + 1;
                            }
                            
                            if( jj + ihdf5_block > M ) jsize = ihdf5_block + 1;
                        }
                        
                    }
                    
                }

            } else {
                throw std::range_error("multiplication error: non-conformable arguments");
            }

        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception multiplication: "<<ex.what()<< " \n";
            return void();
        }

        return void();
    }
    
}

#endif // BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_HPP