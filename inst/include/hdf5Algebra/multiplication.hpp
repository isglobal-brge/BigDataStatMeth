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

            // if( ii + blockSize > maxPosition ) {
            //     isize = blockSize + 1; }
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
        
        // unsigned int ithreads;
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
                
                std::vector<hsize_t> vsizetoRead, vstart,
                                     vsizetoReadM, vstartM,
                                     vsizetoReadK, vstartK;
                
                Rcpp::Rcout<<"\nBlock Size val: "<<block_size;
                
                getBlockPositionsSizes_hdf5( N, block_size, vstart, vsizetoRead );
                getBlockPositionsSizes_hdf5( M, block_size, vstartM, vsizetoReadM );
                getBlockPositionsSizes_hdf5( K, block_size, vstartK, vsizetoReadK );
                
                int ithreads = get_number_threads(threads, R_NilValue);
                // int chunks = vstart.size()/ithreads;
                
                #pragma omp parallel num_threads(ithreads) shared(A, B, C) //..// , chunk) private(tid ) 
                {
                    
                    #pragma omp for schedule (static) // collapse(3)
                    for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                    {
                        for (hsize_t jj = 0; jj < vstartM.size(); jj++)
                        {
                            for (hsize_t kk = 0; kk < vstartK.size(); kk++)
                            {
                                C.block(vstart[ii], vstartM[jj], vsizetoRead[ii], vsizetoReadM[jj]) = 
                                    C.block(vstart[ii], vstartM[jj], vsizetoRead[ii], vsizetoReadM[jj]) + 
                                    ( A.block(vstart[ii], vstartK[kk], vsizetoRead[ii], vsizetoReadK[kk]) * 
                                      B.block(vstartK[kk], vstartM[jj], vsizetoReadK[kk], vsizetoReadM[jj]) );
                            }
                        }
                    }
                }
                
            } else {
                throw std::range_error("multiplication error: non-conformable arguments");
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcerr<< "c++ exception Bblock_matrix_mul_parallel: "<<ex.what()<< " \n";
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
            hsize_t L = dsB->ncols();
            hsize_t M = dsB->nrows();
            
             
             if( hdf5_block.isNotNull()) {
                 ihdf5_block =  Rcpp::as<int>(hdf5_block);
             } else {
                 ihdf5_block =  MAXBLOCKSIZE/3;
             }

            if( K == L )
            {

                std::vector<hsize_t> stride = {1, 1},
                                     block = {1, 1},
                                     vsizetoRead, vstart,
                                     vsizetoReadM, vstartM,
                                     vsizetoReadK, vstartK;
                
                dsC->createDataset( N, M, "real");
                
                if( dsC->getDatasetptr() != nullptr) {
                    
                    getBlockPositionsSizes_hdf5( N, ihdf5_block, vstart, vsizetoRead );
                    getBlockPositionsSizes_hdf5( M, ihdf5_block, vstartM, vsizetoReadM );
                    getBlockPositionsSizes_hdf5( K, ihdf5_block, vstartK, vsizetoReadK );
                    
                    
                    int ithreads = get_number_threads(threads, R_NilValue);
                    int chunks = vstart.size()/ithreads;
                    
                    #pragma omp parallel num_threads(ithreads) shared(dsA, dsB, dsC, chunks, vstart, vsizetoRead)
                    {
                        
                    #pragma omp for ordered 
                        for (hsize_t ii = 0; ii < vstart.size(); ii++)
                        {
                            
                            for (hsize_t jj = 0; jj < vstartM.size(); jj++)
                            {
                                
                                for (hsize_t kk = 0; kk < vstartK.size(); kk++)
                                {
                                    
                                    hsize_t iColsA = vsizetoReadK[kk],
                                            iRowsA = vsizetoRead[ii],
                                            iColsB = vsizetoReadM[jj],
                                            iRowsB = vsizetoReadK[kk];
                                    
                                    std::vector<double> vdA( iRowsA * iColsA );
                                    #pragma omp critical(accessFile)
                                    {
                                        dsA->readDatasetBlock( {vstartK[kk], vstart[ii]}, {vsizetoReadK[kk], vsizetoRead[ii]}, stride, block, vdA.data() );
                                    }
                                    
                                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), vsizetoReadK[kk], vsizetoRead[ii] );
                                    
                                    std::vector<double> vdB( iRowsB * iColsB );
                                    #pragma omp critical(accessFile)
                                    {
                                        dsB->readDatasetBlock( {vstartM[jj], vstartK[kk]}, {vsizetoReadM[jj], vsizetoReadK[kk]}, stride, block, vdB.data() );
                                    }
                                    
                                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), vsizetoReadM[jj], vsizetoReadK[kk] );
                                    
                                    std::vector<double> vdC( vsizetoReadM[jj] * vsizetoRead[ii] );
                                    #pragma omp critical(accessFile)
                                    {
                                        dsC->readDatasetBlock( {vstartM[jj], vstart[ii]}, {vsizetoReadM[jj],  vsizetoRead[ii]}, stride, block, vdC.data() );
                                    }

                                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> C (vdC.data(), B.rows(), A.cols() );

                                    C = C + B * A;

                                    std::vector<hsize_t> offset = {vstartM[jj], vstart[ii]};
                                    std::vector<hsize_t> count = {vsizetoReadM[jj], vsizetoRead[ii] };
                                    
                                    #pragma omp critical(accessFile)
                                    {
                                        dsC->writeDatasetBlock(vdC, offset, count, stride, block);
                                    }
                                }
                            }
                        }
                    }
                } 
                
            } else {
                throw std::range_error("multiplication error: non-conformable arguments");
            }

        }  catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            checkClose_file(dsA, dsB, dsC);
            Rcpp::Rcerr<<"\nc++ c++ exception multiplication (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            checkClose_file(dsA, dsB, dsC);
            Rcpp::Rcerr<<"\nc++ exception multiplication (DataSet IException)\n";
            return void();
        } catch(std::exception &ex) {
            checkClose_file(dsA, dsB, dsC);
            Rcpp::Rcerr<<"\nc++ exception multiplication\n";
            return void();
        }  catch (...) {
            checkClose_file(dsA, dsB, dsC);
            Rcpp::Rcerr<<"\nC++ exception multiplication (unknown reason)";
            return void();
        }

        return void();
    }
    
}

#endif // BIGDATASTATMETH_ALGEBRA_MULTIPLICATION_HPP