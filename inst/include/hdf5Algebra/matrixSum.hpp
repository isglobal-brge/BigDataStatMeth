#ifndef BIGDATASTATMETH_ALGEBRA_SUM_HPP
#define BIGDATASTATMETH_ALGEBRA_SUM_HPP

#include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include <thread>
#include "memAlgebra/memSum.hpp"

namespace BigDataStatMeth {


extern inline BigDataStatMeth::hdf5Dataset*  Rcpp_block_matrix_sum_hdf5( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
        hsize_t hdf5_block, hsize_t mem_block_size, bool bparal, Rcpp::Nullable<int> threads);

extern inline BigDataStatMeth::hdf5Dataset* hdf5_block_matrix_vector_sum_hdf5_transposed( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
        hsize_t hdf5_block, bool bparal, Rcpp::Nullable<int> threads);


// Working directly with C-matrix in hdf5 file
extern inline BigDataStatMeth::hdf5Dataset*  Rcpp_block_matrix_sum_hdf5( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
        hsize_t hdf5_block, bool bparal, Rcpp::Nullable<int> threads  = R_NilValue)
{
    
    try {
        
        unsigned int ithreads;
        hsize_t K = dsA->nrows();
        hsize_t N = dsA->ncols();
        
        if( K == dsB->nrows() && N == dsB->ncols())
        {
            
            hsize_t isize = hdf5_block + 1;
            
            std::vector<hsize_t> stride = {1, 1};
            std::vector<hsize_t> block = {1, 1};
            
            dsC->createDataset( N, K, "real"); 
            
            if(threads.isNotNull()) {
                if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
                    ithreads = Rcpp::as<int> (threads);
                } else {
                    ithreads = getDTthreads(0, true);
                }
            } else {
                ithreads = getDTthreads(0, true);
            }
            
            if( K<=N ) {
                
                #pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(dsA, dsB, dsC)
                {
                #pragma omp for schedule (static)
                    for (hsize_t ii = 0; ii < N; ii += hdf5_block)
                    {
                        
                        if( ii + hdf5_block > N ) {
                            isize = N - ii; }
                        
                        hsize_t sizetoRead = getOptimBlockSize( N, hdf5_block, ii, isize);
                        
                        std::vector<double> vdA( K * sizetoRead ); 
                        dsA->readDatasetBlock( {0, ii}, { K, sizetoRead}, stride, block, vdA.data() );
                        
                        std::vector<double> vdB( K * sizetoRead ); 
                        dsB->readDatasetBlock( {0, ii}, {K, sizetoRead}, stride, block, vdB.data() );
                        
                        std::transform (vdA.begin(), vdA.end(),
                                        vdB.begin(), vdA.begin(), std::plus<double>());
                        
                        std::vector<hsize_t> offset = { 0, ii };
                        std::vector<hsize_t> count = { K, sizetoRead };
                        #pragma omp critical 
                        {
                            dsC->writeDatasetBlock(vdA, offset, count, stride, block, true);
                        }
                        
                        if( ii + hdf5_block > N ) isize = hdf5_block + 1;
                        if( sizetoRead > hdf5_block ) {
                            ii = ii - hdf5_block + sizetoRead; }
                    }
                }
                
            } else {
                
                #pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(dsA, dsB, dsC)
                {
                #pragma omp for schedule (static)
                    for (hsize_t ii = 0; ii < K; ii += hdf5_block)
                    {
                        
                        if( ii + hdf5_block > K ) 
                            isize = K - ii;
                        
                        hsize_t sizetoRead = getOptimBlockSize( K, hdf5_block, ii, isize);
                        
                        std::vector<double> vdA( sizetoRead * N ); 
                        dsA->readDatasetBlock( {ii, 0}, { sizetoRead, N}, stride, block, vdA.data() );
                        
                        std::vector<double> vdB( sizetoRead * N); 
                        dsB->readDatasetBlock( {ii, 0}, {sizetoRead, N}, stride, block, vdB.data() );
                        
                        std::transform (vdA.begin(), vdA.end(),
                                        vdB.begin(), vdA.begin(), std::plus<double>());
                        
                        std::vector<hsize_t> offset = { ii, 0 };
                        std::vector<hsize_t> count = { N, sizetoRead };
                        #pragma omp critical 
                        {
                            dsC->writeDatasetBlock(vdA, offset, count, stride, block, true);
                        }
                        
                        if( ii + hdf5_block > K ) isize = hdf5_block + 1;
                        if( sizetoRead > hdf5_block ) {
                            ii = ii - hdf5_block + sizetoRead; }
                    }
                }
            }
        } else {
            Rcpp::Rcout<<"matrix sum error: non-conformable arguments\n";
        }
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return(dsC);
    }
    
    return(dsC);
}



// This function makes summatori with a matrix and a vector, returns a matrix.
//
// Working directly with C-matrix in hdf5 file
// browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
//             else, indicates that R data is stored in hdf5 as column major
extern inline BigDataStatMeth::hdf5Dataset* hdf5_block_matrix_vector_sum_hdf5_transposed( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
        hsize_t hdf5_block, bool bparal, Rcpp::Nullable<int> threads  = R_NilValue)
{

    // Vector
    hsize_t K = dsA->nrows();
    hsize_t N = dsA->ncols();

    // Matrix
    hsize_t M = dsB->nrows();
    hsize_t L = dsB->ncols();
    
    try {

        if(hdf5_block == 1) {
            hdf5_block = 8192;
        }

        hsize_t isize = hdf5_block + 1;

        std::vector<hsize_t> stride = {1, 1};
        std::vector<hsize_t> block = {1, 1};

        dsC->createDataset( L, M, "real");
        
        std::vector<double> vdA( K * N );
        dsA->readDatasetBlock( {0, 0}, {K, N}, stride, block, vdA.data() );
        Eigen::Map<Eigen::VectorXd> A (vdA.data(), K*N );
        
        if(  K == M )
        { // Sum vector to every col
            
            for (hsize_t ii = 0; ii < L; ii += hdf5_block)
            {
                if( ii + hdf5_block > L )
                    isize = L - ii;
                
                hsize_t sizetoRead = getOptimBlockSize( L, hdf5_block, ii, isize);

                std::vector<double> vdB( K * sizetoRead );
                dsB->readDatasetBlock( {0, ii}, {K, sizetoRead}, stride, block, vdB.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), K, sizetoRead );

                Eigen::MatrixXd C = B.colwise() + A;
                
                std::vector<hsize_t> offset = { 0, ii};
                std::vector<hsize_t> count = { (hsize_t)C.rows(), (hsize_t)C.cols()};

                dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, true);

                if( ii + hdf5_block > L ) isize = hdf5_block + 1;
                if( sizetoRead > hdf5_block ) {
                    ii = ii - hdf5_block + sizetoRead; }
                
            }

        } else if(  K == L ) { // Sum vector to every row
            
            for (hsize_t ii = 0; ii < M; ii += hdf5_block)
            {
                if( ii + hdf5_block > M )
                    isize = M - ii;
                
                hsize_t sizetoRead = getOptimBlockSize( M, hdf5_block, ii, isize);
                
                std::vector<double> vdB( L * sizetoRead );
                dsB->readDatasetBlock( {ii, 0}, {sizetoRead, K}, stride, block, vdB.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), sizetoRead, K );    
                
                Eigen::MatrixXd C = B.rowwise() + A.transpose();

                std::vector<hsize_t> offset = {ii, 0};
                std::vector<hsize_t> count = {  (hsize_t)C.rows(), (hsize_t)C.cols() };

                // NOTE: Transposed to overcome discrepancies between rowMajor - colMajor between Eigen and HDF5
                dsC->writeDatasetBlock(Rcpp::wrap(C.transpose()), offset, count, stride, block, false); 

                if( ii + hdf5_block > M ) isize = hdf5_block + 1;
                if( sizetoRead > hdf5_block ) {
                    ii = ii - hdf5_block + sizetoRead; }
            }

        } else {
            Rcpp::Rcout<< "vector sum error: non-conformable arguments\n";
        }

    } catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception multiplication: "<<ex.what()<< " \n";
        return(dsC);
    }

    return(dsC);

}


}

#endif // BIGDATASTATMETH_ALGEBRA_SUM_HPP