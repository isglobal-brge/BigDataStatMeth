#ifndef BIGDATASTATMETH_ALGEBRA_CROSSPROD_HPP
#define BIGDATASTATMETH_ALGEBRA_CROSSPROD_HPP

#include "hdf5Algebra/multiplication.hpp"

namespace BigDataStatMeth {


extern inline BigDataStatMeth::hdf5Dataset* crossprod( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, 
        BigDataStatMeth::hdf5Dataset* dsC, hsize_t hdf5_block, hsize_t mem_block_size, 
        bool bparal, bool browmajor, Rcpp::Nullable<int> threads  = R_NilValue) 

{
    
    try {

        hsize_t N = dsA->nrows();
        hsize_t K = dsA->ncols();
        hsize_t M = dsB->nrows();
        hsize_t L = dsB->ncols();
        
        if( K == L)
        {
            hsize_t isize = hdf5_block + 1,
                    ksize = hdf5_block + 1,
                    jsize = hdf5_block + 1;

            std::vector<hsize_t> stride = {1, 1};
            std::vector<hsize_t> block = {1, 1};

            dsC->createDataset( N, M, "real");

            for (hsize_t ii = 0; ii < N; ii += hdf5_block)
            {
                hsize_t iRowsA = getOptimBlockSize( N, hdf5_block, ii, isize);
                
                if( ii + hdf5_block > N ) isize = N - ii;
                
                for (hsize_t jj = 0; jj < M; jj += hdf5_block)
                {
                    hsize_t iRowsB = getOptimBlockSize( M, hdf5_block, jj, jsize);
                    
                    if( jj + hdf5_block > M) jsize = M - jj;
                    
                    for(hsize_t kk = 0; kk < K; kk += hdf5_block)
                    {
                        if( kk + hdf5_block > K ) ksize = K - kk;

                        hsize_t iColsA = getOptimBlockSize( K, hdf5_block, kk, ksize),
                            iColsB = iColsA;
                        
                        Eigen::MatrixXd C;
                        Eigen::MatrixXd B;
                        
                        std::vector<double> vdA( iRowsA * iColsA);
                        dsA->readDatasetBlock( {ii, kk}, {iRowsA,iColsA}, stride, block, vdA.data() );
                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), iRowsA, iColsA );

                        {
                            std::vector<double> vdB( iRowsB * iColsB);
                            dsB->readDatasetBlock( {jj, kk}, {iRowsB, iColsB}, stride, block, vdB.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_B (vdB.data(), iRowsB, iColsB);   
                            B = tmp_B.transpose();
                        }
                        
                        {
                            std::vector<double> vdC( iRowsB * iRowsA);
                            dsC->readDatasetBlock( {jj, ii}, {iRowsB, iRowsA}, stride, block, vdC.data() );
                            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmp_C (vdC.data(), iRowsB, iRowsA);
                            C = tmp_C.transpose();
                        }
                           
                        C = C + (A * B);
                        
                        std::vector<hsize_t> offset = {jj,ii};
                        std::vector<hsize_t> count = {iRowsB, iRowsA};
                        
                        dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, false);

                        // Readjust counters
                        if( kk + hdf5_block > K ) ksize = hdf5_block + 1;
                        if( iColsA > hdf5_block ) {
                            kk = kk - hdf5_block + iColsA; }
                    }

                    if( jj + hdf5_block > M ) jsize = hdf5_block + 1;
                    if( iRowsB > hdf5_block ) {
                        jj = jj - hdf5_block + iRowsB; }
                }

                if( ii + hdf5_block > N ) isize = hdf5_block + 1;
                if( iRowsA > hdf5_block ) {
                    ii = ii - hdf5_block + iRowsA; }
            }

        } else {
            throw std::range_error("non-conformable arguments");
        }

    } catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception crossprod: "<<ex.what()<< " \n";
        return(dsC);
    }

    return(dsC);
}

}

#endif // BIGDATASTATMETH_ALGEBRA_CROSSPROD_HPP