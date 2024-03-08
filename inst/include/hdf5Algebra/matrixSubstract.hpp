#ifndef BIGDATASTATMETH_ALGEBRA_SUBSTRACT_HPP
#define BIGDATASTATMETH_ALGEBRA_SUBSTRACT_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
// #include <thread>

namespace BigDataStatMeth {


extern inline BigDataStatMeth::hdf5Dataset*  hdf5_block_matrix_substract_hdf5_indatasets_transposed( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
        hsize_t hdf5_block, hsize_t mem_block_size, bool bparal, bool browmajor,  Rcpp::Nullable<int> threads);

extern inline BigDataStatMeth::hdf5Dataset* hdf5_block_matrix_vector_substract_hdf5_transposed( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
        hsize_t hdf5_block, bool bparal, bool browmajor, Rcpp::Nullable<int> threads);

// IMPORTANT : R data stored as Row-major in hdf5 file (stored transposed data from R)
// Working directly with C-matrix in hdf5 file
// browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
//             else, indicates that R data is stored in hdf5 as column major
extern inline BigDataStatMeth::hdf5Dataset*  hdf5_block_matrix_substract_hdf5_indatasets_transposed( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
        hsize_t hdf5_block, hsize_t mem_block_size, bool bparal, bool browmajor,  Rcpp::Nullable<int> threads  = R_NilValue)
{
    
    try {
        
        hsize_t K = dsA->nrows();
        hsize_t N = dsA->ncols();
        
        if( K == dsB->nrows() && N == dsB->ncols())
        {
            
            hsize_t isize = hdf5_block + 1;
            
            std::vector<hsize_t> stride = {1, 1};
            std::vector<hsize_t> block = {1, 1};
            
            dsC->createDataset( N, K, "real"); 
            
            if( K<N ) {
                for (hsize_t ii = 0; ii < N; ii += hdf5_block)
                {
                    // Rcpp::Rcout<<"\n K<N \n";
                    
                    if( ii + hdf5_block > N ) {
                        isize = N - ii; }
                    
                    hsize_t sizetoRead = getOptimBlockSize( N, hdf5_block, ii, isize);
                    
                    std::vector<double> vdA( K * sizetoRead ); 
                    dsA->readDatasetBlock( {0, ii}, { K, sizetoRead}, stride, block, vdA.data() );
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), K, sizetoRead );
                    
                    std::vector<double> vdB( K * sizetoRead ); 
                    dsB->readDatasetBlock( {0, ii}, {K, sizetoRead}, stride, block, vdB.data() );
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), K, sizetoRead );
                        
                    Eigen::MatrixXd C = A - B;
                    
                    std::vector<hsize_t> offset = { 0, ii };
                    std::vector<hsize_t> count = { sizetoRead, K };
                    
                    dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, true);
                    
                    if( ii + hdf5_block > N ) isize = hdf5_block + 1;
                    if( sizetoRead > hdf5_block ) {
                        ii = ii - hdf5_block + sizetoRead; }
                }
                
            } else {
                for (hsize_t ii = 0; ii < K; ii += hdf5_block)
                {
                    
                    // Rcpp::Rcout<<"\n K>=N \n";
                    if( ii + hdf5_block > K ) 
                        isize = K - ii;
                    
                    hsize_t sizetoRead = getOptimBlockSize( K, hdf5_block, ii, isize);
                    
                    std::vector<double> vdA( sizetoRead * N ); 
                    dsA->readDatasetBlock( {ii, 0}, { sizetoRead, N}, stride, block, vdA.data() );
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), sizetoRead, N );
                    
                    std::vector<double> vdB( sizetoRead * N); 
                    dsB->readDatasetBlock( {ii, 0}, {sizetoRead, N}, stride, block, vdB.data() );
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), sizetoRead, N );
                    
                    Eigen::MatrixXd C = A - B;
                    
                    std::vector<hsize_t> offset = { ii, 0 };
                    std::vector<hsize_t> count = { N, sizetoRead };

                    dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, true);
                    
                    if( ii + hdf5_block > K ) isize = hdf5_block + 1;
                    if( sizetoRead > hdf5_block ) {
                        ii = ii - hdf5_block + sizetoRead; }
                }
            }
            
        } else {
            Rcpp::Rcout<<"matrix substract error: non-conformable arguments\n";
        }
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return(dsC);
    }
    
    return(dsC);
}



// This function makes substract using a matrix and a vector, returns a matrix.
//
// Working directly with C-matrix in hdf5 file
// browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
//             else, indicates that R data is stored in hdf5 as column major
extern inline BigDataStatMeth::hdf5Dataset* hdf5_block_matrix_vector_substract_hdf5_transposed( 
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
        hsize_t hdf5_block, bool bparal, bool browmajor, Rcpp::Nullable<int> threads  = R_NilValue)
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

                Eigen::MatrixXd C = B.colwise() - A;
                
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
                
                Eigen::MatrixXd C = B.rowwise() - A.transpose();

                std::vector<hsize_t> offset = {ii, 0};
                std::vector<hsize_t> count = {  (hsize_t)C.rows(), (hsize_t)C.cols() };

                // NOTE: Transposed to overcome discrepancies between rowMajor - colMajor between Eigen and HDF5
                dsC->writeDatasetBlock(Rcpp::wrap(C.transpose()), offset, count, stride, block, false); 

                if( ii + hdf5_block > M ) isize = hdf5_block + 1;
                if( sizetoRead > hdf5_block ) {
                    ii = ii - hdf5_block + sizetoRead; }
            }

        } else {
            Rcpp::Rcout<< "vector substract error: non-conformable arguments\n";
        }

    } catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception substraction: "<<ex.what()<< " \n";
        return(dsC);
    }

    return(dsC);

}


}

#endif // BIGDATASTATMETH_ALGEBRA_SUBSTRACT_HPP