#ifndef BIGDATASTATMETH_ALGEBRA_MSUM_HPP
#define BIGDATASTATMETH_ALGEBRA_MSUM_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
// #include <thread>

namespace BigDataStatMeth {



    template< typename T>  extern inline Rcpp::RObject Rcpp_matrix_sum ( T  A, T  B);
    template< typename T, typename U>  extern inline Rcpp::RObject Rcpp_matrix_vect_sum ( T  A, U  B);
    template< typename T>  extern inline Rcpp::RObject Rcpp_vector_sum ( T  A, T  B);
    
    template< typename T>  extern inline Rcpp::RObject Rcpp_matrix_blockSum ( T  A, T  B, Rcpp::Nullable<int> threads = R_NilValue);
    template< typename T>  extern inline Rcpp::RObject Rcpp_matrix_vector_blockSum( T  A, T  B, Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> threads);

    template< typename T>
    extern inline Eigen::MatrixXd Rcpp_block_matrix_vector_sum( T  A, T  B, hsize_t block_size, 
                                                  bool bparal, Rcpp::Nullable<int> threads = R_NilValue);



    extern inline void getBlockPositionsSizes( hsize_t maxPosition, hsize_t blockSize, std::vector<hsize_t>& starts, std::vector<hsize_t>& sizes ){
        
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


    template< typename T>
    extern inline Rcpp::RObject Rcpp_matrix_sum ( T  A, T  B)
    {
        
        // static_assert(std::is_same<T, Eigen::MatrixXd >::value ||
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value ||
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value ||
        //               std::is_same<T, Rcpp::NumericMatrix >::value,
        //               "Error - type not allowed");
        
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericMatrix m2 = Rcpp::as<Rcpp::NumericMatrix>(B);
        
        if( m.rows() == m2.rows() && m.cols() == m2.cols()) {
            Rcpp::NumericVector C = m + m2;
            C.attr("dim") = Rcpp::Dimension( m.rows(), m.cols());
            
            return(C);
            
        } else {
            Rcpp::Rcout<<"Error: non-conformable arguments";
        }
        
        return(R_NilValue);
        
    }
    
    
    // Suma per files o columnes depenent si la mida del vector es igual al nombre
    // de files o igual al nombre de columnes
    template< typename T, typename U>
    extern inline Rcpp::RObject Rcpp_matrix_vect_sum ( T  A, U  B)
    {
        
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(B);
        
        if( v.length() == m.rows()) {
            
            Rcpp::NumericMatrix C = Rcpp::no_init( m.rows(), m.cols());
            
            for( int i=0; i<m.cols(); i++) {
                C( Rcpp::_, i) = m( Rcpp::_, i) + v;  
            }    
            return(C);
            
        } else if( v.length() == m.cols()) {
            
            Rcpp::NumericMatrix C = Rcpp::no_init( m.rows(), m.cols());
            
            for( int i=0; i<m.rows(); i++) {
                C( i, Rcpp::_) = m( i, Rcpp::_) + v;  
            }    
            return(C);
            
        } else {
            Rcpp::Rcout<<"Error: non-conformable arguments";
        }
        
        return(R_NilValue);
    }


    template< typename T>
    extern inline Rcpp::RObject Rcpp_vector_sum ( T  A, T  B)
    {
        
        Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(A);
        Rcpp::NumericVector v2 = Rcpp::as<Rcpp::NumericVector>(B);
        
        if(v.size() == v2.size()) {
            Rcpp::NumericVector C = Rcpp::no_init( v.size());
            
            std::transform (v.begin(), v.end(), v2.begin(), C.begin(), std::plus<double>());
            
            C.attr("dim") = Rcpp::Dimension( C.size(), 1); 
            
            return(C);
        }
        
        return(R_NilValue);
        
    }
    
    
    template< typename T>
    extern inline Rcpp::RObject Rcpp_matrix_blockSum ( T  A, T  B, Rcpp::Nullable<int> threads)
    {
        
        // static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value ,
        //               "Error - type not allowed");
        
        Rcpp::NumericMatrix X = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericMatrix Y = Rcpp::as<Rcpp::NumericMatrix>(B);
        unsigned int ithreads;
        
        hsize_t N = X.rows();
        hsize_t M = X.cols();
        
        Rcpp::NumericMatrix C = Rcpp::no_init( N, M);
        
        hsize_t block_size; 
        
        try {
            
            std::vector<hsize_t> vsizetoRead;
            std::vector<hsize_t> vstart;
            
            std::vector<hsize_t> blockSize = getMatrixBlockSize( N, M);
            if(N < M) {
                block_size = blockSize.at(0);    
            } else {
                block_size = blockSize.at(1);
            }
            
            if(block_size > 0 ) {
                
                if( N == Y.rows() && M == Y.cols())
                {
                    hsize_t size = block_size + 1;
                    
                    ithreads = get_number_threads(threads, R_NilValue);
                    
                    getBlockPositionsSizes( N*M, block_size, vstart, vsizetoRead );
                    int chunks = vstart.size()/ithreads;
                    
                    #pragma omp parallel num_threads(ithreads) shared(A, B, C)
                    {
                    #pragma omp for schedule (dynamic, chunks)
                        for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                        {
                            
                            if( vstart[ii] + vsizetoRead[ii] >= N*M ) {
                                std::transform (X.begin() + vstart[ii], X.end(),
                                                Y.begin() + vstart[ii], C.begin() + vstart[ii], std::plus<double>());
                            } else {
                                std::transform (X.begin() + vstart[ii], X.begin() + vstart[ii] + vsizetoRead[ii],
                                                Y.begin() + vstart[ii], C.begin() + vstart[ii], std::plus<double>());   
                            }
                        }
                    }

                } else {
                    Rcpp::Rcout<<"matrix sum error: non-conformable arguments\n";
                    return(R_NilValue);
                }
                
            } else{
                Rcpp::Rcout<<"matrix sum error: Error whent computing block sizes\n";
                return(R_NilValue);
            }
            
        } catch(std::exception &ex) {
            Rcpp::Rcout<< ex.what();
            return(R_NilValue);
        }
        
        C.attr("dim") = Rcpp::Dimension( N, M);
        return(C);
        
    }
    
    
    
    
    template< typename T>
    extern inline Rcpp::RObject Rcpp_matrix_vector_blockSum( T  A, T  B,  
                                 Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> threads)
    {
        
        // NOTA: Per defecte, suma per columnes tal i com raja.... 
        
        // static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
        //               "Error - type not allowed");
        
        bool btransposed = false;
        unsigned int ithreads;
        hsize_t block_size;
        
        Rcpp::NumericMatrix X = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericVector Y = Rcpp::as<Rcpp::NumericVector>(B);
        Rcpp::NumericMatrix C;
        
        // Matrix
        hsize_t M = X.rows(),
                N = X.cols();
        
        // Vector
        hsize_t K = Y.length();
        
        try {
            
            if( K==N || K==M) {
                if ( K == N){
                    // Sum vector to every col
                    btransposed = true;
                    
                    X = Rcpp::transpose(X);
                    
                    hsize_t N = X.rows();
                    hsize_t M = X.cols();
                } 
                
                std::vector<hsize_t> vsizetoRead;
                std::vector<hsize_t> vstart;
                
                ithreads = get_number_threads(threads, bparal);
                
                C = Rcpp::no_init( M, N);
                
                block_size = getMatrixBlockSize( N, M).at(0);
                
                // minimum block size: 2 columns
                if(block_size <= 0 ) {
                    block_size = M*2;
                }
                
                // Mínimum block size: 2 columns
                getBlockPositionsSizes( M*N, block_size, vstart, vsizetoRead );
                int chunks = vstart.size()/ithreads;
                
                #pragma omp parallel num_threads(ithreads) shared(A, B, C, chunks)
                {
                    #pragma omp for schedule (dynamic, chunks) collapse(2)
                    for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                    {
                        // Duplicate vector
                        std::size_t const no_of_duplicates = vsizetoRead[ii] / Y.length();
                        
                        std::vector<double> v = Rcpp::as<std::vector<double> >(Y); 
                        v.reserve(Y.size() * no_of_duplicates);
                        auto end = std::end(v);
                        
                        for(std::size_t i = 1; i < no_of_duplicates; ++i)
                            v.insert(std::end(v), std::begin(v), end);
                        
                        // Sum vector to matrix by columns / rows
                        if( vstart[ii] + vsizetoRead[ii] >= M*N ) {
                            std::transform (X.begin() + vstart[ii], X.end(),
                                            v.begin(), C.begin() + vstart[ii], std::plus<double>());
                        } else {
                            std::transform (X.begin() + vstart[ii], X.begin() + vstart[ii] + vsizetoRead[ii],
                                            v.begin() , C.begin() + vstart[ii], std::plus<double>());   
                        }
                    }
                }
            
            } else {
                
                Rcpp::Rcout<< "vector sum error: non-conformable arguments\n";
                return(R_NilValue);
            }
            
    
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception multiplication: "<<ex.what()<< " \n";
            return(R_NilValue);
        }
        
        if(btransposed == true){
            Rcpp::transpose(C);
        } 
        
        C.attr("dim") = Rcpp::Dimension( M, N);
        return(C);
    
    }





// 
// 
// // IMPORTANT : R data stored as Row-major in hdf5 file (stored transposed data from R)
// // 
// // This function differs from "hdf5_block_matrix_sum_hdf5_transposed" in parameters, 
// // 
// // Working directly with C-matrix in hdf5 file
// // browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
// //             else, indicates that R data is stored in hdf5 as column major
// extern inline BigDataStatMeth::hdf5Dataset*  hdf5_block_matrix_sum_hdf5_indatasets_transposed( 
//         BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
//         hsize_t hdf5_block, hsize_t mem_block_size, bool bparal, bool browmajor,  Rcpp::Nullable<int> threads  = R_NilValue)
// {
//     
//     try {
//         
//         hsize_t K = dsA->nrows();
//         hsize_t N = dsA->ncols();
//         
//         if( K == dsB->nrows() && N == dsB->ncols())
//         {
//             
//             hsize_t isize = hdf5_block + 1;
//             
//             std::vector<hsize_t> stride = {1, 1};
//             std::vector<hsize_t> block = {1, 1};
//             
//             dsC->createDataset( N, K, "real"); 
//             
//             if( K<N ) {
//                 int volta = 0;
//                 for (hsize_t ii = 0; ii < N; ii += hdf5_block)
//                 {
//                     
//                     if( ii + hdf5_block > N ) {
//                         isize = N - ii; }
//                     
//                     hsize_t sizetoRead = getOptimBlockSize( N, hdf5_block, ii, isize);
//                     Rcpp::Rcout<<"\nEstem llegint la mida de block:"<< sizetoRead<<" volta: "<<volta<<"\n";
//                     volta = volta + 1;
//                     
//                     std::vector<double> vdA( K * sizetoRead ); 
//                     dsA->readDatasetBlock( {0, ii}, { K, sizetoRead}, stride, block, vdA.data() );
//                     Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), K, sizetoRead );
//                     
//                     std::vector<double> vdB( K * sizetoRead ); 
//                     dsB->readDatasetBlock( {0, ii}, {K, sizetoRead}, stride, block, vdB.data() );
//                     Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), K, sizetoRead );
//                         
//                     Eigen::MatrixXd C = A + B;
//                     
//                     std::vector<hsize_t> offset = { 0, ii };
//                     std::vector<hsize_t> count = { sizetoRead, K };
//                     
//                     dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, true);
//                     
//                     if( ii + hdf5_block > N ) isize = hdf5_block + 1;
//                     if( sizetoRead > hdf5_block ) {
//                         ii = ii - hdf5_block + sizetoRead; }
//                 }
//                 
//             } else {
//                 
//                 int volta = 0;
//                 
//                 for (hsize_t ii = 0; ii < K; ii += hdf5_block)
//                 {
//                     
//                     if( ii + hdf5_block > K ) 
//                         isize = K - ii;
//                     
//                     hsize_t sizetoRead = getOptimBlockSize( K, hdf5_block, ii, isize);
//                     Rcpp::Rcout<<"\nEstem llegint la mida de block:"<< sizetoRead<<" volta: "<<volta<<"\n";
//                     volta = volta + 1;
//                     
//                     std::vector<double> vdA( sizetoRead * N ); 
//                     dsA->readDatasetBlock( {ii, 0}, { sizetoRead, N}, stride, block, vdA.data() );
//                     Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), sizetoRead, N );
//                     
//                     std::vector<double> vdB( sizetoRead * N); 
//                     dsB->readDatasetBlock( {ii, 0}, {sizetoRead, N}, stride, block, vdB.data() );
//                     Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), sizetoRead, N );
//                     
//                     Eigen::MatrixXd C = A + B;
//                     
//                     std::vector<hsize_t> offset = { ii, 0 };
//                     std::vector<hsize_t> count = { N, sizetoRead };
// 
//                     dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, true);
//                     
//                     if( ii + hdf5_block > K ) isize = hdf5_block + 1;
//                     if( sizetoRead > hdf5_block ) {
//                         ii = ii - hdf5_block + sizetoRead; }
//                 }
//             }
//             
//         } else {
//             Rcpp::Rcout<<"matrix sum error: non-conformable arguments\n";
//         }
//         
//     } catch(std::exception &ex) {
//         Rcpp::Rcout<< ex.what();
//         return(dsC);
//     }
//     
//     return(dsC);
// }
// 
// 
// 
// // This function makes summatori with a matrix and a vector, returns a matrix.
// //
// // Working directly with C-matrix in hdf5 file
// // browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
// //             else, indicates that R data is stored in hdf5 as column major
// extern inline BigDataStatMeth::hdf5Dataset* hdf5_block_matrix_vector_sum_hdf5_transposed( 
//         BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
//         hsize_t hdf5_block, bool bparal, bool browmajor, Rcpp::Nullable<int> threads  = R_NilValue)
// {
// 
//     // Vector
//     hsize_t K = dsA->nrows();
//     hsize_t N = dsA->ncols();
// 
//     // Matrix
//     hsize_t M = dsB->nrows();
//     hsize_t L = dsB->ncols();
//     
//     try {
// 
//         if(hdf5_block == 1) {
//             hdf5_block = 8192;
//         }
// 
//         hsize_t isize = hdf5_block + 1;
// 
//         std::vector<hsize_t> stride = {1, 1};
//         std::vector<hsize_t> block = {1, 1};
// 
//         dsC->createDataset( L, M, "real");
//         
//         std::vector<double> vdA( K * N );
//         dsA->readDatasetBlock( {0, 0}, {K, N}, stride, block, vdA.data() );
//         Eigen::Map<Eigen::VectorXd> A (vdA.data(), K*N );
//         
//         if(  K == M )
//         { // Sum vector to every col
//             
//             for (hsize_t ii = 0; ii < L; ii += hdf5_block)
//             {
//                 if( ii + hdf5_block > L )
//                     isize = L - ii;
//                 
//                 hsize_t sizetoRead = getOptimBlockSize( L, hdf5_block, ii, isize);
// 
//                 std::vector<double> vdB( K * sizetoRead );
//                 dsB->readDatasetBlock( {0, ii}, {K, sizetoRead}, stride, block, vdB.data() );
//                 Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), K, sizetoRead );
// 
//                 Eigen::MatrixXd C = B.colwise() + A;
//                 
//                 std::vector<hsize_t> offset = { 0, ii};
//                 std::vector<hsize_t> count = { (hsize_t)C.rows(), (hsize_t)C.cols()};
// 
//                 dsC->writeDatasetBlock(Rcpp::wrap(C), offset, count, stride, block, true);
// 
//                 if( ii + hdf5_block > L ) isize = hdf5_block + 1;
//                 if( sizetoRead > hdf5_block ) {
//                     ii = ii - hdf5_block + sizetoRead; }
//                 
//             }
// 
//         } else if(  K == L ) { // Sum vector to every row
//             
//             for (hsize_t ii = 0; ii < M; ii += hdf5_block)
//             {
//                 if( ii + hdf5_block > M )
//                     isize = M - ii;
//                 
//                 hsize_t sizetoRead = getOptimBlockSize( M, hdf5_block, ii, isize);
//                 
//                 std::vector<double> vdB( L * sizetoRead );
//                 dsB->readDatasetBlock( {ii, 0}, {sizetoRead, K}, stride, block, vdB.data() );
//                 Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> B (vdB.data(), sizetoRead, K );    
//                 
//                 Eigen::MatrixXd C = B.rowwise() + A.transpose();
// 
//                 std::vector<hsize_t> offset = {ii, 0};
//                 std::vector<hsize_t> count = {  (hsize_t)C.rows(), (hsize_t)C.cols() };
// 
//                 // NOTE: Transposed to overcome discrepancies between rowMajor - colMajor between Eigen and HDF5
//                 dsC->writeDatasetBlock(Rcpp::wrap(C.transpose()), offset, count, stride, block, false); 
// 
//                 if( ii + hdf5_block > M ) isize = hdf5_block + 1;
//                 if( sizetoRead > hdf5_block ) {
//                     ii = ii - hdf5_block + sizetoRead; }
//             }
// 
//         } else {
//             Rcpp::Rcout<< "vector sum error: non-conformable arguments\n";
//         }
// 
//     } catch(std::exception& ex) {
//         Rcpp::Rcout<< "c++ exception multiplication: "<<ex.what()<< " \n";
//         return(dsC);
//     }
// 
//     return(dsC);
// 
// }


}

#endif // BIGDATASTATMETH_ALGEBRA_MSUM_HPP