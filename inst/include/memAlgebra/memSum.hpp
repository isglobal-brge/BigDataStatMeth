#ifndef BIGDATASTATMETH_ALGEBRA_MSUM_HPP
#define BIGDATASTATMETH_ALGEBRA_MSUM_HPP

#include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include <thread>

namespace BigDataStatMeth {



    template< typename T>  extern inline Rcpp::RObject Rcpp_matrix_sum ( T  A, T  B);
    template< typename T, typename U>  extern inline Rcpp::RObject Rcpp_matrix_vect_sum ( T  A, U  B);
    template< typename T>  extern inline Rcpp::RObject Rcpp_vector_sum ( T  A, T  B);
    
    template< typename T>  extern inline Rcpp::RObject Rcpp_matrix_blockSum ( T  A, T  B, Rcpp::Nullable<int> threads = R_NilValue);



// 
// template< typename T> extern inline Eigen::MatrixXd Rcpp_block_matrix_sum ( T  A, T  B, 
//                                 hsize_t block_size, bool bparal, Rcpp::Nullable<int> threads = R_NilValue);


template< typename T>
extern inline Eigen::MatrixXd Rcpp_block_matrix_vector_sum( T  A, T  B, hsize_t block_size, 
                                              bool bparal, Rcpp::Nullable<int> threads = R_NilValue);

// 
// template< typename T>
// extern inline Eigen::MatrixXd Rcpp_block_matrix_sum ( T  A, T  B, hsize_t block_size, 
//                             bool bparal, Rcpp::Nullable<int> threads)
// {
//     
//     static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
//                   std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
//                   std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
//                   "Error - type not allowed");
//     
//     Eigen::Map<Eigen::MatrixXd > X = A;
//     Eigen::Map<Eigen::MatrixXd > Y = B;
//     
//     hsize_t K = A.rows();
//     hsize_t N = A.cols();
//     
//     Eigen::MatrixXd C = Eigen::MatrixXd::Zero(K,N);
//     
//     try {
//         
//         
//         
//         if( N == Y.cols() && K == Y.rows())
//         {
//             return(A+B);
//              
//             //  hsize_t isize = block_size + 1;
//             // 
//             // if( K<N ) {
//             //     #pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C) 
//             //     {
//             //         #pragma omp for schedule (dynamic) 
//             //         for (hsize_t ii = 0; ii < N; ii += block_size)
//             //         {
//             //             
//             //             if( ii + block_size > N ) {
//             //                 isize = N - ii; }
//             //             
//             //             hsize_t sizetoRead = getOptimBlockSize( N, block_size, ii, isize);
//             //             
//             //             C.block(0, ii, K, sizetoRead) = A.block(0, ii, K, sizetoRead) + B.block(0, ii, K, sizetoRead);
//             //             
//             //             if( ii + block_size > N ) isize = block_size + 1;
//             //             if( sizetoRead > block_size ) {
//             //                 ii = ii - block_size + sizetoRead; }
//             //         }
//             //     }
//             //     
//             // } else {
//             //     
//             //     #pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C) 
//             //     {
//             //         #pragma omp for schedule (dynamic) 
//             //         
//             //         for (hsize_t ii = 0; ii < K; ii += block_size)
//             //         {
//             //             
//             //             if( ii + block_size > K ) 
//             //                 isize = K - ii;
//             //             
//             //             hsize_t sizetoRead = getOptimBlockSize( K, block_size, ii, isize);
//             //             
//             //             C.block(ii, 0, sizetoRead, N) = A.block(ii, 0, sizetoRead, N) + B.block(ii, 0, sizetoRead, N);
//             //             
//             //             if( ii + block_size > K ) isize = block_size + 1;
//             //             if( sizetoRead > block_size ) {
//             //                 ii = ii - block_size + sizetoRead; }
//             //         }
//             //     }
//             // }
//             
//         } else {
//             Rcpp::Rcout<<"matrix sum error: non-conformable arguments\n";
//         }
//         
//     } catch(std::exception &ex) {
//         Rcpp::Rcout<< ex.what();
//         return(C);
//     }
//     
//     return(C);
// 
// }



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
    
    
    
    template< typename T, typename U>
    extern inline Rcpp::RObject Rcpp_matrix_vect_sum ( T  A, U  B)
    {
        
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(B);
        
        if( v.length() == m.cols()) {
            Rcpp::NumericMatrix C = Rcpp::no_init( m.rows(), m.cols());
            
            for( int i=0; i<m.cols(); i++) {
                C( Rcpp::_, i) = m( Rcpp::_, i) + v;  
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
            
            std::transform (v.begin(), v.end(), v2.begin(), v.begin(), std::plus<double>());
            
            // Rcpp::NumericVector vC = v + v2;
            
            // vC.attr("dim") = Rcpp::Dimension( v.size(), 1); 
            
            return( v);
        }
        
        return(R_NilValue);
        
    }


    
    
    // 
    // template< typename T>
    // extern inline Rcpp::RObject Rcpp_matrix_blockSum ( T  A, T  B, Rcpp::Nullable<int> threads)
    // {
    //     
    //     // static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
    //     //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
    //     //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value ,
    //     //               "Error - type not allowed");
    //     
    //     Rcpp::NumericMatrix X = Rcpp::as<Rcpp::NumericMatrix>(A);
    //     Rcpp::NumericMatrix Y = Rcpp::as<Rcpp::NumericMatrix>(B);
    //     
    //     hsize_t N = X.rows();
    //     hsize_t M = X.cols();
    //     
    //     std::vector<hsize_t> block_size; 
    //     
    //     Eigen::MatrixXd C;
    //     
    //     try {
    //         
    //         block_size = getMatrixBlockSize( N, M); 
    //         
    //         if(block_size[0] > 0 && block_size[1] > 0) {
    //             
    //             C = Eigen::MatrixXd::Zero( N, M);
    //             
    //             if( N == Y.rows() && M == Y.cols())
    //             {
    //                 hsize_t xsize = block_size[0] + 1;
    //                 hsize_t ysize = block_size[1] + 1;
    //                 
    //                 #pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C)
    //                 {
    //                 #pragma omp for schedule (dynamic)
    //                     for (hsize_t ii = 0; ii < N; ii += block_size[0])
    //                     {
    //                         int tid = omp_get_thread_num();
    //                         printf("Hello world from omp thread %d\n", tid);
    //                         
    //                         if( ii + block_size[0] > N ) {
    //                             xsize = N - ii; }
    //                         
    //                         hsize_t sizetoRead_x = getOptimBlockSize( N, block_size[0], ii, xsize);
    //                         
    //                         for (hsize_t jj = 0; jj < M; jj += block_size[1]) {
    //                             
    //                             if( jj + block_size[1] > M ) {
    //                                 ysize = M - jj; }
    //                             
    //                             hsize_t sizetoRead_y = getOptimBlockSize( M, block_size[1], jj, ysize);
    //                             
    //                             Rcpp::Rcout<<"Llegim posicions: \n\tFiles: "<<ii<<" x "<<ii+sizetoRead_x-1<<"\n\tColumnes: "<<jj<<" x "<<jj+sizetoRead_y-1<<"\n\tMida Block: "<<sizetoRead_x<<" x "<<sizetoRead_y<<"\n";
    //                             Rcpp::Rcout<<"Escrivim a la posició començant: \n\tCoordenades: "<<ii<<" x "<<jj<<"\n\tMida Block: "<<sizetoRead_x<<" x "<<sizetoRead_y<<"\n";
    //                             Rcpp::Rcout<<"\n Una mica de parafernalia per saltar una mica més per si de cas....\n";
    //                             
    //                             //// CREC QUE AIXÒ ES PODRIA FER AMB VECTORS D'UNA FORMA MOLT MES FÀCIL I RÀPIDA... PARTINT ELS VECTORS DIRECTAMENT
    //                             //// EN BLOCKS I FENT EL SUMATORI 1 A 1.... 
    //                             //// FINALMENT-> CONVERTINT-HO DE NOU A NUMERICMATRIX (TAL I COM FAIG AL FER EL SUMATORI DE VECTORS!!)
    //                             
    //                             C.block(ii,jj, sizetoRead_x, sizetoRead_y) = Rcpp::as<Eigen::MatrixXd>(
    //                                 Rcpp_matrix_sum( Rcpp::wrap(X( Rcpp::Range(ii, ii+sizetoRead_x-1 ), Rcpp::Range(jj, jj+sizetoRead_y-1))), 
    //                                                  Rcpp::wrap(Y( Rcpp::Range(ii, ii+sizetoRead_x-1 ), Rcpp::Range(jj, jj+sizetoRead_y-1)))));
    //                             
    //                             if( jj + block_size[1] > M ) ysize = block_size[1] + 1;
    //                             if( sizetoRead_y > block_size[1] ) {
    //                                 jj = jj - block_size[1] + sizetoRead_y; }
    //                         }
    //                         
    //                         if( ii + block_size[0] > N ) xsize = block_size[0] + 1;
    //                         if( sizetoRead_x > block_size[0] ) {
    //                             ii = ii - block_size[0] + sizetoRead_x; }
    //                     }
    //                 }
    //                 
    //             } else {
    //                 Rcpp::Rcout<<"matrix sum error: non-conformable arguments\n";
    //                 return(R_NilValue);
    //             }
    //             
    //         } else{
    //             Rcpp::Rcout<<"matrix sum error: Error whent computing block sizes\n";
    //             return(R_NilValue);
    //         }
    //         
    //     } catch(std::exception &ex) {
    //         Rcpp::Rcout<< ex.what();
    //         return(R_NilValue);
    //     }
    //     
    //     return(Rcpp::wrap(C));
    //     
    // }
    // 
    
    
    
    
    
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
            
            // block_size = getVectorBlockSize( N*M); 
            
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
                    
                    if(threads.isNotNull()) {
                        if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
                            ithreads = Rcpp::as<int> (threads);
                        } else {
                            ithreads = getDTthreads(0, true);
                        }
                    } else {
                        ithreads = getDTthreads(0, true);
                    }
                    
                    #pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(A, B, C)
                    {
                    #pragma omp for schedule (dynamic)
                        for (hsize_t ii = 0; ii < N*M; ii += block_size)
                        {
                            
                            if( ii + block_size > N*M ) {
                                size = N*M - ii; }
                            
                            hsize_t sizetoRead = getOptimBlockSize( N*M, block_size, ii, size);
                            
                            if( ii + sizetoRead >= N*M ) {
                                std::transform (X.begin() + ii, X.end(),
                                                Y.begin() + ii, C.begin() + ii, std::plus<double>());
                            } else {
                                std::transform (X.begin() + ii, X.begin() + ii + sizetoRead,
                                                Y.begin() + ii, C.begin() + ii, std::plus<double>());   
                            }
                            
                            if( ii + block_size > N*M ) size = block_size + 1;
                            if( sizetoRead > block_size ) {
                                ii = ii - block_size + sizetoRead; }
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
extern inline Eigen::MatrixXd Rcpp_block_matrix_vector_sum( T  A, T  B, hsize_t block_size, 
                            bool bparal, Rcpp::Nullable<int> threads)
{
    
    // static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
    //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
    //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
    //               "Error - type not allowed");
    
    Eigen::Map<Eigen::MatrixXd > X = A;
    Eigen::Map<Eigen::MatrixXd > Y = B;
    
    // Vector
    hsize_t K = A.rows();
    hsize_t N = A.cols();
    
    // Matrix
    hsize_t M = B.rows();
    hsize_t L = B.cols();
    
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M,L);
    
    try {

        if(block_size == 1) {
            block_size = 8192;
        }

        hsize_t isize = block_size + 1;

        /*
        if(  K == M )
        { // Sum vector to every col

            for (hsize_t ii = 0; ii < L; ii += block_size)
            {
                if( ii + block_size > L )
                    isize = L - ii;

                hsize_t sizetoRead = getOptimBlockSize( L, block_size, ii, isize);

                C.block(0, ii, K, sizetoRead) = B.block(0, ii, K, sizetoRead).colwise() + A;

                if( ii + block_size > L ) isize = block_size + 1;
                if( sizetoRead > block_size ) {
                    ii = ii - block_size + sizetoRead; }

            }

        } else if(  K == L ) { // Sum vector to every row

            for (hsize_t ii = 0; ii < M; ii += block_size)
            {
                if( ii + block_size > M )
                    isize = M - ii;

                hsize_t sizetoRead = getOptimBlockSize( M, block_size, ii, isize);

                C.block(ii, 0, sizetoRead, K) = B.block(ii, 0, sizetoRead, K).colwise() + A.transpose();

                if( ii + block_size > M ) isize = block_size + 1;
                if( sizetoRead > block_size ) {
                    ii = ii - block_size + sizetoRead; }
            }

        } else {
            Rcpp::Rcout<< "vector sum error: non-conformable arguments\n";
        }
         */

    } catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception multiplication: "<<ex.what()<< " \n";
        return(C);
    }

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