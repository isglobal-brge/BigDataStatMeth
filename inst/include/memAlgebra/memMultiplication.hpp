#ifndef BIGDATASTATMETH_ALGEBRA_MEM_MULTIPLICATION_HPP
#define BIGDATASTATMETH_ALGEBRA_MEM_MULTIPLICATION_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
// #include <thread>

namespace BigDataStatMeth {

    // extern inline Eigen::MatrixXd Rcpp_block_matrix_mul( Eigen::MatrixXd A, Eigen::MatrixXd B, Rcpp::Nullable<int>  iblock_size);
    template<typename T, typename U> extern inline Eigen::MatrixXd Rcpp_block_matrix_mul( T X, U Y, Rcpp::Nullable<int>  iblock_size);

    
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
    
    
    extern inline void getBlockPositionsSizes_mat( hsize_t maxPosition, hsize_t blockSize, std::vector<hsize_t>& starts, std::vector<hsize_t>& sizes ){
        
        hsize_t isize = blockSize + 1;
        
        for (hsize_t ii = 0; ii < maxPosition; ii += blockSize)
        {
            if( ii + blockSize > maxPosition ) {
                isize = maxPosition - ii; }
            
            hsize_t sizetoRead = getOptimBlockSize( maxPosition, blockSize, ii, isize);
            
            starts.push_back(ii);
            sizes.push_back(sizetoRead);
            
            if( sizetoRead > blockSize ) {
                ii = ii - blockSize + sizetoRead; }
        }
        
    }
    
    
    // // In-memory execution - Serial version by Blocks
    // extern inline Eigen::MatrixXd Rcpp_block_matrix_mul( Eigen::MatrixXd A, Eigen::MatrixXd B, Rcpp::Nullable<int>  iblock_size)
    // {
    //     
    //     Eigen::MatrixXd C;
    //     
    //     try{
    //         
    //         int M = A.rows(),
    //             K = A.cols(),
    //             N = B.cols(),
    //             block_size;
    //         
    //         if( iblock_size.isNotNull()) {
    //             block_size =  Rcpp::as<int>(iblock_size);
    //         } else {
    //             block_size =  MAXBLOCKSIZE/3;
    //         }
    //         
    //         if( K == B.rows())
    //         {
    //             C = Eigen::MatrixXd::Zero(M,N) ; 
    //             
    //             int isize = block_size+1,
    //                 ksize = block_size+1,
    //                 jsize = block_size+1;
    //             
    //             for (int ii = 0; ii < M; ii += block_size)
    //             {
    //                 if( ii + block_size > M ) isize = M - ii;
    //                 for (int jj = 0; jj < N; jj += block_size)
    //                 {
    //                     if( jj + block_size > N) jsize = N - jj;
    //                     for(int kk = 0; kk < K; kk += block_size)
    //                     {
    //                         if( kk + block_size > K ) ksize = K - kk;
    //                         
    //                         hsize_t minii = std::min(block_size,isize),
    //                                 minjj = std::min(block_size,jsize),
    //                                 minkk = std::min(block_size,ksize);
    //                         
    //                         C.block(ii, jj, minii, minjj) =  C.block(ii, jj, minii, minjj) + 
    //                             (A.block(ii, kk, minii, minkk) * B.block(kk, jj, minkk, minjj));
    //                         
    //                         if( kk + block_size > K ) ksize = block_size+1;
    //                     }
    //                     if( jj + block_size > N ) jsize = block_size+1;
    //                 }
    //                 if( ii + block_size > M ) isize = block_size+1;
    //             }
    //             
    //         } else {
    //             throw std::range_error("non-conformable arguments");
    //         }    
    //         
    //     } catch(std::exception &ex) {
    //         Rcpp::Rcout<<"c++ error : Bblock_matrix_mul : " <<ex.what();
    //         return(C);
    //         
    //     } catch(...) { 
    //         ::Rf_error("c++ exception in Bblock_matrix_mul (unknown reason)"); 
    //     }
    //     
    //     return(C);
    // }
    // 
    // 
   
   // In-memory execution - Serial version by Blocks
   template<typename T, typename U>
   extern inline Eigen::MatrixXd Rcpp_block_matrix_mul( T X, U Y, Rcpp::Nullable<int>  iblock_size)
   {
       
       Eigen::MatrixXd C;
       
       try{
           
           // int M = A.rows(),
           //     K = A.cols(),
           //     N = B.cols(),
           //     block_size;
           // 
           // if( iblock_size.isNotNull()) {
           //     block_size =  Rcpp::as<int>(iblock_size);
           // } else {
           //     block_size =  MAXBLOCKSIZE/3;
           // }
           // 
           // if( K == B.rows())
           // {
           //     C = Eigen::MatrixXd::Zero(M,N) ; 
           //     
           //     int isize = block_size+1,
           //         ksize = block_size+1,
           //         jsize = block_size+1;
           //     
           //     for (int ii = 0; ii < M; ii += block_size)
           //     {
           //         if( ii + block_size > M ) isize = M - ii;
           //         for (int jj = 0; jj < N; jj += block_size)
           //         {
           //             if( jj + block_size > N) jsize = N - jj;
           //             for(int kk = 0; kk < K; kk += block_size)
           //             {
           //                 if( kk + block_size > K ) ksize = K - kk;
           //                 
           //                 hsize_t minii = std::min(block_size,isize),
           //                     minjj = std::min(block_size,jsize),
           //                     minkk = std::min(block_size,ksize);
           //                 
           //                 C.block(ii, jj, minii, minjj) =  C.block(ii, jj, minii, minjj) + 
           //                     (A.block(ii, kk, minii, minkk) * B.block(kk, jj, minkk, minjj));
           //                 
           //                 if( kk + block_size > K ) ksize = block_size+1;
           //             }
           //             if( jj + block_size > N ) jsize = block_size+1;
           //         }
           //         if( ii + block_size > M ) isize = block_size+1;
           //     }
           //     
           // } else {
           //     throw std::range_error("non-conformable arguments");
           // }    
           
           
           static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
                         std::is_same<T, Eigen::Map< Eigen::MatrixXd >>::value || 
                         std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                         std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
                         "Error - type not allowed");
           
           static_assert(std::is_same<U, Eigen::MatrixXd >::value || 
                         std::is_same<U, Eigen::Map< Eigen::MatrixXd >>::value || 
                         std::is_same<U, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                         std::is_same<U, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
                         "Error - type not allowed");
           
           Eigen::MatrixXd A = X,
                           B = Y;
           
           // if(transpX == true){ A = X.transpose(); }
           // if(transpY == true){ B = Y.transpose(); }
           
           // int chunks;//, tid;
           hsize_t block_size;
           
           std::vector<hsize_t> vsizetoReadN, vstartN,
                                vsizetoReadM, vstartM,
                                vsizetoReadK, vstartK;
           
           unsigned int ithreads;
           hsize_t M = A.rows();
           hsize_t K = A.cols();
           hsize_t N = B.cols();
           
           if( K == B.rows()) {
               
               if( iblock_size.isNotNull()) {
                   block_size =  Rcpp::as<int>(iblock_size);  
               } else {
                   block_size =  MAXBLOCKSIZE/3;  
               }
               
               C = Eigen::MatrixXd::Zero(M,N) ;
               if(block_size > std::min( N, std::min(M,K)) )
                   block_size = std::min( N, std::min(M,K)); 
               
               getBlockPositionsSizes( N, block_size, vstartN, vsizetoReadN );
               getBlockPositionsSizes( M, block_size, vstartM, vsizetoReadM );
               getBlockPositionsSizes( K, block_size, vstartK, vsizetoReadK );
               
               for (int ii = 0; ii < vstartM.size(); ii++)
               {
                   for (int jj = 0; jj < vstartN.size(); jj++)
                   {
                       for(int kk = 0; kk < vstartK.size(); kk++)
                       {
                           C.block(vstartM[ii], vstartN[jj], vsizetoReadM[ii], vsizetoReadN[jj]) =  C.block(vstartM[ii], vstartN[jj], vsizetoReadM[ii], vsizetoReadN[jj]) + 
                               (A.block(vstartM[ii], vstartK[kk], vsizetoReadM[ii], vsizetoReadK[kk]) * B.block(vstartK[kk], vstartN[jj], vsizetoReadK[kk], vsizetoReadN[jj]));
                       }
                   }
               }
               
           } else {
               throw std::range_error("non-conformable arguments");
               }
           
           
       } catch(std::exception &ex) {
           Rcpp::Rcout<<"c++ error : Bblock_matrix_mul : " <<ex.what();
           return(C);
           
       } catch(...) { 
           ::Rf_error("c++ exception in Bblock_matrix_mul (unknown reason)"); 
       }
       
       return(C);
   }
    
    
    
    
    
    
    
    // In-memory execution - Parallel version - by Blocks
    template<typename T, typename U>
    extern inline Eigen::MatrixXd Rcpp_block_matrix_mul_parallel( T X, U Y, 
                                                                  bool transpX,
                                                                  bool transpY,
                                                                  Rcpp::Nullable<int>  iblock_size, 
                                                                  Rcpp::Nullable<int> threads  = R_NilValue)
    {
        
        Eigen::MatrixXd C;
        try{
            static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
                          std::is_same<T, Eigen::Map< Eigen::MatrixXd >>::value || 
                          std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                          std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
                          "Error - type not allowed");
            
            static_assert(std::is_same<U, Eigen::MatrixXd >::value || 
                          std::is_same<U, Eigen::Map< Eigen::MatrixXd >>::value || 
                          std::is_same<U, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                          std::is_same<U, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
                          "Error - type not allowed");
            
            Eigen::MatrixXd A = X,
                            B = Y;
            
            if(transpX == true){ A = X.transpose(); }
            if(transpY == true){ B = Y.transpose(); }
            
            // int chunks;//, tid;
            hsize_t block_size;
            
            std::vector<hsize_t> vsizetoReadN, vstartN,
                                 vsizetoReadM, vstartM,
                                 vsizetoReadK, vstartK;
            
            unsigned int ithreads;
            hsize_t M = A.rows();
            hsize_t K = A.cols();
            hsize_t N = B.cols();
            
            if( iblock_size.isNotNull()) {
                block_size =  Rcpp::as<int>(iblock_size);  
            } else {
                block_size =  MAXBLOCKSIZE/3;  
            }
            
            C = Eigen::MatrixXd::Zero(M,N) ;
            if(block_size > std::min( N, std::min(M,K)) )
                block_size = std::min( N, std::min(M,K)); 
            
            ithreads = get_number_threads(threads, R_NilValue);
            
            getBlockPositionsSizes( N, block_size, vstartN, vsizetoReadN );
            getBlockPositionsSizes( M, block_size, vstartM, vsizetoReadM );
            getBlockPositionsSizes( K, block_size, vstartK, vsizetoReadK );
            
            // 
            // 
            // Rcpp::Rcout<<"\nInicis i posicions de M: ";
            // for (hsize_t ii = 0; ii <vstartM.size(); ii ++) {
            //     Rcpp::Rcout<<"\n\t"<<vstartM[ii]<<" - "<<vsizetoReadM[ii];
            // }
            // 
            // Rcpp::Rcout<<"\nInicis i posicions de N: ";
            // for (hsize_t ii = 0; ii <vstartN.size(); ii ++) {
            //     Rcpp::Rcout<<"\n\t"<<vstartN[ii]<<" - "<<vsizetoReadN[ii];
            // }
            // 
            // Rcpp::Rcout<<"\nInicis i posicions de K: ";
            // for (hsize_t ii = 0; ii <vstartK.size(); ii ++) {
            //     Rcpp::Rcout<<"\n\t"<<vstartK[ii]<<" - "<<vsizetoReadK[ii];
            // }
            // 
            
            // chunks = vstartM.size()/ithreads;
            
            #pragma omp parallel num_threads(ithreads) shared(A, B, C)// chunks) // private(tid ) 
            {
                
            #pragma omp for schedule (static) // collapse(3)
                for (int ii = 0; ii < vstartM.size(); ii++)
                {
                    for (int jj = 0; jj < vstartN.size(); jj++)
                    {
                        for(int kk = 0; kk < vstartK.size(); kk++)
                        {
                            C.block(vstartM[ii], vstartN[jj], vsizetoReadM[ii], vsizetoReadN[jj]) =  C.block(vstartM[ii], vstartN[jj], vsizetoReadM[ii], vsizetoReadN[jj]) + 
                                (A.block(vstartM[ii], vstartK[kk], vsizetoReadM[ii], vsizetoReadK[kk]) * B.block(vstartK[kk], vstartN[jj], vsizetoReadK[kk], vsizetoReadN[jj]));
                        }
                    }
                }
            }

        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception Rcpp_block_matrix_mul_parallel: "<<ex.what()<< " \n";
        }
        
        return(C);
    }
    
    template< typename T, typename U>
    extern inline Rcpp::RObject Rcpp_matrix_vect_mult ( T  A, U  B)
    {
        
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(B);
        
        if( v.length() == m.rows()) {
            
            Rcpp::NumericMatrix C = Rcpp::no_init( m.rows(), m.cols());
            
            for( int i=0; i<m.cols(); i++) {
                C( Rcpp::_, i) = m( Rcpp::_, i) * v;  
            }    
            return(C);
            
        } else if( v.length() == m.cols()) {
            
            Rcpp::NumericMatrix C = Rcpp::no_init( m.rows(), m.cols());
            
            for( int i=0; i<m.rows(); i++) {
                C( i, Rcpp::_) = m( i, Rcpp::_) * v;  
            }    
            return(C);
            
        } else {
            Rcpp::Rcout<<"Error: non-conformable arguments";
        }
        
        return(R_NilValue);
    }
    
    
    template< typename T>
    extern inline Rcpp::RObject Rcpp_vector_mult ( T  A, T  B)
    {
        
        Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(A);
        Rcpp::NumericVector v2 = Rcpp::as<Rcpp::NumericVector>(B);
        
        if(v.size() == v2.size()) {
            Rcpp::NumericVector C = Rcpp::no_init( v.size());
            
            std::transform (v.begin(), v.end(), v2.begin(), C.begin(), std::multiplies<double>());
            
            C.attr("dim") = Rcpp::Dimension( C.size(), 1); 
            
            return(C);
        }
        
        return(R_NilValue);
        
    }
    
    
    template< typename T>
    extern inline Rcpp::RObject Rcpp_matrix_vector_blockMult( T  A, T  B, Rcpp::Nullable<bool> bparal, 
                            Rcpp::Nullable<int> iblock_size, Rcpp::Nullable<int> threads)
    {
        
        // NOTA: Per defecte, multiplica per columnes tal i com raja.... 

        bool btransposed = false;
        unsigned int ithreads;
        hsize_t block_size;
        int chunks;
        
        Rcpp::NumericMatrix X = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericVector Y = Rcpp::as<Rcpp::NumericVector>(B);
        Rcpp::NumericMatrix C;
        
        // Matrix
        hsize_t M = X.rows(), N = X.cols();
        // Vector
        hsize_t K = Y.length();
        
        try {
            
            if( K==N || K==M) {
                if ( K == N){
                    // multiplies vector to every col
                    btransposed = true;
                    
                    X = Rcpp::transpose(X);
                    
                    N = X.rows();
                    M = X.cols();
                } 
                
                std::vector<hsize_t> vsizetoRead;
                std::vector<hsize_t> vstart;
                
                ithreads = get_number_threads(threads, bparal);
                
                C = Rcpp::no_init( M, N);
                
                if( iblock_size.isNotNull()) {
                    block_size =  Rcpp::as<int>(iblock_size);  
                } else {
                    block_size = getMatrixBlockSize( N, M).at(0);
                }
                
                // minimum block size: 2 columns
                if(block_size <= 0 ) {
                    block_size = M*2;
                }
                
                // Mínimum block size: 2 columns
                getBlockPositionsSizes( M*N, block_size, vstart, vsizetoRead );
                
                chunks = vstart.size()/ithreads;
                
                #pragma omp parallel num_threads(ithreads) shared(A, B, C, chunks)
                {
                #pragma omp for schedule (static)
                    for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                    {
                        // Duplicate vector
                        std::size_t const no_of_duplicates = vsizetoRead[ii] / Y.length();
                        
                        std::vector<double> v = Rcpp::as<std::vector<double> >(Y); 
                        v.reserve(Y.size() * no_of_duplicates);
                        auto end = std::end(v);
                        
                        for(std::size_t i = 1; i < no_of_duplicates; ++i)
                            v.insert(std::end(v), std::begin(v), end);
                        
                        // Mult vector to matrix by columns / rows
                        if( vstart[ii] + vsizetoRead[ii] >= M*N ) {
                            std::transform (X.begin() + vstart[ii], X.end(),
                                            v.begin(), C.begin() + vstart[ii], std::multiplies<double>());
                        } else {
                            std::transform (X.begin() + vstart[ii], X.begin() + vstart[ii] + vsizetoRead[ii],
                                            v.begin() , C.begin() + vstart[ii], std::multiplies<double>());   
                        }
                    }
                }

            } else {
                
                Rcpp::Rcout<< "vector sum error: non-conformable arguments\n";
                return(R_NilValue);
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception Rcpp_matrix_vector_blockMult: "<<ex.what()<< " \n";
            return(R_NilValue);
        }
        
        if(btransposed == true){
            Rcpp::transpose(C);
        } 
        
        C.attr("dim") = Rcpp::Dimension( M, N);
        
        return(C);
    }
    

}

#endif // BIGDATASTATMETH_ALGEBRA_MEM_MULTIPLICATION_HPP