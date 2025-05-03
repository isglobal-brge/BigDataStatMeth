#ifndef BIGDATASTATMETH_ALGEBRA_MSUBSTRACT_HPP
#define BIGDATASTATMETH_ALGEBRA_MSUBSTRACT_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include "memAlgebra/memOtherFunctions.hpp"
// #include <thread>

namespace BigDataStatMeth {


    
    template< typename T>  extern inline Rcpp::RObject Rcpp_matrix_substract ( T  A, T  B);
    template< typename T, typename U>  extern inline Rcpp::RObject Rcpp_matrix_vect_substract ( T  A, U  B);
    template< typename T>  extern inline Rcpp::RObject Rcpp_vector_substract ( T  A, T  B);
    
    template< typename T>  extern inline Rcpp::RObject Rcpp_matrix_blockSubstract ( T  A, T  B, Rcpp::Nullable<int> threads = R_NilValue);
    template< typename T>  extern inline Rcpp::RObject Rcpp_matrix_vector_blockSubstract( T  A, T  B, Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> threads);
    
    template< typename T>
    extern inline Eigen::MatrixXd Rcpp_block_matrix_vector_substract( T  A, T  B, hsize_t block_size, 
                                                                bool bparal, Rcpp::Nullable<int> threads = R_NilValue);
    
    
    
    // extern inline void getBlockPositionsSizes( hsize_t maxPosition, hsize_t blockSize, std::vector<hsize_t>& starts, std::vector<hsize_t>& sizes ){
    // 
    //     hsize_t isize = blockSize + 1;
    // 
    //     for (hsize_t ii = 0; ii < maxPosition; ii += blockSize)
    //     {
    //         if( ii + blockSize > maxPosition ) {
    //             isize = maxPosition - ii; }
    // 
    //         hsize_t sizetoRead = getOptimBlockSize( maxPosition, blockSize, ii, isize);
    // 
    //         starts.push_back(ii);
    //         sizes.push_back(sizetoRead);
    // 
    //         if( ii + blockSize > maxPosition ) isize = blockSize + 1;
    //         if( sizetoRead > blockSize ) {
    //             ii = ii - blockSize + sizetoRead; }
    //     }
    // 
    // }

    
    template< typename T>
    extern inline Rcpp::RObject Rcpp_matrix_substract ( T  A, T  B)
    {
        
        
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericMatrix m2 = Rcpp::as<Rcpp::NumericMatrix>(B);
        
        if( m.rows() == m2.rows() && m.cols() == m2.cols()) {
            Rcpp::NumericVector C = m - m2;
            // C.attr("dim") = Rcpp::Dimension( m.rows(), m.cols());
            
            return(C);
            
        } else {
            Rcpp::Rcout<<"Error: non-conformable arguments";
        }
        
        return(R_NilValue);
        
    }
    
    
    // Resta per files o columnes depenent si la mida del vector es igual al nombre
    // de files o igual al nombre de columnes
    template< typename T, typename U>
    extern inline Rcpp::RObject Rcpp_matrix_vect_substract ( T  A, U  B)
    {
        
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(B);
        
        if( v.length() == m.rows()) {
            
            Rcpp::NumericMatrix C = Rcpp::no_init( m.rows(), m.cols());
            
            for( int i=0; i<m.cols(); i++) {
                C( Rcpp::_, i) = m( Rcpp::_, i) - v;  
            }    
            return(C);
            
        } else if( v.length() == m.cols()) {
            
            Rcpp::NumericMatrix C = Rcpp::no_init( m.rows(), m.cols());
            
            for( int i=0; i<m.rows(); i++) {
                C( i, Rcpp::_) = m( i, Rcpp::_) - v;  
            }    
            return(C);
            
        } else {
            Rcpp::Rcout<<"Error: non-conformable arguments";
        }
        
        return(R_NilValue);
    }
    
    
    template< typename T>
    extern inline Rcpp::RObject Rcpp_vector_substract ( T  A, T  B)
    {
        
        Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(A);
        Rcpp::NumericVector v2 = Rcpp::as<Rcpp::NumericVector>(B);
        
        if(v.size() == v2.size()) {
            
            Rcpp::NumericVector C = Rcpp::no_init( v.size());
            
            std::transform (v.begin(), v.end(), v2.begin(), C.begin(), std::minus<double>());
            
            C.attr("dim") = Rcpp::Dimension( C.size(), 1); 
            
            return(C);
        }
        
        return(R_NilValue);
        
    }
    
    
    template< typename T>
    extern inline Rcpp::RObject Rcpp_matrix_blockSubstract ( T  A, T  B, Rcpp::Nullable<int> threads)
    {
        
        // static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value ,
        //               "Error - type not allowed");
        
        Rcpp::NumericMatrix X = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericMatrix Y = Rcpp::as<Rcpp::NumericMatrix>(B);
        
        hsize_t N = X.rows();
        hsize_t M = X.cols();
        
        Rcpp::NumericMatrix C = Rcpp::no_init( N, M);
        
        try {
            
            unsigned int ithreads;
            hsize_t block_size; 
            
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
                    
                    ithreads = get_number_threads(threads, R_NilValue);
                    
                    getBlockPositionsSizes( N*M, block_size, vstart, vsizetoRead );
                    int chunks = vstart.size()/ithreads;
                    
                    #pragma omp parallel num_threads(ithreads) shared(A, B, C, chunks)
                    {
                    #pragma omp for schedule (dynamic)
                        for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                        {
                            
                            if( vstart[ii] + vsizetoRead[ii] >= N*M ) {
                                std::transform (X.begin() + vstart[ii], X.end(),
                                                Y.begin() + vstart[ii], C.begin() + vstart[ii], std::minus<double>());
                            } else {
                                std::transform (X.begin() + vstart[ii], X.begin() + vstart[ii] + vsizetoRead[ii],
                                                Y.begin() + vstart[ii], C.begin() + vstart[ii], std::minus<double>());   
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
    extern inline Rcpp::RObject Rcpp_matrix_vector_blockSubstract( T  A, T  B,  
                                Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> threads)
    {
        
        // NOTA: Per defecte, suma per columnes tal i com raja.... 
        
        // static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
        //               std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
        //               "Error - type not allowed");
        
        bool btransposed = false;
        
        Rcpp::NumericMatrix X = Rcpp::as<Rcpp::NumericMatrix>(A);
        Rcpp::NumericVector Y = Rcpp::as<Rcpp::NumericVector>(B);
        Rcpp::NumericMatrix C;
        
        // Matrix
        hsize_t M = X.rows(),
                N = X.cols();
        
        // Vector
        hsize_t K = Y.length();
        
        try {
            
            unsigned int ithreads;
            hsize_t block_size;
            
            if( K==N || K==M) {
                if ( K == N){
                    // Sum vector to every col
                    btransposed = true;
                    
                    X = Rcpp::transpose(X);
                    
                    //.. Revisar-ho comentat 2024/04/06 ..// hsize_t N = X.rows();
                    //.. Revisar-ho comentat 2024/04/06 ..// hsize_t M = X.cols();
                    N = X.rows();
                    M = X.cols();
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
                #pragma omp for schedule (dynamic) // collapse(2)
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
                                            v.begin(), C.begin() + vstart[ii], std::minus<double>());
                        } else {
                            std::transform (X.begin() + vstart[ii], X.begin() + vstart[ii] + vsizetoRead[ii],
                                            v.begin() , C.begin() + vstart[ii], std::minus<double>());   
                        }
                    }
                }
    
            } else {
                
                Rcpp::Rcout<< "vector sum error: non-conformable arguments\n";
                return(R_NilValue);
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception Rcpp_matrix_vector_blockSubstract: "<<ex.what()<< " \n";
            return(R_NilValue);
        }
        
        if(btransposed == true){
            Rcpp::transpose(C);
        } 
        
        C.attr("dim") = Rcpp::Dimension( M, N);
        
        return(C);
        
    }



}

#endif // BIGDATASTATMETH_ALGEBRA_MSUBSTRACT_HPP