#ifndef BIGDATASTATMETH_ALGEBRA_SUM_HPP
#define BIGDATASTATMETH_ALGEBRA_SUM_HPP

// #include <RcppEigen.h>
// #include "Utilities/openme-utils.hpp"
// #include <thread>
#include "memAlgebra/memSum.hpp"

namespace BigDataStatMeth {


    extern inline BigDataStatMeth::hdf5Dataset*  Rcpp_block_matrix_sum_hdf5( 
            BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
            hsize_t hdf5_block, hsize_t mem_block_size, bool bparal, Rcpp::Nullable<int> threads);
    
    extern inline BigDataStatMeth::hdf5Dataset* Rcpp_block_matrix_vector_sum_hdf5( 
            BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
            hsize_t hdf5_block, bool bparal, Rcpp::Nullable<int> threads);
    
    
    
    // Working directly with C-matrix in hdf5 file
    extern inline BigDataStatMeth::hdf5Dataset*  Rcpp_block_matrix_sum_hdf5( 
            BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
            hsize_t hdf5_block, bool bparal, Rcpp::Nullable<int> threads  = R_NilValue)
    {
        
        try {
            
            hsize_t K = dsA->nrows();
            hsize_t N = dsA->ncols();
            
            if( K == dsB->nrows() && N == dsB->ncols())
            {
                
                // Parallellization and Block variables 
                unsigned int ithreads;
                std::vector<hsize_t> vstart, vsizetoRead;
                std::vector<hsize_t> stride = {1, 1};
                std::vector<hsize_t> block = {1, 1};
                hsize_t isize = hdf5_block + 1;
                
                dsC->createDataset( N, K, "real"); 
                
                // if(bparal == false) {
                //     ithreads = 1;
                // } else {
                //     if(threads.isNotNull()) {
                //         if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
                //             ithreads = Rcpp::as<int> (threads);
                //         } else {
                //             ithreads = getDTthreads(0, true);
                //         }
                //     } else {
                //         ithreads = getDTthreads(0, true);
                //     }    
                // }
                
                ithreads = get_threads(bparal, threads);
                
                if( K<=N ) {
                    
                    getBlockPositionsSizes( N, hdf5_block, vstart, vsizetoRead );
                    
                    #pragma omp parallel num_threads(ithreads) shared(dsA, dsB, dsC)
                    {
                    #pragma omp for schedule (static)
                        for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                        {
                            
                            std::vector<double> vdA( K * vsizetoRead[ii] ); 
                            // #pragma omp critical 
                            // {
                                dsA->readDatasetBlock( {0, vstart[ii]}, { K, vsizetoRead[ii]}, stride, block, vdA.data() );
                            // }
                            
                            std::vector<double> vdB( K * vsizetoRead[ii] ); 
                            // #pragma omp critical 
                            // {
                                dsB->readDatasetBlock( {0, vstart[ii]}, {K, vsizetoRead[ii]}, stride, block, vdB.data() );
                            // }
                            std::transform (vdA.begin(), vdA.end(),
                                            vdB.begin(), vdA.begin(), std::plus<double>());
                            
                            std::vector<hsize_t> offset = { 0, vstart[ii] };
                            std::vector<hsize_t> count = { K, vsizetoRead[ii] };
                            #pragma omp critical 
                            {
                                dsC->writeDatasetBlock(vdA, offset, count, stride, block);
                            }
                        }
                    }
                    
                } else {
                    
                    getBlockPositionsSizes( K, hdf5_block, vstart, vsizetoRead );
                    #pragma omp parallel num_threads(ithreads) shared(dsA, dsB, dsC)
                    {
                    #pragma omp for schedule (static)
                        for (hsize_t ii = 0; ii < vstart.size(); ii++)
                        {
                            std::vector<double> vdA( vsizetoRead[ii] * N ); 
                            // #pragma omp critical 
                            // {
                                dsA->readDatasetBlock( {vstart[ii], 0}, { vsizetoRead[ii], N}, stride, block, vdA.data() );
                            // }
                            
                            std::vector<double> vdB( vsizetoRead[ii] * N); 
                            // #pragma omp critical 
                            // {
                                dsB->readDatasetBlock( {vstart[ii], 0}, {vsizetoRead[ii], N}, stride, block, vdB.data() );
                            // }
                            
                            std::transform (vdA.begin(), vdA.end(),
                                            vdB.begin(), vdA.begin(), std::plus<double>());
                            
                            std::vector<hsize_t> offset = { vstart[ii], 0 };
                            std::vector<hsize_t> count = { vsizetoRead[ii], N };
                            #pragma omp critical 
                            {
                                dsC->writeDatasetBlock(vdA, offset, count, stride, block);
                            }
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
    extern inline BigDataStatMeth::hdf5Dataset* Rcpp_block_matrix_vector_sum_hdf5( 
            BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsC,
            hsize_t hdf5_block, bool bparal, Rcpp::Nullable<int> threads  = R_NilValue)
    {
    
        try {
            
            // Vector
            hsize_t K = dsA->nrows();
            hsize_t N = dsA->ncols();
            
            // Matrix
            hsize_t M = dsB->nrows();
            hsize_t L = dsB->ncols();
            
            std::vector<hsize_t> vstart, vsizetoRead;
            unsigned int ithreads;
    
            if(hdf5_block == 1) {
                hdf5_block = ceil(MAXBLOCKSIZE/(K*N));
            }
    
            hsize_t isize = hdf5_block + 1;
            std::vector<hsize_t> stride = {1, 1},
                                 block = {1, 1};
    
            
            dsC->createDataset( L, M, "real");
            
            std::vector<double> vdA( K * N );
            dsA->readDatasetBlock( {0, 0}, {K, N}, stride, block, vdA.data() );
            
            
            ithreads = get_threads(bparal, threads);
            
            if(  K == M )
            { // Sum vector to every col
                
                getBlockPositionsSizes( L, hdf5_block, vstart, vsizetoRead );
                
                #pragma omp parallel num_threads(ithreads) shared(dsA, dsB, dsC)
                {
                    #pragma omp for schedule (static)
                    for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                    {
                        std::vector<double> vdB( K * vsizetoRead[ii] );
                        dsB->readDatasetBlock( {0, vstart[ii]}, {K, vsizetoRead[ii]}, stride, block, vdB.data() );
                        
                        // Duplicate vector
                        std::size_t const no_of_duplicates = (K * vsizetoRead[ii]) / vdA.size();
                        
                        std::vector<double> v = vdA; 
                        v.reserve(vdA.size() * no_of_duplicates);
                        auto end = std::end(v);
                        
                        for(std::size_t i = 1; i < no_of_duplicates; ++i)
                            v.insert(std::end(v), std::begin(v), end);
                        
                        // Sum vector to matrix by columns / rows
                        std::transform (vdB.begin(), vdB.end(),
                                        v.begin(), vdB.begin(), std::plus<double>());
                        
                        
                        std::vector<hsize_t> offset = { 0, vstart[ii]};
                        std::vector<hsize_t> count = { K, vsizetoRead[ii]};
                        
                        #pragma omp critical 
                        {
                            dsC->writeDatasetBlock( vdB, offset, count, stride, block);
                        }
                    }
                }
                
    
            } else if(  K == L ) { // Sum vector to every row
                
                getBlockPositionsSizes( M, hdf5_block, vstart, vsizetoRead );
                
                #pragma omp parallel num_threads(ithreads) shared(dsA, dsB, dsC)
                {
                    #pragma omp for schedule (static)
                    for (hsize_t ii = 0; ii < vstart.size(); ii ++)
                    {
                        std::vector<double> vdB( L * vsizetoRead[ii] );
                        dsB->readDatasetBlock( {vstart[ii], 0}, {vsizetoRead[ii], K}, stride, block, vdB.data() );
                        Rcpp::NumericMatrix B (vsizetoRead[ii], K, vdB.begin());
                        
                        Rcpp::transpose(B);
                        
                        // Duplicate vector
                        std::size_t const no_of_duplicates = (vsizetoRead[ii] * K) / vdA.size();
                        
                        std::vector<double> v = vdA; 
                        v.reserve(v.size() * no_of_duplicates);
                        auto end = std::end(v);
                        
                        for(std::size_t i = 1; i < no_of_duplicates; ++i)
                            v.insert(std::end(v), std::begin(v), end);
                        
                        // Sum vector to matrix by columns / rows
                        std::transform (B.begin(), B.end(),
                                        v.begin(), B.begin() , std::plus<double>());
                    
                        std::vector<hsize_t> offset = { vstart[ii], 0};
                        std::vector<hsize_t> count = { vsizetoRead[ii], K};
                        
                        #pragma omp critical 
                        {
                            // dsC->writeDatasetBlock( Rcpp::transpose(B), offset, count, stride, block, false); 
                            dsC->writeDatasetBlock( Rcpp::as<std::vector<double> >(B), offset, count, stride, block);
                        }
                    }
                }
            } else {
                Rcpp::Rcout<< "vector sum error: non-conformable arguments\n";
            }
    
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception Sum: "<<ex.what()<< " \n";
            return(dsC);
        }
        return(dsC);
    }

}

#endif // BIGDATASTATMETH_ALGEBRA_SUM_HPP