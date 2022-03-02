#ifndef parallelBlockMultSparse_hdf5
#define parallelBlockMultSparse_hdf5

    #include <RcppEigen.h>
    #include <Eigen/Sparse>
    #include "BigDataStatMeth.h"
    // #include <omp.h>
    #include "pkg_omp.h"
    #include <thread>
    #include <cstdlib>
    #include <cmath>
    #include "rhdf5Utils.h"
    // #include "hdf5_to_Eigen.h"
    
    
    
    // C++ functions 
    Eigen::SparseMatrix<double> matmult_sparse_parallel ( Eigen::SparseMatrix<double>  A, 
                                                        Eigen::SparseMatrix<double> B,
                                                        Rcpp::Nullable<int> threads );
    
    Eigen::SparseMatrix<double> matmult_sparse( Eigen::Map<Eigen::SparseMatrix<double> > A, 
                                              Eigen::Map<Eigen::SparseMatrix<double> > B);
    
    // R functions
    Rcpp::RObject bdblockmult_sparse(Rcpp::RObject A, Rcpp::RObject B, 
                                 Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads );

#endif
