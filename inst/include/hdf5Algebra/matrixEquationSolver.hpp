#ifndef BIGDATASTATMETH_HDF5_MATRIXEQUATIONSOLVER_HPP
#define BIGDATASTATMETH_HDF5_MATRIXEQUATIONSOLVER_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

#include "Utilities/openme-utils.hpp"
#include "Utilities/Utilities.hpp"

namespace BigDataStatMeth {

    // Symbols in the LAPACK library files : 
    // DGESV computes the solution to a real system of linear equations : A * X = B, where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    extern "C" {
        extern int dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
    }
    
    // DSYSV computes the solution to a real system of linear equations :   A * X = B, where A is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices.
    extern "C" {
        extern int dsysv_( char*, int*, int*, double*, int*, int*, double*, int*, double*, int*, int*);
    }


    extern inline Eigen::MatrixXd RcppSolve(Eigen::Map<Eigen::MatrixXd> a, Eigen::Map<Eigen::MatrixXd> b)
    {
        
        try {
            
            char Uchar = 'U';
            int info = 0;
            
            // Declare matrix variables
            int n = a.rows();
            int nrhs = b.cols();
            int lwork = std::max( 1, n );
            int lda = std::max( 1, n );
            int ldb = std::max( 1, n );
            std::vector<int> ipiv(n);
            std::vector<double> work(lwork);
            
            // Solve matrix equation
            if( a == a.transpose()  )
            {
                // dsysv_( char* UPLO, int* N , int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, double* WORK, int* LWORK, int* INFO);
                dsysv_( &Uchar, &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, work.data(), &lwork, &info);
            } else {
                
                // dgesv( int N, int NRHS, double A, int LDA, int IPIV, double B, int LDB, int INFO);
                dgesv_( &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, &info );
            }
            
        } catch(std::exception &ex) {
            Rcpp::Rcout<< ex.what();
            return(Eigen::MatrixXd::Zero(2,2));
        }
        
        return(b);
        
    }
    
    
    extern inline void RcppSolveHdf5(BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, BigDataStatMeth::hdf5Dataset* dsX )
    {
        
        try {
            
            std::vector<hsize_t> offset = {0,0},
                countA = {dsA->nrows(), dsA->ncols()},
                countB = {dsB->nrows(), dsB->ncols()},
                stride = {1,1},
                block = {1,1};
            
            std::vector<double> vdB( countB[0] * countB[1] );
            dsB->readDatasetBlock( {offset[0], offset[1]}, {countB[0], countB[1]}, stride, block, vdB.data() );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> b (vdB.data(), countB[0], countB[1]);
            
            {
                char Uchar = 'U';
                int info = 0;
                
                // Declare matrix variables
                int n = dsA->nrows();
                int nrhs = dsB->nrows();
                int lwork = std::max( 1, n );
                int lda = std::max( 1, n );
                int ldb = std::max( 1, n );
                std::vector<int> ipiv(n);
                std::vector<double> work(lwork);
                
                std::vector<double> vdA( countA[0] * countA[1] );
                dsA->readDatasetBlock( {offset[0], offset[1]}, {countA[0], countA[1]}, stride, block, vdA.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> a (vdA.data(), countA[0], countA[1]);
                
                // Solve matrix equation
                if( a == a.transpose()  )
                {
                    // dsysv_( char* UPLO, int* N , int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, double* WORK, int* LWORK, int* INFO);
                    dsysv_( &Uchar, &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, work.data(), &lwork, &info);
                } else {
                    // dgesv( int N, int NRHS, double A, int LDA, int IPIV, double B, int LDB, int INFO);
                    dgesv_( &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, &info );
                }
                
            }
            
            dsX->writeDataset(b.data());
            
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"Error in RcppSolveHdf5";
            Rcpp::Rcout<< ex.what();
            return void();
        }
        
        return void();
        
    }



}

#endif // BIGDATASTATMETH_HDF5_MATRIXEQUATIONSOLVER_HPP