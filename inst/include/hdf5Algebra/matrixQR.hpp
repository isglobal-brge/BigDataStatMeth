#ifndef BIGDATASTATMETH_HDF5_MATRIXQR_HPP
#define BIGDATASTATMETH_HDF5_MATRIXQR_HPP

// #include <RcppEigen.h>
// #include "H5Cpp.h"

#include "Utilities/openme-utils.hpp"
#include "Utilities/Utilities.hpp"
#include "hdf5Algebra/multiplication.hpp"

namespace BigDataStatMeth {

// Function declaration
template< typename M> extern inline strQR RcppQR ( M X, bool bthin );
extern inline void RcppQRHdf5( BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsQ, BigDataStatMeth::hdf5Dataset* dsR, 
                               bool bthin, Rcpp::Nullable<int> block_size, Rcpp::Nullable<int> threads );

// Function definition
template< typename M>
extern inline strQR RcppQR ( M X, bool bthin )
{
    
    static_assert(std::is_same<M, Eigen::MatrixXd >::value || 
                  std::is_same<M, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                  std::is_same<M, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
                  "Error - type not allowed");
    
    Eigen::MatrixXd A = X;
    int m = A.rows(), n = A.cols();
    int irank;
    strQR vQR;
    
    Eigen::MatrixXd R;
    Eigen::MatrixXd Q;
    
    Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    
    qr.compute(A);
    irank = lu_decomp.rank();
    
    if (irank == m + 1 || irank == n + 1 )
    {
        vQR.R = qr.matrixQR().template triangularView<Eigen::Upper>();
    } else {
        vQR.R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>(); 
    }
    
    if (bthin == false)
    {
        vQR.Q =  qr.householderQ();       // Full decomposition
    } else {
        
        vQR.Q = Eigen::MatrixXd::Identity(m,n);
        vQR.Q = qr.householderQ() * vQR.Q;    // Thin decomposition
    }
    
    return(vQR);
}




// extern strQR RcppQR( Eigen::MatrixXd & A, bool bthin)
// {
//     
//     int m = A.rows(), n = A.cols();
//     int irank;
//     strQR vQR;
//     
//     Eigen::MatrixXd R;
//     Eigen::MatrixXd Q;
//     
//     Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
//     Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
//     
//     qr.compute(A);
//     irank = lu_decomp.rank();
//     
//     if (irank == m + 1 || irank == n + 1 )
//     {
//         vQR.R = qr.matrixQR().template triangularView<Eigen::Upper>();
//     } else {
//         vQR.R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>(); 
//     }
//     
//     if (bthin == false)
//     {
//         vQR.Q =  qr.householderQ();       // Full decomposition
//     } else {
//         
//         vQR.Q = Eigen::MatrixXd::Identity(m,n);
//         vQR.Q = qr.householderQ() * vQR.Q;    // Thin decomposition
//     }
//     
//     return(vQR);
//     
// }



extern inline void RcppQRHdf5( BigDataStatMeth::hdf5Dataset* dsA, 
                        BigDataStatMeth::hdf5Dataset* dsQ, 
                        BigDataStatMeth::hdf5Dataset* dsR, 
                        bool bthin,
                        Rcpp::Nullable<int> block_size = R_NilValue,
                        Rcpp::Nullable<int> threads = R_NilValue )
{
    
    
    try {

        int irank;
        std::vector<hsize_t> offset = {0,0},
            count = {dsA->nrows(), dsA->ncols()},
            stride = {1,1},
            block = {1,1};
        Eigen::MatrixXd Q;
        
        
        std::vector<double> vdA( count[0] * count[1] ); 
        dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
        Eigen::MatrixXd A = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdA.data(), count[0], count[1] );
        A.transposeInPlace();
        
        
        Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
        
        qr.compute(A);
        irank = lu_decomp.rank();
        
        Rcpp::Rcout<<"\nA veure on estem.... - 1";
        
        if (irank == count[0] + 1 || irank == count[1] + 1 )
        {
            Rcpp::Rcout<<"\nA veure on estem.... - 2.1";
            Eigen::MatrixXd R = qr.matrixQR().template triangularView<Eigen::Upper>();
            dsR->createDataset( R.rows(), R.cols(), "real" );
            dsR->writeDataset(Rcpp::wrap(R));
        } else {
            Rcpp::Rcout<<"\nA veure on estem.... - 2.2";
            Eigen::MatrixXd R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>();
            dsR->createDataset( R.rows(), R.cols(), "real" );
            dsR->writeDataset(Rcpp::wrap(R));
        }
        
        Rcpp::Rcout<<"\nA veure on estem.... - 3";
        
        if (bthin == false)
        {
            Rcpp::Rcout<<"\nA veure on estem.... - 4.1";
            //..// vQR.Q =  qr.householderQ();       // Full decomposition
            Eigen::MatrixXd Q = qr.householderQ();
            dsQ->createDataset( Q.rows(), Q.cols(), "real" );
            dsQ->writeDataset( Rcpp::wrap(Q) );
        } else {
            Rcpp::Rcout<<"\nA veure on estem.... - 4.2";
            int iblock_size = BigDataStatMeth::getMaxBlockSize( qr.householderQ().rows() , qr.householderQ().cols(), count[0], count[1], block_size);
            
            Eigen::MatrixXd Qthin = qr.householderQ() * Eigen::MatrixXd::Identity(count[0], count[1]);
            dsQ->createDataset( Qthin.rows(), Qthin.cols(), "real" );
            dsQ->writeDataset( Rcpp::wrap(Qthin));
            
        }
        
        Rcpp::Rcout<<"\nA veure on estem.... - 5";
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception RcppQRHdf5 (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception RcppQRHdf5 (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception RcppQRHdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception RcppQRHdf5" << ex.what();
        return void();
    }
    
    Rcpp::Rcout<<"\nA veure on estem.... - 6";
    return void();
    
}



}

#endif // BIGDATASTATMETH_HDF5_MATRIXQR_HPP