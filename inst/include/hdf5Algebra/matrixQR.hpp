#ifndef BIGDATASTATMETH_HDF5_MATRIXQR_HPP
#define BIGDATASTATMETH_HDF5_MATRIXQR_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

#include "Utilities/openme-utils.hpp"
#include "Utilities/Utilities.hpp"
#include "hdf5Algebra/multiplication.hpp"

namespace BigDataStatMeth {



extern strQR RcppQR( Eigen::MatrixXd & A, bool bthin)
{
    
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



extern void RcppQRHdf5( BigDataStatMeth::hdf5Dataset* dsA, 
                        BigDataStatMeth::hdf5Dataset* dsQ, 
                        BigDataStatMeth::hdf5Dataset* dsR, 
                        bool bthin,
                        Rcpp::Nullable<int> block_size = R_NilValue,
                        Rcpp::Nullable<int> threads = R_NilValue )
{
    
    
    try {

                
        // int m = dsA->nrows(), n = dsA->ncols();
        int irank;
        std::vector<hsize_t> offset = {0,0},
            count = {dsA->nrows(), dsA->ncols()},
            stride = {1,1},
            block = {1,1};
        // strQR vQR;
        // Eigen::MatrixXd R;
        Eigen::MatrixXd Q;
        
        
        std::vector<double> vdA( count[0] * count[1] ); 
        dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
        Eigen::MatrixXd A = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdA.data(), count[0], count[1] );
        A.transposeInPlace();
        
        // Rcpp::Rcout<<"\nMatriu A llegida:\n"<<A<<"\n";
        
        Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
        
        qr.compute(A);
        irank = lu_decomp.rank();
        
        if (irank == count[0] + 1 || irank == count[1] + 1 )
        {
            
            // vQR.R = qr.matrixQR().template triangularView<Eigen::Upper>();
            Eigen::MatrixXd R = qr.matrixQR().template triangularView<Eigen::Upper>();
            dsR->createDataset( R.rows(), R.cols(), "real" );
            dsR->writeDataset(Rcpp::wrap(R));
        } else {
            
            // vQR.R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>(); 
            Eigen::MatrixXd R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>();
            dsR->createDataset( R.rows(), R.cols(), "real" );
            dsR->writeDataset(Rcpp::wrap(R));
        }
        
        if (bthin == false)
        {
            //..// vQR.Q =  qr.householderQ();       // Full decomposition
            Eigen::MatrixXd Q = qr.householderQ();
            dsQ->createDataset( Q.rows(), Q.cols(), "real" );
            dsQ->writeDataset( Rcpp::wrap(Q) );
        } else {
            
            //..// vQR.Q = Eigen::MatrixXd::Identity(count[0], count[1]);
            //..// vQR.Q = qr.householderQ() * vQR.Q;    // Thin decomposition
            
            int iblock_size = BigDataStatMeth::getMaxBlockSize( qr.householderQ().rows() , qr.householderQ().cols(), count[0], count[1], block_size);
            
            Eigen::MatrixXd Qthin = qr.householderQ() * Eigen::MatrixXd::Identity(count[0], count[1]);
            dsQ->createDataset( Qthin.rows(), Qthin.cols(), "real" );
            dsQ->writeDataset( Rcpp::wrap(Qthin));
            
        }
        
        
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
    
    return void();
    
}



}

#endif // BIGDATASTATMETH_HDF5_MATRIXQR_HPP