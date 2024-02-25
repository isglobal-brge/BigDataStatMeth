#ifndef BIGDATASTATMETH_ALGEBRA_VECTORMATRIX_HPP
#define BIGDATASTATMETH_ALGEBRA_VECTORMATRIX_HPP

#include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include <thread>

namespace BigDataStatMeth {


// by Rows

extern inline Eigen::MatrixXd Rcpp_matrixVectorMultiplication_byRow(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().colwise() * v.array();
    return(X);
}

extern inline Eigen::MatrixXd Rcpp_matrixVectorSum_byRow(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().colwise() + v.array();
    return(X);
}

extern inline Eigen::MatrixXd Rcpp_matrixVectorSubstract_byRow(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().colwise() - v.array();
    return(X);
}

extern inline Eigen::MatrixXd Rcpp_matrixVectorDivision_byRow(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().colwise() / v.array();
    return(X);
}


// by Columns

extern inline Eigen::MatrixXd Rcpp_matrixVectorMultiplication_byCol(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().rowwise() * v.transpose().array();    
    return(X);
}

extern inline Eigen::MatrixXd Rcpp_matrixVectorSum_byCol(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().rowwise() + v.transpose().array();    
    return(X);
}

extern inline Eigen::MatrixXd Rcpp_matrixVectorSubstract_byCol(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().rowwise() - v.transpose().array();    
    return(X);
}

extern inline Eigen::MatrixXd Rcpp_matrixVectorDivision_byCol(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().rowwise() / v.transpose().array();    
    return(X);
}



/* *****************
 * Function:
 *  1: '+'
 *  2: '-'
 *  3: '*'
 *  4: '/'
***************** */

extern inline BigDataStatMeth::hdf5Dataset* hdf5_matrixVector_calculus(
        BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, 
        BigDataStatMeth::hdf5Dataset* dsC, int function, bool bbyrows, 
        bool bparal, Rcpp::Nullable<int> threads  = R_NilValue)
{
    
    try{
        
        std::vector<hsize_t> stride = {1, 1};
        std::vector<hsize_t> block = {1, 1};
        
        int blocksize = 0;
        
        // Define blocksize atending number of elements in rows and cols
        if( bbyrows == false) {
            if( dsA->ncols() > MAXELEMSINBLOCK ) {
                blocksize = 1;
            } else {
                hsize_t maxsize = std::max<hsize_t>(  dsA->nrows(),  dsA->ncols());
                blocksize = std::ceil( MAXELEMSINBLOCK / maxsize);
            }
            
        } else {
            if( dsA->nrows() > MAXELEMSINBLOCK) {
                blocksize = 1;
            } else {
                hsize_t maxsize = std::max<hsize_t>( dsA->nrows(), dsA->ncols());
                blocksize = std::ceil( MAXELEMSINBLOCK / maxsize);
            }
        }
        
        
        dsC->createDataset( dsA->nrows(), dsA->ncols(), "real");    
        
        std::vector<double> vdB( dsB->nrows() * dsB->ncols());
        dsB->readDatasetBlock( {0, 0}, { dsB->nrows(), dsB->ncols()}, stride, block, vdB.data() );
        
        
        if( bbyrows == false) {
            
            Eigen::RowVectorXd vWeights = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(vdB.data(), vdB.size());
            
            for(hsize_t i=0; (i * blocksize) <= dsA->nrows(); i++)
            {
                hsize_t sizetoread = 0;
                if((i+1)*blocksize<dsA->nrows()) {
                    sizetoread = blocksize;
                } else {
                    sizetoread = dsA->nrows()-(i*blocksize);
                }
                
                std::vector<hsize_t> offset = { i*blocksize, 0};
                std::vector<hsize_t> count = {sizetoread, dsA->ncols()};
                
                // Compute and write data
                std::vector<double> vdA( count[0] * count[1]);
                dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdA.data(), count[0], count[1] );    
                                
                if( function == 0 ) {
                    X = Rcpp_matrixVectorSum_byCol(X, vWeights); 
                } else if( function == 1 ) {
                    X = Rcpp_matrixVectorSubstract_byCol(X, vWeights);
                }else if( function == 2 ) {
                    X = Rcpp_matrixVectorMultiplication_byCol(X, vWeights);
                }else if( function == 3 ) {
                    X = Rcpp_matrixVectorDivision_byCol(X, vWeights);
                }
                
                dsC->writeDatasetBlock(Rcpp::wrap(X.transpose()), offset, count, stride, block, false);
            }
        } else {
            
            Eigen::VectorXd vWeights = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(vdB.data(), vdB.size());
            
            for(hsize_t i=0; i*blocksize <= dsA->ncols() ; i++) {
                hsize_t sizetoread = 0;
                if( (i+1)*blocksize < dsA->ncols() ){
                    sizetoread = blocksize;
                } else {
                    sizetoread = dsA->ncols()-(i*blocksize);
                }
                
                std::vector<hsize_t> offset = {0, i*blocksize};
                std::vector<hsize_t> count = { dsA->nrows(), sizetoread};
                
                // Compute and write data
                std::vector<double> vdA( count[0] * count[1]);
                dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdA.data(), count[0], count[1] );
                
                if( function == 0 ) {
                    X = Rcpp_matrixVectorSum_byRow(X, vWeights); 
                } else if( function == 1 ) {
                    X = Rcpp_matrixVectorSubstract_byRow(X, vWeights);
                }else if( function == 2 ) {
                    X = Rcpp_matrixVectorMultiplication_byRow(X, vWeights);
                }else if( function == 3 ) {
                    X = Rcpp_matrixVectorDivision_byRow(X, vWeights);
                }
                
                offset = {i*blocksize, 0};
                
                dsC->writeDatasetBlock(Rcpp::wrap(X.transpose()), offset, count, stride, block, true);
            }
        }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        error.printErrorStack();
        ::Rf_error( "c++ exception hdf5_matrixVector_calculus (File IException)" );
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        error.printErrorStack();
        ::Rf_error( "c++ exception hdf5_matrixVector_calculus (DataSet IException)" );
    } catch(std::exception& error) {
        Rcpp::Rcout<< "c++ exception vector-matrix functions: "<<error.what()<< " \n";
        return(dsC);
    }
    
    return(dsC);
    
}



}

#endif // BIGDATASTATMETH_ALGEBRA_VECTORMATRIX_HPP