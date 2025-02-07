#ifndef BIGDATASTATMETH_HDF5_MATRIXSDMEAN_HPP
#define BIGDATASTATMETH_HDF5_MATRIXSDMEAN_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include "hdf5Utilities/hdf5Utilities.hpp"
#include <cmath>
// #include "H5Cpp.h"

namespace BigDataStatMeth {


// Returns de block size to be used in sd and mean calculus 
extern inline hsize_t get_block_size( Rcpp::Nullable<int> wsize, hsize_t reference_size, hsize_t alternative_size) {
    
    hsize_t bsize = 0;
    
    if( wsize.isNull()) {
        if( reference_size > MAXELEMSINBLOCK ){
            bsize = 1;
        } else {
            hsize_t maxsize = std::max( alternative_size, reference_size);
            bsize = std::ceil( MAXELEMSINBLOCK / maxsize);
        }
    } else {
        if(reference_size > MAXELEMSINBLOCK){
            bsize = 1;
        } else {
            bsize = Rcpp::as<int> (wsize);
        }
    }
    
    return(bsize);
    
}



// Get mean and corrected sd from each row in dataset 
// void get_HDF5_mean_sd_by_row_ptr( BigDataStatMeth::hdf5Dataset* dsA, Eigen::MatrixXd& normalize, hsize_t block_size )
extern inline void get_HDF5_mean_sd_by_row( BigDataStatMeth::hdf5Dataset* dsA, Eigen::MatrixXd& normalize, Rcpp::Nullable<int> wsize )
{
    
    try
    {
        
        hsize_t block_size = 0;
        hsize_t* dims_out = dsA->dim();

        std::vector<hsize_t> stride = {1, 1},
                             block = {1, 1},
                             offset = {0, 0},
                             count = {0, 0};
        
        block_size = get_block_size(wsize, dims_out[0], dims_out[1]);

        count[0] = dims_out[0];
        if( block_size < dims_out[1] ) {
            count[1] = block_size;
        } else{
            count[1] = dims_out[1];
        }

        // Read data in blocks of 500 columns
        for( hsize_t i=0; (i <= floor(dims_out[1]/block_size)) || i==0 ; i++)
        {
            
            // if( i>0 ) {
                

            if( offset[1] + block_size <= dims_out[1] ) {
                count[1] = block_size;
            } else {
                count[1] = dims_out[1] - offset[1];
            }
            // }

            std::vector<double> vdA( count[0] * count[1] ); 
            dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdA.data(), count[0], count[1] );

            Eigen::RowVectorXd mean = X.colwise().mean();
            Eigen::RowVectorXd sd = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
            
            normalize.block( 0, offset[1], 1, mean.size()) = mean;
            normalize.block( 1, offset[1], 1, sd.size()) = sd;
            
            offset[1] = offset[1] + block_size;

        }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        error.printErrorStack();
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_row (File IException)" );
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        error.printErrorStack();
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_row (DataSet IException)" );
    } catch(std::exception& error) {
        Rcpp::Rcout<< "c++ exception get_HDF5_mean_sd_by_row function: "<<error.what()<< " \n";
        return void();
    }
    
    return void(); 
    
}


// Get mean and corrected sd from each column in dataset in the case of n<<m, this information is used
// to normalize data, center or scale.
extern inline void get_HDF5_mean_sd_by_column( BigDataStatMeth::hdf5Dataset* dsA, Eigen::MatrixXd& normalize, Rcpp::Nullable<int> wsize )
{
    
    // IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
    
    try
    {

        hsize_t block_size = 0;
        hsize_t* dims_out = dsA->dim();
        
        std::vector<hsize_t> stride = {1, 1},
                             block = {1, 1},
                             offset = {0, 0},
                             count = {0, 0};
        
        
        block_size = get_block_size(wsize, dims_out[1], dims_out[0]);

        count[1] = dims_out[1];
        if( block_size < dims_out[0] )
            count[0] = block_size;
        else
            count[0] = dims_out[0];
        
        Rcpp::Rcout<<"\n--------- Estem a les modificacions??  -------------\n";
        
        // Read data in blocks of 500 columns
        for(hsize_t i=0; (i <= floor(dims_out[0]/block_size)) || i==0; i++)
        {

            if( offset[0] + block_size <= dims_out[0] ) {
                count[0] = block_size;
            }else {
                count[0] = dims_out[0] - offset[0];
            }
            
            std::vector<double> vdA( count[0] * count[1] ); 
            dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdA.data(), count[0], count[1] );

            Eigen::VectorXd mean = X.rowwise().mean();
            Eigen::VectorXd sd = ((X.colwise() - mean).array().square().rowwise().sum() / (X.cols() - 1)).sqrt();

            normalize.block( 0, offset[0], 1, mean.size()) = mean.transpose();
            normalize.block( 1, offset[0], 1, sd.size()) = sd.transpose();
            
            offset[0] = offset[0] + block_size;

        }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        error.printErrorStack();
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column (File IException)" );
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        error.printErrorStack();
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column (DataSet IException)" );
    } catch(std::exception& error) {
        Rcpp::Rcout<< "c++ exception get_HDF5_mean_sd_by_column function: "<<error.what()<< " \n";
        return void();
    }
    
    return void();  // successfully terminated
    
}


}

#endif // BIGDATASTATMETH_HDF5_MATRIXSDMEAN_HPP

