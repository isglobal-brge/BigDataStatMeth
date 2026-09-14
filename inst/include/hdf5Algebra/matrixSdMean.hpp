/**
 * @file matrixSdMean.hpp
 * @brief Mean and standard deviation calculations for HDF5 matrices
 * @details This header file provides implementations for computing mean and
 * standard deviation statistics for matrices stored in HDF5 format. The
 * implementation includes:
 * 
 * Key features:
 * - Row-wise statistics
 * - Column-wise statistics
 * - Block-based computation
 * - Memory-efficient algorithms
 * - Parallel processing support
 * 
 * Supported operations:
 * - Mean calculation
 * - Standard deviation calculation
 * - Corrected standard deviation
 * - Block-based processing
 * - Large matrix support
 * 
 * Performance features:
 * - Cache-friendly algorithms
 * - Dynamic block sizing
 * - Multi-threaded processing
 * - I/O optimization
 * - Memory management
 * 
 * The implementation uses:
 * - Efficient statistical algorithms
 * - Block-based computation
 * - HDF5 chunked storage
 * - Parallel I/O
 * - Vectorized operations
 */

#ifndef BIGDATASTATMETH_HDF5_MATRIXSDMEAN_HPP
#define BIGDATASTATMETH_HDF5_MATRIXSDMEAN_HPP


#include <RcppEigen.h>
#include "H5Cpp.h"
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace BigDataStatMeth {

/**
 * @brief Maximum number of example indices listed in a zero-variance report
 */
const int ZERO_VARIANCE_MAX_EXAMPLES = 5;

/**
 * @brief Build the message reported for constant (zero-variance) rows/columns
 * @details Wording follows the precedent set by base R's prcomp(), which stops
 * with "cannot rescale a constant/zero column to unit variance", and adds the
 * information a user needs to fix the input: how many rows/columns are
 * affected and a few of their positions.
 *
 * @param examples Positions of the offending rows/columns, 1-based, as seen
 *        from R, at most ZERO_VARIANCE_MAX_EXAMPLES of them
 * @param nzero Total number of offending rows/columns
 * @param axis Noun used in the message, "column" or "row"
 * @return The complete message
 */
inline std::string zeroVarianceMessage( const std::vector<hsize_t>& examples,
                                        hsize_t nzero,
                                        const std::string& axis )
{
    std::ostringstream msg;

    msg << "cannot rescale a constant/zero " << axis << " to unit variance: "
        << nzero << " " << axis << (nzero == 1 ? " has" : "s have")
        << " zero or non-finite variance (";

    if( (hsize_t)examples.size() < nzero ) { msg << "e.g. "; }
    msg << axis << (examples.size() == 1 ? " " : "s ");

    for( std::size_t i = 0; i < examples.size(); i++ ) {
        if( i > 0 ) { msg << ", "; }
        msg << examples[i];
    }
    msg << ")";

    return msg.str();
}

/**
 * @brief Guard a vector of standard deviations against zero variance
 * @details Throws when any standard deviation is exactly zero or not finite,
 * i.e. when scaling would divide by zero. Used by the SVD/PCA flow, where a
 * silent 0/0 turns the whole decomposition into NaN and every singular value
 * is reported as 0 with no indication of what went wrong.
 *
 * @param sd Standard deviations, one per row/column of the matrix as seen
 *        from R, in R order
 * @param axis Noun used in the message, "column" or "row"
 * @throws std::runtime_error when at least one value is zero or not finite
 */
inline void checkNonZeroVarianceSd( const Eigen::RowVectorXd& sd,
                                    const std::string& axis = "column" )
{
    std::vector<hsize_t> examples;
    hsize_t nzero = 0;

    for( Eigen::Index i = 0; i < sd.size(); i++ ) {
        if( !std::isfinite(sd(i)) || sd(i) == 0.0 ) {
            nzero++;
            if( (int)examples.size() < ZERO_VARIANCE_MAX_EXAMPLES ) {
                examples.push_back( (hsize_t)(i + 1) );   // 1-based, R convention
            }
        }
    }

    if( nzero > 0 ) {
        throw std::runtime_error( zeroVarianceMessage(examples, nzero, axis) );
    }

    return void();
}

/**
 * @brief Guard pre-computed mean/sd statistics against zero variance
 * @details Convenience overload for the 2 x n layout produced by
 * get_HDF5_mean_sd_by_row() and get_HDF5_mean_sd_by_column(): row 0 holds the
 * means and row 1 the standard deviations.
 *
 * @param normalize 2 x n matrix of pre-computed statistics
 * @param axis Noun used in the message, "column" or "row"
 * @throws std::runtime_error when at least one standard deviation is zero or
 *         not finite
 */
inline void checkNonZeroVariance( const Eigen::MatrixXd& normalize,
                                  const std::string& axis = "column" )
{
    if( normalize.rows() < 2 ) { return void(); }
    checkNonZeroVarianceSd( normalize.row(1), axis );
    return void();
}

/**
 * @brief Guard an in-memory matrix against constant columns before scaling
 * @details Recomputes the column standard deviations exactly as
 * RcppNormalizeColwise() does - centred when @p bcenter is true, root mean
 * square otherwise, which is also what base::scale() does - so the guard fires
 * precisely when the subsequent division would be by zero.
 *
 * @tparam M Eigen matrix or mapped matrix type
 * @param X Matrix whose columns are about to be scaled
 * @param bcenter Whether the data will also be centred
 * @param axis Noun used in the message, "column" or "row", i.e. what a column
 *        of @p X is from the point of view of the R caller
 * @throws std::runtime_error when at least one column has zero or non-finite
 *         standard deviation
 */
template< typename M>
inline void checkNonZeroVarianceColwise( const M& X, bool bcenter,
                                         const std::string& axis = "column" )
{
    if( X.rows() < 2 ) { return void(); }   // sd undefined, nothing to assert

    Eigen::RowVectorXd sd;

    if( bcenter ) {
        Eigen::RowVectorXd mean = X.colwise().mean();
        sd = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
    } else {
        sd = (X.array().square().colwise().sum() / (X.rows() - 1)).sqrt();
    }

    checkNonZeroVarianceSd( sd, axis );

    return void();
}

/**
 * @brief Calculate optimal block size for processing
 * @details Determines the optimal block size for processing based on matrix
 * dimensions and memory constraints.
 *
 * @param wsize User-specified block size (optional)
 * @param reference_size Primary dimension size
 * @param alternative_size Secondary dimension size
 * @return Optimal block size for processing
 */
inline hsize_t get_block_size( Rcpp::Nullable<int> wsize, hsize_t reference_size, hsize_t alternative_size) {
    
    hsize_t bsize = 0;
    
    if( wsize.isNull()) {
        if( reference_size > MAXELEMSINBLOCK ){
            bsize = 1;
        } else {
            hsize_t maxsize = std::max( alternative_size, reference_size);
            bsize = std::ceil( MAXELEMSINBLOCK / maxsize);
        }
    } else {
        //.. 20260224 ..// if(reference_size > MAXELEMSINBLOCK){
        //.. 20260224 ..//     bsize = 1;
        //.. 20260224 ..// } else {
        //.. 20260224 ..//     bsize = Rcpp::as<int> (wsize);
        //.. 20260224 ..// }
        
        bsize = Rcpp::as<int>(wsize); 
    }
    
    return(bsize);
    
}

/**
 * @brief Calculate row-wise mean and standard deviation
 * @details Computes mean and standard deviation for each row of the matrix
 * using block-based processing for memory efficiency.
 * 
 * @param dsA Input matrix dataset
 * @param normalize Output matrix for mean and std values
 * @param wsize Block size for processing
 */
inline void get_HDF5_mean_sd_by_row( BigDataStatMeth::hdf5Dataset* dsA, 
                                     Eigen::MatrixXd& normalize, 
                                     bool bsd, bool bmean, 
                                     Rcpp::Nullable<int> wsize )
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
        //.. 20260224 ..// if( block_size < dims_out[1] ) {
        //.. 20260224 ..//     count[1] = block_size;
        //.. 20260224 ..// } else{
        //.. 20260224 ..//     count[1] = dims_out[1];
        //.. 20260224 ..// }

        // Read data in blocks of 500 columns
        //.. 20260224 ..// for( hsize_t i=0; (i <= floor(dims_out[1]/block_size)) || i==0 ; i++)
        for( hsize_t i=0; offset[1] < dims_out[1]; i++)
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
            normalize.block( 0, offset[1], 1, mean.size()) = mean;
            
            if(bsd){
                Eigen::RowVectorXd sd = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
                normalize.block( 1, offset[1], 1, sd.size()) = sd;
            }
            
            
            
            offset[1] = offset[1] + block_size;

        }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        // error.printErrorStack();
        // checkClose_file(dsA);
        throw std::runtime_error("get_HDF5_mean_sd_by_row (File IException)");
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        // error.printErrorStack();
        // checkClose_file(dsA);
        throw std::runtime_error("get_HDF5_mean_sd_by_row (DataSet IException)");
    } catch(std::exception& error) {
        // checkClose_file(dsA);
        throw std::runtime_error(std::string("get_HDF5_mean_sd_by_row function: ") + error.what());
    }
    
    return void(); 
    
}

/**
 * @brief Calculate column-wise mean and standard deviation
 * @details Computes mean and standard deviation for each column of the matrix
 * using block-based processing for memory efficiency. Optimized for cases
 * where n << m (rows much fewer than columns).
 * 
 * @param dsA Input matrix dataset
 * @param normalize Output matrix for mean and std values
 * @param bsd compute sd
 * @param bmean compute mean
 * @param wsize Block size for processing
 */
inline void get_HDF5_mean_sd_by_column( BigDataStatMeth::hdf5Dataset* dsA,
                                        Eigen::MatrixXd& normalize, 
                                        bool bsd, bool bmean, 
                                        Rcpp::Nullable<int> wsize )
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
        //.. 20260224 ..// if( block_size < dims_out[0] )
        //.. 20260224 ..//     count[0] = block_size;
        //.. 20260224 ..//  else
        //.. 20260224 ..//      count[0] = dims_out[0];
        
        // Read data in blocks of 500 columns
        //.. 20260224 ..// for(hsize_t i=0; (i <= floor(dims_out[0]/block_size)) || i==0; i++)
        for(hsize_t i=0; offset[0] < dims_out[0]; i++)
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
            normalize.block( 0, offset[0], 1, mean.size()) = mean.transpose();
            
            if(bsd) {
                Eigen::VectorXd sd = ((X.colwise() - mean).array().square().rowwise().sum() / (X.cols() - 1)).sqrt();
                normalize.block( 1, offset[0], 1, sd.size()) = sd.transpose();
            }
            
            offset[0] = offset[0] + block_size;

        }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        // error.printErrorStack();
        // checkClose_file(dsA);
        throw std::runtime_error("get_HDF5_mean_sd_by_column (File IException)");
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        // error.printErrorStack();
        // checkClose_file(dsA);
        throw std::runtime_error("get_HDF5_mean_sd_by_column (DataSet IException)");
    } catch(std::exception& error) {
        // checkClose_file(dsA);
        throw std::runtime_error(std::string("get_HDF5_mean_sd_by_column function: ") + error.what());
        // return void();
    }
    
    return void();  // successfully terminated
    
}


}

#endif // BIGDATASTATMETH_HDF5_MATRIXSDMEAN_HPP

