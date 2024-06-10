#ifndef BIGDATASTATMETH_UTIL_SPLIT_DATASETS_HPP
#define BIGDATASTATMETH_UTIL_SPLIT_DATASETS_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include "memAlgebra/memMultiplication.hpp"
// #include <thread>

namespace BigDataStatMeth {


// Internal call 
void RcppSplit_matrix_hdf5 ( BigDataStatMeth::hdf5Dataset* dstosplit, bool bycols, 
                            std::string stroutgroup, std::string stroutdataset, 
                            int blocksize, int irows, int icols )
{
    
    try {
        
        int blocks;
        hsize_t inrows = irows, 
            incols = icols,
            ii = 0,
            kk = 0;
        
        std::vector<hsize_t> stride = {1, 1},
            block = {1, 1};
        
        std::string newDatasetName = "";
        
        if( bycols == true ) {
            blocks = (icols + blocksize - 1) / blocksize;
            incols = blocksize;
        } else {
            blocks = (irows + blocksize - 1) / blocksize;
            inrows = blocksize;
        }
        
        for ( int i=0; i<blocks; i++)
        {
            // newDatasetName = stroutgroup + "/" + stroutdataset + "." + std::to_string(i), incols;
            newDatasetName = stroutgroup + "/" + stroutdataset + "." + std::to_string(i);
            
            if( bycols == true) { 
                kk = i * blocksize;
                if( kk + blocksize > icols) { incols = icols - kk; }
            } else 
            {
                ii = i * blocksize;
                if( ii + blocksize > irows) { inrows = irows - ii; }
            }
            
            // Get block from complete matrix
            // Eigen::MatrixXd Block = GetCurrentBlock_hdf5( file, dstosplit, kk, ii, incols, inrows);
            // Transform to rowmajor
            // Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mapBlock(Block.transpose().data(), inrows, incols);
            
            std::vector<double> vdts( inrows * incols );
            dstosplit->readDatasetBlock( {kk, ii}, {incols, inrows}, stride, block, vdts.data() );
            //.. REALMENT HO NECESSITEM ?? .. // Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdts.data(), inrows, incols );

            
            // ESTIC AQUÍ !!!
            // HE DE CREAR ELS RESULTATS AMB ELS SPLITS
            // FALTARIA ESCRIURE EL BLOCK AL FITXER + REVISAR ELS CATCH
            // I PODRIA SER QUE REVISAR-HO ABSOLUTAMENT TOT PERQUÈ DE 
            // SEGUR QUE ESTÀ PROGRAMAT AMB EL CUL !!!
            
            BigDataStatMeth::hdf5Dataset* dsOut = new BigDataStatMeth::hdf5Dataset(dstosplit->getFileName(), stroutgroup, stroutdataset, false);
            dsOut->createDataset( inrows, incols, "real"); 
            dsOut->writeDataset(vdts.data());
            
            // write_HDF5_matrix_from_R_ptr(file, newDatasetName, Rcpp::wrap(mapBlock), false);
            
        }
        
    }catch( H5::FileIException& error ) {
        ::Rf_error( "c++ exception (File IException )" );
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the dstosplit operations
        ::Rf_error( "c++ exception (dstosplit IException )" );
        return void();
    } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception (DataSpace IException )" );
        return void();
    } 
    
    return void();
    
    
}

}

#endif // BIGDATASTATMETH_UTIL_SPLIT_DATASETS_HPP