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
            newDatasetName = stroutgroup + "/" + stroutdataset + "." + std::to_string(i);
            
            if( bycols == true) { 
                kk = i * blocksize;
                if( kk + blocksize > icols) { incols = icols - kk; }
            } else 
            {
                ii = i * blocksize;
                if( ii + blocksize > irows) { inrows = irows - ii; }
            }
            
            std::vector<double> vdts( inrows * incols );
            dstosplit->readDatasetBlock( {kk, ii}, {incols, inrows}, stride, block, vdts.data() );

            //. ORIGINAL A 2025/01/15.// BigDataStatMeth::hdf5Dataset* dsOut = new BigDataStatMeth::hdf5Dataset(dstosplit->getFileName(), stroutgroup, stroutdataset, false);
            
            BigDataStatMeth::hdf5Dataset* dsOut;
            
            try {
                
                dsOut = new BigDataStatMeth::hdf5Dataset(dstosplit->getFileName(), newDatasetName, false);
                dsOut->createDataset( inrows, incols, "real"); 
                dsOut->writeDataset(vdts.data());
                
            } catch( H5::FileIException& error ) {
                delete dsOut;
                ::Rf_error( "c++ exception (File IException )" );
                return void();
            } catch( H5::DataSetIException& error ) { // catch failure caused by the dstosplit operations
                delete dsOut;
                ::Rf_error( "c++ exception (dstosplit IException )" );
                return void();
            } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
                delete dsOut;
                ::Rf_error( "c++ exception (DataSpace IException )" );
                return void();
            } 
            
            delete dsOut;
            
        }
        
    } catch( H5::FileIException& error ) {
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