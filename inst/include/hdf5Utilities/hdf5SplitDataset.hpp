#ifndef BIGDATASTATMETH_UTIL_SPLIT_DATASETS_HPP
#define BIGDATASTATMETH_UTIL_SPLIT_DATASETS_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
// #include "memAlgebra/memMultiplication.hpp"
// #include <thread>

namespace BigDataStatMeth {

    
    // Internal call 
    extern inline void RcppSplit_matrix_hdf5 ( BigDataStatMeth::hdf5Dataset* dstosplit, bool bycols, 
                                std::string stroutgroup, std::string stroutdataset, 
                                int blocksize, int irows, int icols )
    {
        
        BigDataStatMeth::hdf5Dataset* dsOut = nullptr;
        
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
                } else  {
                    ii = i * blocksize;
                    if( ii + blocksize > irows) { inrows = irows - ii; }
                }
                
                std::vector<double> vdts( inrows * incols );
                dstosplit->readDatasetBlock( {kk, ii}, {incols, inrows}, stride, block, vdts.data() );
                    
                dsOut = new BigDataStatMeth::hdf5Dataset(dstosplit->getFileName(), newDatasetName, true);
                dsOut->createDataset( inrows, incols, "real"); 
                
                if( dsOut->getDatasetptr() != nullptr ){
                    dsOut->writeDataset(vdts.data());
                }
                
                delete dsOut; dsOut = nullptr;
                
            }
            
        } catch( H5::FileIException& error ) {
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"c++ exception RcppSplit_matrix_hdf5(File IException )";
            return void();
        } catch( H5::DataSetIException& error ) { 
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"c++ exception RcppSplit_matrix_hdf5 (DataSet IException )";
            return void();
        } catch( H5::DataSpaceIException& error ) { 
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"c++ exception RcppSplit_matrix_hdf5 (DataSpace IException )";
            return void();
        } catch(std::exception &ex) {
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"\nC++ exception RcppSplit_matrix_hdf5 : "<< ex.what();
        } catch (...) {
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"\nC++ exception RcppSplit_matrix_hdf5 (unknown reason)";
            return void();
        }
        
        return void();
        
    }
    
    
    // Works with C objects not R (transposed legnths) - to work with R object transposed
    // lengths use RcppSplit_matrix_hdf5
    extern inline void RcppSplit_matrix_hdf5_internal ( BigDataStatMeth::hdf5Dataset* dstosplit,
                                       std::string stroutgroup, std::string stroutdataset,
                                       bool bycols, int nblocks, int iblocksize, 
                                       int irows, int icols )
    {
        
        BigDataStatMeth::hdf5Dataset* dsOut = nullptr;
        
        try {
            
            std::string newDatasetName = "";
            hsize_t inrows = irows, 
                    incols = icols,
                    ii = 0,
                    kk = 0;
            int blocks = nblocks, 
                blocksize = iblocksize;;
            std::vector<hsize_t> stride = {1, 1},
                                 block = {1, 1};
            
            if( nblocks <= 0 && iblocksize <= 0 ){
                checkClose_file(dstosplit);
                Rcpp::Rcerr<<"\n Block size or number of blocks are needed to proceed with matrix split. Please, review parameters";
                return void();
                
            } else if (nblocks > 0 && iblocksize > 0 ) {
                checkClose_file(dstosplit);
                Rcpp::Rcerr<<"\nBlock size and number of blocks are defined, please define only one option, split by number of blocks or by block size";
                return void();
                
            } else if( nblocks > 0 ) {
                double module;
                
                if(bycols == true) {
                    blocksize = icols / blocks;
                    module = icols % blocksize;
                } else {
                    blocksize = irows / blocks;
                    module = irows % blocksize;
                }
                if (module > 0) { blocksize = blocksize + 1; }
            } 
            
            if(blocksize > 0) {
                
                if( bycols == true ) {
                    blocks = (icols + blocksize - 1) / blocksize;
                    incols = blocksize;
                } else {
                    blocks = (irows + blocksize - 1) / blocksize;
                    inrows = blocksize;
                }    
            }
            
            for ( int i=0; i<blocks; i++)
            {
                newDatasetName = stroutgroup + "/" + stroutdataset + "." + std::to_string(i);
                
                if( bycols == true) {
                    kk = i * blocksize;
                    if( kk + blocksize > icols) { incols = icols - kk; }
                } else  {
                    ii = i * blocksize;
                    if( ii + blocksize > irows) { inrows = irows - ii; }
                }
    
                std::vector<double> vdts( inrows * incols );
                dstosplit->readDatasetBlock( {ii, kk}, {inrows, incols}, stride, block, vdts.data() );
                
                dsOut = new BigDataStatMeth::hdf5Dataset(dstosplit->getFileName(), newDatasetName, true);
                dsOut->createDataset(  incols, inrows, "real"); 
                
                if( dsOut->getDatasetptr() != nullptr ){
                    dsOut->writeDataset(vdts.data());
                }
                
                delete dsOut; dsOut = nullptr;
                
            }
            
        } catch( H5::FileIException& error ) {
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"c++ exception RcppSplit_matrix_hdf5_internal(File IException )";
            return void();
        } catch( H5::DataSetIException& error ) { 
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"c++ exception RcppSplit_matrix_hdf5_internal (DataSet IException )";
            return void();
        } catch( H5::DataSpaceIException& error ) { 
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"c++ exception RcppSplit_matrix_hdf5_internal (DataSpace IException )";
            return void();
        } catch(std::exception &ex) {
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"\nC++ exception RcppSplit_matrix_hdf5_internal : "<< ex.what();
        } catch (...) {
            checkClose_file(dstosplit, dsOut);
            Rcpp::Rcerr<<"\nC++ exception RcppSplit_matrix_hdf5_internal (unknown reason)";
            return void();
        }
        
        return void();
        
    }


}

#endif // BIGDATASTATMETH_UTIL_SPLIT_DATASETS_HPP