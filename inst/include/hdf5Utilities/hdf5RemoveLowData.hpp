#ifndef BIGDATASTATMETH_UTIL_QC_BASICS_HPP
#define BIGDATASTATMETH_UTIL_QC_BASICS_HPP

#include "hdf5Utilities.hpp"

namespace BigDataStatMeth {


    // Removes row or column with high missing data percentage
    extern inline int Rcpp_Remove_Low_Data_hdf5( BigDataStatMeth::hdf5Dataset* dsIn, BigDataStatMeth::hdf5Dataset* dsOut, bool bycols, double pcent)
    {
        
        int itotrem = 0;
        
        try{
        
            int ilimit,
                blocksize = 1000;
            
            bool bcreated = false;
            
            std::vector<hsize_t> offset = {0,0},
                                 count = {0,0},
                                 stride = {1,1},
                                 block = {1,1},
                                 newoffset = {0,0};
            
            hsize_t* dims_out = dsIn->dim();
            
            // id bycols == true : read all rows by group of columns ; else : all columns by group of rows
            if (bycols == true) {
                ilimit = dims_out[0];
                count[1] = dims_out[1];
                offset[1] = 0;
            } else {
                ilimit = dims_out[1];
                count[0] = dims_out[0];
                offset[0] = 0;
            };
            
            for( int i=0; i<=(ilimit/blocksize); i++) 
            {
                int iread;
                int iblockrem = 0;
                
                if( (i+1)*blocksize < ilimit) iread = blocksize;
                else iread = ilimit - (i*blocksize);
                
                if(bycols == true) {
                    count[0] = iread; 
                    offset[0] = i*blocksize;
                } else {
                    count[1] = iread; 
                    offset[1] = i*blocksize;
                }
                
                // read block
                std::vector<double> vdCurDataset( count[0] * count[1] ); 
                dsIn->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdCurDataset.data() );
                Eigen::MatrixXd data = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdCurDataset.data(), count[0], count[1] );
                
                
                if(bycols == true) // We have to do it by rows
                {
                    int readedrows = data.rows();
                    
                    for( int row = readedrows-1 ; row>=0; row--)  // COMPLETE EXECUTION
                    {
                        if((data.row(row).array() == 3).count()/(double)count[1]>= pcent )
                        {
                            removeRow(data, row);
                            iblockrem = iblockrem + 1;
                        } 
                    }
                    
                } else {
                    
                    int readedcols = data.cols();
                    
                    for( int col = readedcols-1 ; col>=0; col--)  // COMPLETE EXECUTION
                    { 
                        if((data.col(col).array() == 3).count()/(double)count[0]>=pcent )
                        {
                            removeColumn(data, col);
                            iblockrem = iblockrem + 1;
                        } 
                    }
                }
                
                int extendcols = data.cols();
                int extendrows = data.rows();
                
                if( bcreated == false) {
                    
                    if( extendcols > 0 && extendrows > 0){
                        dsOut->createUnlimitedDataset(extendrows, extendcols, "real");
                        dsOut->openDataset();
                        bcreated = true;
                    } 
                    
                } else {
                    if(bycols == true){
                        dsOut->extendUnlimitedDataset(extendrows, 0);
                    }else{
                        dsOut->extendUnlimitedDataset( 0, extendcols);
                    }
                }

                if( bcreated == true) {
                    std::vector<hsize_t> countblock = {(unsigned long long)extendrows, (unsigned long long)extendcols};
                    dsOut->writeDatasetBlock( Rcpp::wrap(data), newoffset, countblock, stride, block, false);
                    
                    if(bycols == true)
                        newoffset[0] =  newoffset[0] + extendrows;
                    else
                        newoffset[1] =  newoffset[1] + extendcols;
                }
                
                itotrem = itotrem - iblockrem;
                
            }
            
            if( bcreated == false ) {
                Rcpp::warning("All data removed - please adjust pcent parameter or review data");
            }
            
        } catch( H5::FileIException& error) { // catch failure caused by the H5File operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_Low_Data_hdf5 (File IException)" );
            return -1;
        } catch( H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_Low_Data_hdf5 (DataSet IException)" );
            return -1;
        } catch( H5::GroupIException& error) { // catch failure caused by the Group operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_Low_Data_hdf5 (Group IException)" );
            return -1;
        } catch( H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_Low_Data_hdf5 (DataSpace IException)" );
            return -1;
        } catch( H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_Low_Data_hdf5 (Data TypeIException)" );
            return -1;
        }
        
        return(itotrem);
    }



}

#endif // BIGDATASTATMETH_UTIL_QC_BASICS_HPP