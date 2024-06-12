#ifndef BIGDATASTATMETH_UTIL_BIND_DATASETS_HPP
#define BIGDATASTATMETH_UTIL_BIND_DATASETS_HPP



namespace BigDataStatMeth {


    // 
    // Possible functions: 
    //  if oper = 
    //      0: Bind by Cols
    //      1: Bind by Rows
    //      2: Bind by Rows taking in to account an index
    //      
    //  IMPORTANT: Take in to account that in R data is transposed when stored 
    //          in hdf5 data files
    //          
    
    void RcppBind_datasets_hdf5( std::string filename, std::string group, 
                                 Rcpp::StringVector datasets, 
                                 BigDataStatMeth::hdf5Dataset* dsOut,  
                                 int func, bool binternal )
    {
        
        try {
            
            
            hsize_t* dims_out;
            std::vector<hsize_t> stride = {1, 1},
                                 block = {1, 1},
                                 offset = {0, 0},
                                 count = {0, 0};
            
            // Seek all datasets to perform calculus
            for( int i=0; i < datasets.size(); i++ ) 
            {
                
                std::string strdataset = group +"/" + datasets(i);
                
                BigDataStatMeth::hdf5Dataset* dsIn = new BigDataStatMeth::hdf5Dataset(filename, strdataset, false);
                dsIn->openDataset();
                
                // Real data set dimension
                dims_out =   dsIn->dim();
                
                // Get block from complete matrix
                std::vector<double> vdIn( dims_out[0] * dims_out[1] ); 
                dsIn->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, stride, block, vdIn.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> original (vdIn.data(), dims_out[0], dims_out[1] );    
                
                delete dsIn;
                
                
                if( func == 0 || func == 1) {
                    
                    if( func == 0 ) { // byCols
                        
                        // Test if dimmensions are correct
                        if( original.cols() != count[0] && i!=0) {
                            
                            if( count[0] > original.cols()) {
                                // Append needed cols to merge by cols
                                int iappend = count[0] - original.cols();
                                //..// original.conservativeResize(original.rows(), original.cols() + iappend);
                                original.resize(original.rows(), original.cols() + iappend);
                                
                            } else {
                                std::string strmessage = "Can't bind current dataset, number of rows differ on size";
                                Rcpp::message(Rcpp::wrap(strmessage));
                                return void();
                            }
                        }
                        
                        offset[0] = offset[0] + count[1];
                        
                    } else { // byRows
                        
                        // Check if dimmensions are correct
                        if( original.rows() != count[0]  && i!=0) {
                            
                            if( count[0] > original.rows()) {
                                // Append needed rows to merge by rows
                                int iappend = count[0] - original.rows();
                                original.resize(original.rows() + iappend, original.cols());
                                
                            } else {
                                std::string strmessage = "Can't bind current dataset, number of columns differ on size";
                                Rcpp::message(Rcpp::wrap(strmessage));
                                return void();
                            }
                            
                        }
                        offset[1] = offset[1] + count[0];
                    }
                    
                    count[0] = original.cols();
                    count[1] = original.rows();
                    
                    if(i == 0) 
                        dsOut->createUnlimitedDataset(count[0], count[1], "real");
                    dsOut->openDataset();
                    
                    if( func == 1 && i!=0) {
                        dsOut->extendUnlimitedDataset(count[0], 0 );
                    } else if ( func == 0 && i!=0) {
                        dsOut->extendUnlimitedDataset( 0, count[1] );
                    }
                    
                    dsOut->writeDatasetBlock( Rcpp::wrap(original), offset, count, stride, block, true);
                    
                } else {
                    
                    Rcpp::Rcout<<"\nDebug else";
                    
                    delete dsIn;
                    // pdataset->close();
                    // file->close();
                    Rcpp::Rcout<<"Group not exists, create the input datasets before proceed";
                    return void();
                }
            }
            
            
        } catch( H5::FileIException& error ) {
            ::Rf_error( "c++ exception RcppBind_datasets_hdf5 (File IException )" );
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the dstosplit operations
            ::Rf_error( "c++ exception RcppBind_datasets_hdf5 (dstosplit IException )" );
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception RcppBind_datasets_hdf5 (DataSpace IException )" );
            return void();
        } 
        
        
        return void();
        
        
    }


}

#endif // BIGDATASTATMETH_UTIL_BIND_DATASETS_HPP