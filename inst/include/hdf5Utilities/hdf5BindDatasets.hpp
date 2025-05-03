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
    
    extern inline void RcppBind_datasets_hdf5( std::string filename, std::string group, 
                                 Rcpp::StringVector datasets, 
                                 BigDataStatMeth::hdf5Dataset* dsOut,  
                                 int func, bool binternal )
    {
        
        
        BigDataStatMeth::hdf5Dataset* dsIn = nullptr;
        
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
                
                dsIn = new BigDataStatMeth::hdf5Dataset(filename, strdataset, false);
                dsIn->openDataset();
                
                // Real data set dimension
                dims_out =   dsIn->dim();
                
                // Get block from complete matrix
                std::vector<double> vdIn( dims_out[0] * dims_out[1] ); 
                dsIn->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, stride, block, vdIn.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> original (vdIn.data(), dims_out[0], dims_out[1] );    
                
                delete dsIn; dsIn = nullptr;
                
                
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
                    delete dsIn; dsIn = nullptr;
                    Rcpp::Rcerr<<"Group not exists, create the input datasets before proceed";
                    return void();
                }
            }
            
        } catch( H5::FileIException& error ) {
            checkClose_file(dsIn, dsOut);
            Rcpp::Rcerr<<"\nc++ exception RcppBind_datasets_hdf5 (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the dstosplit operations
            checkClose_file(dsIn, dsOut);
            Rcpp::Rcerr<<"\nc++ exception RcppBind_datasets_hdf5 (dstosplit IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            checkClose_file(dsIn, dsOut);
            Rcpp::Rcerr<<"\nc++ exception RcppBind_datasets_hdf5 (DataSpace IException)\n";
            return void();
        } 
        
        return void();
        
    }




    extern inline void RcppBind_datasets_hdf5( std::string filename, std::string group, Rcpp::StringVector datasets, 
                               std::string outgroup, std::string outdataset, std::string func, 
                               bool binternal, Rcpp::Nullable<bool> overwrite = false )
    {
        
        BigDataStatMeth::hdf5Dataset* dsOut = nullptr;
        
        try
        {
            
            Rcpp::NumericVector oper = {0, 1, 2};
            oper.names() = Rcpp::CharacterVector({ "bindCols", "bindRows", "bindRowsbyIndex"});
            
            bool boverwrite;
            
            if( overwrite.isNull()) { boverwrite = false; } 
            else {   boverwrite = Rcpp::as<bool>(overwrite); }
            
            if (func.compare("bindCols") != 0 && func.compare("bindRows") != 0  && func.compare("bindRowsbyIndex") != 0 ) {
                throw std::range_error( "Function to apply must be \"bindRows\", \"bindCols\" or \"bindRowsbyIndex\" other values are not allowed" );
                return void();
            }
            
            int bindFunction = oper.findName( func );
            
            dsOut = new BigDataStatMeth::hdf5Dataset(filename, outgroup, outdataset, boverwrite);
            
            RcppBind_datasets_hdf5( filename, group, datasets, dsOut, bindFunction, binternal);
            
            delete dsOut; dsOut = nullptr;
            
        } catch( H5::FileIException& error ) { 
            checkClose_file(dsOut);
            Rcpp::Rcerr<<"\nc++ exception RcppBind_datasets_hdf5_ (File IException)";
            return void();
        } catch( H5::GroupIException & error ) { 
            checkClose_file(dsOut);
            Rcpp::Rcerr <<"\nc++ exception RcppBind_datasets_hdf5_ (Group IException)";
            return void();
        } catch( H5::DataSetIException& error ) { 
            checkClose_file(dsOut);
            Rcpp::Rcerr <<"\nc++ exception RcppBind_datasets_hdf5_ (DataSet IException)";
            return void();
        } catch(std::exception& ex) {
            checkClose_file(dsOut);
            Rcpp::Rcerr <<"\nc++ exception RcppBind_datasets_hdf5_" << ex.what();
            return void();
        }
        
        return void();
        
    }


}

#endif // BIGDATASTATMETH_UTIL_BIND_DATASETS_HPP