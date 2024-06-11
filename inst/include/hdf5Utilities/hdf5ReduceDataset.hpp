#ifndef BIGDATASTATMETH_UTIL_REDUCE_DATASETS_HPP
#define BIGDATASTATMETH_UTIL_REDUCE_DATASETS_HPP



namespace BigDataStatMeth {


    // NOTE: 
    //      Set 
    //          - binternal = true for internal calls, for example if we are
    //          performing a reduction as a intermediate step using C++ interface,
    //          - binternal = false only for R interface (get transposed results)
    void RcppReduce_dataset_hdf5 ( std::string filename, 
                                   std::string stringroup, 
                                   std::string stroutgroup, 
                                   std::string stroutdataset, 
                                   std::string strreducefunction, 
                                   bool boverwrite,
                                   bool bremove,
                                   bool binternal)
    {
        
        
        // BigDataStatMeth::hdf5File* file = new BigDataStatMeth::hdf5File(filename, false);
        
        hsize_t* dims_out;
        std::vector<hsize_t> stride = {1, 1},
                             block = {1, 1},
                             offset = {0, 0};
        
        // DataSet* ploaddataset = nullptr;
        Eigen::MatrixXd fullReduced;
        Eigen::MatrixXd newRead;
        int ndatasets;
        
        try {
            
            BigDataStatMeth::hdf5File* objFile = new BigDataStatMeth::hdf5File(filename, false);
            objFile->openFile("r");
            
            // Get dataset names without prefix, all datasets inside the group
            Rcpp::StringVector joindata =  objFile->getDatasetNames(stringroup, "");
            
            delete objFile; // Close file 
            
            ndatasets = joindata.size();
            
            for ( int i=0; i< ndatasets; i++)
            {
                BigDataStatMeth::hdf5Dataset* dsIn = new BigDataStatMeth::hdf5Dataset(filename, stringroup + "/" + joindata[i], false);
                dsIn->openDataset();
                
                dims_out =   dsIn->dim();
                
                std::vector<double> vdIn( dims_out[0] * dims_out[1] ); 
                dsIn->readDatasetBlock( {offset[0], offset[1]}, {dims_out[0], dims_out[1]}, stride, block, vdIn.data() );
                
                if( i == 0 ) {
                    // fullReduced = GetCurrentBlock_hdf5_Original( file, ploaddataset, 0, 0, dims_out[0],dims_out[1]);
                    if(binternal == true)
                        fullReduced = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdIn.data(), dims_out[0], dims_out[1] );
                    else
                        fullReduced = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> (vdIn.data(), dims_out[0], dims_out[1] );
                    
                } else {
                    
                    if(binternal == true)
                        newRead = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>  (vdIn.data(), dims_out[0], dims_out[1] );
                    else
                        newRead = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>  (vdIn.data(), dims_out[0], dims_out[1] );
                    
                    if( newRead.rows() != fullReduced.rows()){
                        
                        int difference = std::abs(fullReduced.rows() - newRead.rows());
                        if( newRead.rows() > fullReduced.rows()) {
                            newRead.resize( newRead.rows() + difference, Eigen::NoChange);    
                        } else {
                            fullReduced.resize( fullReduced.rows() + difference, Eigen::NoChange);    
                        }
                    }
                    
                    if( newRead.cols() != fullReduced.cols()){
                        
                        int difference = std::abs(fullReduced.cols() - newRead.cols());
                        if( newRead.cols() > fullReduced.cols()){
                            newRead.resize( Eigen::NoChange, newRead.cols() + difference );
                        } else {
                            fullReduced.resize( Eigen::NoChange, fullReduced.cols() + difference );
                        }
                    }
                    
                    // Reduce matrix
                    if( strreducefunction.compare("+")==0) {
                        fullReduced = fullReduced + newRead;
                    } else if (strreducefunction.compare("-")==0) {
                        fullReduced = fullReduced - newRead;
                    } 
                }
                
                if( bremove == true){
                    dsIn->remove();
                }
                
                delete dsIn;
                
            }
            
            BigDataStatMeth::hdf5Dataset* dsOut = new BigDataStatMeth::hdf5Dataset(filename, stroutgroup, stroutdataset, boverwrite);
            
            if(binternal == true) {
                dsOut->createDataset( fullReduced.rows() , fullReduced.cols(), "real");
                dsOut->writeDataset(fullReduced.data());    
            } else {
                dsOut->createDataset( fullReduced.cols() , fullReduced.rows(), "real");
                dsOut->writeDataset(fullReduced.transpose().data());
            }
            
            delete dsOut;
            // delete file;
            
        }catch( H5::FileIException& error ) {
            ::Rf_error( "c++ exception RcppReduce_dataset_hdf5 (File IException )" );
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the dstosplit operations
            ::Rf_error( "c++ exception RcppReduce_dataset_hdf5 (dstosplit IException )" );
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception RcppReduce_dataset_hdf5 (DataSpace IException )" );
            return void();
        } 
        
        return void();
        
        
        
        
        
        // DataSet* ploaddataset = nullptr;
        // Eigen::MatrixXd fullReduced;
        // IntegerVector dims_out;
        // int ndatasets;
        // std::string strdatasetout = stroutgroup + "/" + stroutdataset;
        // 
        // try {
        //     
        //     // Get dataset names without prefix, all datasets inside the group
        //     StringVector joindata =  get_dataset_names_from_group(file, strgroup, "");
        //     
        //     ndatasets = joindata.size();
        //     
        //     for ( int i=0; i< ndatasets; i++)
        //     {
        //         ploaddataset = new DataSet(file->openDataSet(strgroup + "/" + joindata[i]));
        //         dims_out = get_HDF5_dataset_size(*ploaddataset);
        //         
        //         if( i == 0 ) {
        //             fullReduced = GetCurrentBlock_hdf5_Original( file, ploaddataset, 0, 0, dims_out[0],dims_out[1]);
        //             
        //         } else {
        //             
        //             // If readed block is smaller than full matrix adds rows or columns
        //             Eigen::MatrixXd newRead = GetCurrentBlock_hdf5_Original( file, ploaddataset, 0, 0, dims_out[0],dims_out[1]);
        //             if( newRead.rows() != fullReduced.rows()){
        //                 int difference = fullReduced.rows() - newRead.rows();
        //                 newRead.conservativeResize(newRead.rows() + difference, newRead.cols()); 
        //             }
        //             
        //             if( newRead.cols() != fullReduced.cols()){
        //                 int difference = fullReduced.cols() - newRead.cols();
        //                 newRead.conservativeResize(newRead.rows(), newRead.cols() + di
        //                 fference); 
        //             }
        //             
        //             // Reduce matrix
        //             if( strreducefunction.compare("+")==0) {
        //                 fullReduced = fullReduced + newRead;
        //             } else if (strreducefunction.compare("-")==0) {
        //                 fullReduced = fullReduced - newRead;
        //             } 
        //         }
        //         
        //         ploaddataset->close();
        //         
        //         if( bremove == true){
        //             // Remove used dataset
        //             remove_HDF5_element_ptr(file, strgroup + "/" + joindata[i]);
        //         }
        //         
        //     }
        //     
        //     // Transform to rowmajor
        //     Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mapBlock(fullReduced.data(), 
        //                                                                                               fullReduced.cols() , fullReduced.rows());
        //     
        //     write_HDF5_matrix_from_R_ptr(file, strdatasetout, Rcpp::wrap(mapBlock.transpose()), false);
        //     
        // }catch( FileIException& error ) {
        //     ploaddataset->close();
        //     file->close();
        //     ::Rf_error( "c++ exception (File IException )" );
        //     return -1;
        // } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        //     ploaddataset->close();
        //     file->close();
        //     ::Rf_error( "c++ exception (DataSet IException )" );
        //     return -1;
        // } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        //     ploaddataset->close();
        //     file->close();
        //     ::Rf_error( "c++ exception (DataSpace IException )" );
        //     return -1;
        // } 
        // 
        // return(0);
        
        
    }




}

#endif // BIGDATASTATMETH_UTIL_REDUCE_DATASETS_HPP