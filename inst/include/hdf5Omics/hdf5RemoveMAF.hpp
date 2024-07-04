#ifndef BIGDATASTATMETH_OMICS_RM_MAF_HPP
#define BIGDATASTATMETH_OMICS_RM_MAF_HPP

#include "hdf5Utilities/hdf5Utilities.hpp"
#include "hdf5Omics/hdf5OmicsUtils.hpp"



namespace BigDataStatMeth {

    // Removes row or column with high missing data percentage
    int Rcpp_Remove_MAF_hdf5( BigDataStatMeth::hdf5Dataset* dsIn, BigDataStatMeth::hdf5Dataset* dsOut, bool bycols, double pcent, int blocksize)
    {
        
        int itotrem = 0;
        
        try{
        
            bool bcreated = false;    
            int ilimit;
            
            std::vector<hsize_t> offset = {0,0},
                                 count = {0,0},
                                 stride = {1,1},
                                 block = {1,1},
                                 newoffset = {0,0};
            
            // Real data set dimension
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
            
            for( int i=0; i <= (ilimit/blocksize); i++) 
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
                // Eigen::MatrixXd data = GetCurrentBlock_hdf5(file, dataset, offset[0], offset[1], count[0], count[1]);
                std::vector<double> vdCurDataset( count[0] * count[1] ); 
                dsIn->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdCurDataset.data() );
                Eigen::MatrixXd data = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdCurDataset.data(), count[0], count[1] );
                
                if(bycols == true) // We have to do it by rows
                {
                    int readedrows = data.rows();
                    
                    for( int row = readedrows-1 ; row>=0; row--)
                    {
                        if( calc_freq(Rcpp::wrap(data.row(row))) <= pcent ) {
                            removeRow(data, row);
                            iblockrem = iblockrem + 1;
                        } 
                    }
                    
                } else {
                    
                    int readedcols = data.cols();
                    
                    for( int col = readedcols-1 ; col>=0; col--)
                    { 
                        if( calc_freq(Rcpp::wrap(data.col(col))) <= pcent ) {
                            removeColumn(data, col);
                            iblockrem = iblockrem + 1;
                        } 
                        
                    }
                }
                
                int extendcols = data.cols();
                int extendrows = data.rows();
                
                if( extendrows>0 && extendcols>0)
                {
                    if(bcreated == false) {
                        // create_HDF5_unlimited_matrix_dataset_ptr(file, stroutdata, extendrows, extendcols, "numeric");
                        // unlimDataset = new DataSet(file->openDataSet(stroutdata));
                        // bcreated = true;
                        
                        dsOut->createUnlimitedDataset(extendrows, extendcols, "real");
                        dsOut->openDataset();
                        bcreated = true;
                        
                    } else {
                        if(bycols == true){
                            // extend_HDF5_matrix_subset_ptr(file, unlimDataset, extendrows, 0);
                            dsOut->extendUnlimitedDataset(extendrows, 0);
                        }else{
                            // extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, extendcols);
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
                    
                    // IntegerVector countblock = IntegerVector::create(extendrows, extendcols);
                    // write_HDF5_matrix_subset_v2(file, unlimDataset, newoffset, countblock, stride, block, wrap(data) );
                    
                    // if(bycols == true)
                    //     newoffset[0] =  newoffset[0] + extendrows;
                    // else
                    //     newoffset[1] =  newoffset[1] + extendcols;
                    
                    itotrem = itotrem - iblockrem;
                }
            }
            
            // if (bcreated == true) {
            //     unlimDataset->close();
            // }
            
        } catch( H5::FileIException& error) { // catch failure caused by the H5File operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (File IException)" );
            return -1;
        } catch( H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (DataSet IException)" );
            return -1;
        } catch( H5::GroupIException& error) { // catch failure caused by the Group operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (Group IException)" );
            return -1;
        } catch( H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (DataSpace IException)" );
            return -1;
        } catch( H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (Data TypeIException)" );
            return -1;
        } catch(std::exception &ex) {
            Rcpp::Rcout<< ex.what();
            return -1;
        }
        
        
        return(itotrem);
        
        
        
        
        
        // int itotrem = 0;
        // 
        // 
        // 
        // try{
        //     
        //     bool bcreated = false;    
        //     int ilimit;
        //     
        //     std::vector<hsize_t> offset = {0,0},
        //         count = {0,0},
        //         stride = {1,1},
        //         block = {1,1},
        //         newoffset = {0,0};
        //     
        //     // Real data set dimension
        //     IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
        //     
        //     // id bycols == true : read all rows by group of columns ; else : all columns by group of rows
        //     if (bycols == true) {
        //         ilimit = dims_out[0];
        //         count[1] = dims_out[1];
        //         offset[1] = 0;
        //     } else {
        //         ilimit = dims_out[1];
        //         count[0] = dims_out[0];
        //         offset[0] = 0;
        //     };
        //     
        //     for( int i=0; i<=(ilimit/blocksize); i++) 
        //     {
        //         int iread;
        //         int iblockrem = 0;
        //         
        //         if( (i+1)*blocksize < ilimit) iread = blocksize;
        //         else iread = ilimit - (i*blocksize);
        //         
        //         if(bycols == true) {
        //             count[0] = iread; 
        //             offset[0] = i*blocksize;
        //         } else {
        //             count[1] = iread; 
        //             offset[1] = i*blocksize;
        //         }
        //         
        //         // read block
        //         Eigen::MatrixXd data = GetCurrentBlock_hdf5(file, dataset, offset[0], offset[1], count[0], count[1]);
        //         
        //         if(bycols == true) // We have to do it by rows
        //         {
        //             int readedrows = data.rows();
        //             
        //             for( int row = readedrows-1 ; row>=0; row--)
        //             {
        //                 if( calc_freq(wrap(data.row(row))) <= pcent ) {
        //                     removeRow(data, row);
        //                     iblockrem = iblockrem + 1;
        //                 } 
        //             }
        //             
        //         } else {
        //             
        //             int readedcols = data.cols();
        //             
        //             for( int col = readedcols-1 ; col>=0; col--)
        //             { 
        //                 if( calc_freq(wrap(data.col(col))) <= pcent ) {
        //                     removeColumn(data, col);
        //                     iblockrem = iblockrem + 1;
        //                 } 
        //                 
        //             }
        //         }
        //         
        //         
        //         int extendcols = data.cols();
        //         int extendrows = data.rows();
        //         
        //         if( extendrows>0 && extendcols>0)
        //         {
        //             
        //             if(bcreated == false) {
        //                 create_HDF5_unlimited_matrix_dataset_ptr(file, stroutdata, extendrows, extendcols, "numeric");
        //                 unlimDataset = new DataSet(file->openDataSet(stroutdata));
        //                 bcreated = true;
        //             }else {
        //                 if(bycols == true){
        //                     extend_HDF5_matrix_subset_ptr(file, unlimDataset, extendrows, 0);
        //                 }else{
        //                     extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, extendcols);
        //                 }
        //             }
        //             
        //             IntegerVector countblock = IntegerVector::create(extendrows, extendcols);
        //             write_HDF5_matrix_subset_v2(file, unlimDataset, newoffset, countblock, stride, block, wrap(data) );
        //             
        //             if(bycols == true)
        //                 newoffset[0] =  newoffset[0] + extendrows;
        //             else
        //                 newoffset[1] =  newoffset[1] + extendcols;
        //         }
        //         
        //         
        //         itotrem = itotrem - iblockrem;
        //         
        //     }
        //     
        //     if (bcreated == true) {
        //         unlimDataset->close();
        //     }
        //     
        // } catch(FileIException& error) { // catch failure caused by the H5File operations
        //     unlimDataset->close();
        //     ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (File IException)" );
        //     return -1;
        // } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
        //     unlimDataset->close();
        //     ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (DataSet IException)" );
        //     return -1;
        // } catch(GroupIException& error) { // catch failure caused by the Group operations
        //     unlimDataset->close();
        //     ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (Group IException)" );
        //     return -1;
        // } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
        //     unlimDataset->close();
        //     ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (DataSpace IException)" );
        //     return -1;
        // } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
        //     unlimDataset->close();
        //     ::Rf_error( "c++ exception Rcpp_Remove_MAF_hdf5 (Data TypeIException)" );
        //     return -1;
        // }
        // 
        // 
        // return(itotrem);
    }




}

#endif // BIGDATASTATMETH_OMICS_RM_MAF_HPP