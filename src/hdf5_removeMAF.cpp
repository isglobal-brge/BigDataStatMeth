#include <BigDataStatMeth.hpp>
#include "hdf5Omics/hdf5RemoveMAF.hpp"





//' Remove SNPs in hdf5 omic dataset with low data
//'
//' Remove SNPs in hdf5 omic dataset with low data
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param dataset, character array indicating the input dataset to be imputed
//' @param outgroup, character array indicating group where the data set will be saved after remove data with if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
//' @param maf, by default maf = 0.05. Numeric indicating the percentage to be considered to remove SNPs, SNPS with higest MAF will be removed from data
//' @param bycols, boolean by default = true, if true, indicates that SNPs are in cols, if SNPincols = false indicates that SNPs are in rows.
//' @param blocksize, integer, block size dataset to read/write and calculate MAF, by default this operations is made in with 100 rows if byrows = true or 100 cols if byrows = false.
//' @param overwrite, optional boolean if true, previous results in same location 
//' inside hdf5 will be overwritten, by default overwrite = false, data was not overwritten.
//' @return Original hdf5 data file with imputed data
//' @export
// [[Rcpp::export]]
void bdRemoveMAF_hdf5( std::string filename, std::string group, std::string dataset, 
                       std::string outgroup, std::string outdataset, 
                       Rcpp::Nullable<double> maf, Rcpp::Nullable<bool> bycols, 
                       Rcpp::Nullable<int> blocksize, Rcpp::Nullable<bool> overwrite = R_NilValue )
{
    
    // H5File* file = nullptr;
    int iremoved = 0;
    
    try
    {
        
        BigDataStatMeth::hdf5Dataset* dsIn;
        BigDataStatMeth::hdf5DatasetInternal* dsOut;
        
        bool bcols, bforce;
        double dpcent;
        int iblocksize = 100;
        
        std::string stroutdata = outgroup +"/" + outdataset;
        std::string strdataset = group +"/" + dataset;
        
        if(bycols.isNull()){  
            bcols = false ;
        }else{    
            bcols = Rcpp::as<bool>(bycols);
        }
        
        if(maf.isNull()){ dpcent = 0.05 ; } 
        else { dpcent = Rcpp::as<double>(maf); }
        
        if(!blocksize.isNull()){  
            iblocksize = Rcpp::as<int>(blocksize);
        }
        
        if(overwrite.isNull()) { bforce = false ; }
        else { bforce = Rcpp::as<bool>(overwrite); }
        
        
        // if(!ResFileExist_filestream(filename)){
        //     throw std::range_error("File not exits, create file before access to dataset");
        // }
        // 
        // 
        // file = new H5File( filename, H5F_ACC_RDWR );
        // 
        // if(exists_HDF5_element_ptr(file, strdataset)) 
        // {
        //     
        //     DataSet* pdataset = nullptr;
        //     
        //     pdataset = new DataSet(file->openDataSet(strdataset));
            
            
            if( strdataset.compare(stroutdata)!= 0)
            {
                dsIn = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
                dsIn->openDataset();
                
                dsOut = new BigDataStatMeth::hdf5DatasetInternal(filename, outgroup, outdataset, bforce);
                
            } else {
                throw std::range_error("Input and output dataset must be different");  
                return void();
            }
            
            
            // if( strdataset.compare(stroutdata)!= 0)
            // {
            //     
            //     // If output is different from imput --> Remve possible existing dataset and create new
            //     if(exists_HDF5_element_ptr(file, stroutdata))
            //         remove_HDF5_element_ptr(file, stroutdata);
            //     
            //     // Create group if not exists
            //     if(!exists_HDF5_element_ptr(file, outgroup))
            //         file->createGroup(outgroup);
            //     
            // } else {
            //     throw std::range_error("Input and output dataset must be different");  
            // }
            
            iremoved = Rcpp_Remove_MAF_hdf5( dsIn, dsOut, bcols, dpcent, iblocksize);
            
            Rcpp::Function warning("warning");
            if (!bcols )
                warning( std::to_string(iremoved) + " Rows have been removed");
            else
                warning( std::to_string(iremoved) + " Columns have been removed");
            
            // pdataset->close();
            
        // } else{
        //     //.commented 20201120 - warning check().// pdataset->close();
        //     file->close();
        //     throw std::range_error("Dataset does not exits");  
        // }
        
        delete dsIn;
        delete dsOut;
        
    } catch( H5::FileIException& error ){ // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception bdRemoveMAF_hdf5 (File IException)" );
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception bdRemoveMAF_hdf5 (DataSet IException)" );
        return void();
    } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception bdRemoveMAF_hdf5 (DataSpace IException)" );
        return void();
    } catch( H5::DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception bdRemoveMAF_hdf5 (DataType IException)" );
        return void();
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return void();
    }
    
    // file->close();
    return void();
    
    
    // H5File* file = nullptr;
    // int iremoved = 0;
    // int iblocksize = 100;
    // 
    // try
    // {
    //     bool bcols;
    //     double dpcent;
    //     
    //     std::string stroutdata = outgroup +"/" + outdataset;
    //     std::string strdataset = group +"/" + dataset;
    //     
    //     if(bycols.isNull()){  
    //         bcols = false ;
    //     }else{    
    //         bcols = Rcpp::as<bool>(bycols);
    //     }
    //     
    //     if(maf.isNull()){  
    //         dpcent = 0.05 ;
    //     }else{    
    //         dpcent = Rcpp::as<double>(maf);
    //     }
    //     
    //     if(!blocksize.isNull()){  
    //         iblocksize = Rcpp::as<int>(blocksize);
    //     }
    //     
    //     
    //     
    //     if(!ResFileExist_filestream(filename)){
    //         throw std::range_error("File not exits, create file before access to dataset");
    //     }
    //     
    //     
    //     file = new H5File( filename, H5F_ACC_RDWR );
    //     
    //     if(exists_HDF5_element_ptr(file, strdataset)) 
    //     {
    //         
    //         DataSet* pdataset = nullptr;
    //         
    //         pdataset = new DataSet(file->openDataSet(strdataset));
    //         
    //         if( strdataset.compare(stroutdata)!= 0)
    //         {
    //             
    //             // If output is different from imput --> Remve possible existing dataset and create new
    //             if(exists_HDF5_element_ptr(file, stroutdata))
    //                 remove_HDF5_element_ptr(file, stroutdata);
    //             
    //             // Create group if not exists
    //             if(!exists_HDF5_element_ptr(file, outgroup))
    //                 file->createGroup(outgroup);
    //             
    //         } else {
    //             throw std::range_error("Input and output dataset must be different");  
    //         }
    //         
    //         iremoved = Remove_MAF_HDF5( file, pdataset, bcols, stroutdata, dpcent, iblocksize);
    //         
    //         Function warning("warning");
    //         if (!bycols )
    //             warning( std::to_string(iremoved) + " Rows have been removed");
    //         else
    //             warning( std::to_string(iremoved) + " Columns have been removed");
    //         
    //         pdataset->close();
    //         
    //     } else{
    //         //.commented 20201120 - warning check().// pdataset->close();
    //         file->close();
    //         throw std::range_error("Dataset does not exits");  
    //     }
    //     
    //     
    // }catch( FileIException& error ){ // catch failure caused by the H5File operations
    //     file->close();
    //     ::Rf_error( "c++ exception bdremove_maf_hdf5 (File IException)" );
    //     return(wrap(-1));
    // } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    //     file->close();
    //     ::Rf_error( "c++ exception bdremove_maf_hdf5 (DataSet IException)" );
    //     return(wrap(-1));
    // } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
    //     file->close();
    //     ::Rf_error( "c++ exception bdremove_maf_hdf5 (DataSpace IException)" );
    //     return(wrap(-1));
    // } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
    //     file->close();
    //     ::Rf_error( "c++ exception bdremove_maf_hdf5 (DataType IException)" );
    //     return(wrap(-1));
    // }catch(std::exception &ex) {
    //     file->close();
    //     Rcpp::Rcout<< ex.what();
    //     return(wrap(-1));
    // }
    // 
    // file->close();
    // return(wrap(iremoved));
    
}

