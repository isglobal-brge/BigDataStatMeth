#include <BigDataStatMeth.hpp>
#include "hdf5Utilities/hdf5ImputeData.hpp"


//' Impute SNPs in hdf5 omic dataset 
//'
//' Impute SNPs in hdf5 omic dataset 
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param dataset, character array indicating the input dataset to be imputed
//' @param bycols, boolean by default = true, true indicates that the imputation will be done by columns, otherwise, the imputation will be done by rows
//' @param outgroup, optional character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, optional character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
//' @param overwrite, optional boolean if true, previous results in same location 
//' inside hdf5 will be overwritten, by default overwrite = false, data was not overwritten.
//' @return Original hdf5 data file with imputed data
//' @examples
//' print('see vignette')
//' @export
// [[Rcpp::export]]
void bdImputeSNPs_hdf5(std::string filename, std::string group, std::string dataset, 
                        Rcpp::Nullable<std::string> outgroup = R_NilValue, Rcpp::Nullable<std::string> outdataset = R_NilValue, 
                        Rcpp::Nullable<bool> bycols = true, Rcpp::Nullable<bool> overwrite = R_NilValue )
{
    
    
    // H5File* file = nullptr;
    // DataSet* pdataset = nullptr;
    
    try
    {
        BigDataStatMeth::hdf5Dataset* dsIn;
        BigDataStatMeth::hdf5DatasetInternal* dsOut;
        
        std::string strdataset = group +"/" + dataset;
        std::string stroutgroup, stroutdataset, stroutdata;
        // std::string strdatasetout;
        
        //.commented 20201120 - warning check().// int res;
        bool bcols, bforce;
        
        
        if(bycols.isNull()){  bcols = true ; }
        else{  bcols = Rcpp::as<bool>(bycols); }
        
        if(outgroup.isNull())  stroutgroup = group ;
        else    stroutgroup = Rcpp::as<std::string>(outgroup);
        
        if(outdataset.isNull())  stroutdataset = dataset ;
        else    stroutdataset = Rcpp::as<std::string>(outdataset);
        
        if(overwrite.isNull()) { bforce = false ; }
        else { bforce = Rcpp::as<bool>(overwrite); }
        
        stroutdata = stroutgroup + "/" + stroutdataset;
        
        dsIn = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
        dsIn->openDataset();
        
        dsOut = new BigDataStatMeth::hdf5DatasetInternal(filename, stroutgroup, stroutdataset, bforce);

                
        Rcpp_Impute_snps_hdf5( dsIn, dsOut, bcols, stroutdataset);
        
                
        delete dsIn;
        delete dsOut;
        
    } catch( H5::FileIException& error ){ // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception bdImputeSNPs_hdf5 (File IException)" );
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception bdImputeSNPs_hdf5 (DataSet IException)" );
        return void();
    } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception bdImputeSNPs_hdf5 (DataSpace IException)" );
        return void();
    } catch( H5::DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception bdImputeSNPs_hdf5 (DataType IException)" );
        return void();
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return void();
    }
    
    Rcpp::Rcout<<"SNPs with missing values has been imputed\n";
    return void();
    
    
    
    
    // H5File* file = nullptr;
    // DataSet* pdataset = nullptr;
    // 
    // try
    // {
    //     std::string strdataset = group +"/" + dataset;
    //     std::string stroutgroup, stroutdataset, stroutdata;
    //     std::string strdatasetout;
    //     
    //     //.commented 20201120 - warning check().// int res;
    //     bool bcols;
    //     
    //     
    //     if(bycols.isNull()){  bcols = true ; }
    //     else{  bcols = Rcpp::as<bool>(bycols); }
    //     
    //     if(outgroup.isNull())  stroutgroup = group ;
    //     else    stroutgroup = Rcpp::as<std::string>(outgroup);
    //     
    //     if(outdataset.isNull())  stroutdataset = dataset ;
    //     else    stroutdataset = Rcpp::as<std::string>(outdataset);
    //     
    //     if(overwrite.isNull()) { bforce = false ; }
    //     else { bforce = Rcpp::as<bool>(overwrite); }
    //     
    //     stroutdata = stroutgroup +"/" + stroutdataset;
    //     
    //     
    //     if(!ResFileExist_filestream(filename)){
    //         throw std::range_error("File not exits, create file before access to dataset");
    //     }
    //     
    //     file = new H5File( filename, H5F_ACC_RDWR );
    //     
    //     
    //     if(exists_HDF5_element_ptr(file, strdataset)) 
    //     {
    //         try
    //         {
    //             pdataset = new DataSet(file->openDataSet(strdataset));
    //             
    //             if( strdataset.compare(stroutdata)!= 0)
    //             {
    //                 // If output is different from imput --> Remve possible existing dataset and create new
    //                 if(exists_HDF5_element_ptr(file, stroutdata))
    //                     remove_HDF5_element_ptr(file, stroutdata);
    //                 
    //                 // Create group if not exists
    //                 if(!exists_HDF5_element_ptr(file, stroutgroup))
    //                     file->createGroup(stroutgroup);
    //                 
    //             } else {
    //                 stroutdata = "";
    //             }
    //             
    //             Impute_snp_HDF5( file, pdataset, bcols, stroutdata);
    //             
    //         }catch(FileIException& error) {
    //             pdataset->close(); //.created 20201120 - warning check().//
    //             file->close();
    //         }
    //         
    //         
    //     } else{
    //         //.commented 20201120 - warning check().// pdataset->close();
    //         file->close();
    //         throw std::range_error("Dataset not exits");  
    //     }
    //     
    //     
    // }
    // catch( FileIException& error ) { // catch failure caused by the H5File operations
    //     //.commented 20201120 - warning check().// pdataset->close();
    //     file->close();
    //     ::Rf_error( "c++ exception (File IException)" );
    //     return void();
    // }
    // 
    // pdataset->close();
    // file->close();
    // Rcpp::Rcout<<"SNPs with missing values has been imputed\n";
    // return void();
    
}

