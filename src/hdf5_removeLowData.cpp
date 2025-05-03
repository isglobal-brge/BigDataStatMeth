#include <BigDataStatMeth.hpp>
// #include "hdf5Utilities/hdf5RemoveLowData.hpp"


//' Remove SNPs in hdf5 omic dataset with low data
//'
//' Remove SNPs in hdf5 omic dataset with low data
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param dataset, character array indicating the input dataset to be imputed
//' @param outgroup, character array indicating group where the data set will be 
//' saved after remove data with if `outgroup` is NULL, output dataset is stored 
//' in the same input group. 
//' @param outdataset, character array indicating dataset to store the resulting 
//' data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
//' @param pcent, by default pcent = 0.5. Numeric indicating the percentage to be 
//' considered to remove SNPs, SNPS with percentage equal or higest will be removed from data
//' @param bycols, boolean by default = true, if true, indicates that SNPs are in 
//' cols, if SNPincols = false indicates that SNPs are in rows.
//' @param overwrite, optional boolean if true, previous results in same location 
//' inside hdf5 will be overwritten, by default overwrite = false, data was not overwritten.
//' @return Original hdf5 data file without cols/rows with low represented snps
//' @examples
//' print('see vignette')
//' @export
// [[Rcpp::export]]
void bdRemovelowdata_hdf5( std::string filename, std::string group, std::string dataset, std::string outgroup, std::string outdataset, 
                               Rcpp::Nullable<double> pcent, Rcpp::Nullable<bool> bycols, Rcpp::Nullable<bool> overwrite = R_NilValue)
{
    
    try
    {
        
        BigDataStatMeth::hdf5Dataset* dsIn;
        BigDataStatMeth::hdf5DatasetInternal* dsOut;
        
        bool bcols, bforce;
        double dpcent;
        int iremoved = 0;
        
        std::string stroutdata = outgroup +"/" + outdataset;
        std::string strdataset = group +"/" + dataset;
        
        if(bycols.isNull()){  bcols = true ; }
        else{  bcols = Rcpp::as<bool>(bycols); }
        
        if(pcent.isNull()){   dpcent = 0.5 ; }
        else{     dpcent = Rcpp::as<double>(pcent); }
        
        if(overwrite.isNull()) { bforce = false ; }
        else { bforce = Rcpp::as<bool>(overwrite); }
        

        if( strdataset.compare(stroutdata)!= 0)
        {
            dsIn = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
            dsIn->openDataset();
            
            dsOut = new BigDataStatMeth::hdf5DatasetInternal(filename, outgroup, outdataset, bforce);
            
        } else {
            throw std::range_error("Input and output dataset must be different");  
            return void();
        }
        
        iremoved = Rcpp_Remove_Low_Data_hdf5( dsIn, dsOut, bcols, dpcent);
        
        Rcpp::Function warning("warning");
        if (bycols )
            warning( std::to_string(iremoved) + " Columns have been removed");
        else
            warning( std::to_string(iremoved) + " Rows have been removed");
    
        delete dsIn;
        delete dsOut;
        
    } catch( H5::FileIException& error ){ // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception bdRemovelowdata_hdf5 (File IException)" );
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception bdRemovelowdata_hdf5 (DataSet IException)" );
        return void();
    } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception bdRemovelowdata_hdf5 (DataSpace IException)" );
        return void();
    } catch( H5::DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception bdRemovelowdata_hdf5 (DataType IException)" );
        return void();
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return void();
    }
    
    return void();
}
