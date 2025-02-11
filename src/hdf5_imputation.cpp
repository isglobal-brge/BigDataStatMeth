#include <BigDataStatMeth.hpp>
#include "hdf5Utilities/hdf5ImputeData.hpp"


//' Impute SNPs in hdf5 omic dataset 
//'
//' Impute SNPs in hdf5 omic dataset 
//' 
//' @inheritParams bdblockmult_hdf5
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
                        Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                        Rcpp::Nullable<std::string> outdataset = R_NilValue, 
                        Rcpp::Nullable<bool> bycols = true, 
                        Rcpp::Nullable<bool> paral = R_NilValue,
                        Rcpp::Nullable<int> threads = R_NilValue, 
                        Rcpp::Nullable<bool> overwrite = R_NilValue )
{
    
    BigDataStatMeth::hdf5Dataset* dsIn = nullptr;
    BigDataStatMeth::hdf5DatasetInternal* dsOut = nullptr;
    
    try
    {
        
        std::string strdataset = group +"/" + dataset;
        std::string stroutgroup, stroutdataset, stroutdata;
        
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
        
        if( dsIn->getDatasetptr() != nullptr ) {
            dsOut = new BigDataStatMeth::hdf5DatasetInternal(filename, stroutgroup, stroutdataset, bforce);
            Rcpp_Impute_snps_hdf5( dsIn, dsOut, bcols, stroutdataset, threads);
            delete dsOut; dsOut = nullptr;
        }
        
        delete dsIn; dsIn = nullptr;
        
    } catch( H5::FileIException& error ){
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdImputeSNPs_hdf5 (File IException)\n";
        return void();
    } catch( H5::DataSetIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdImputeSNPs_hdf5 (DataSet IException)\n";
        return void();
    } catch( H5::DataSpaceIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdImputeSNPs_hdf5 (DataSpace IException)\n";
        return void();
    } catch( H5::DataTypeIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdImputeSNPs_hdf5 (DataType IException)\n";
        return void();
    } catch(std::exception &ex) {
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdImputeSNPs_hdf5: "<< ex.what()<<"\n";
        return void();
    }  catch (...) {
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nC++ exception bdImputeSNPs_hdf5 (unknown reason)";
        return void();
    }
    
    // Rcpp::Rcout<<"SNPs with missing values has been imputed\n";
    return void();
    
}

