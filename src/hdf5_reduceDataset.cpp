#include <BigDataStatMeth.hpp>
#include "hdf5Utilities/hdf5ReduceDataset.hpp"


//' Reduce hdf5 dataset
//'
//' Reduce hdf5 datasets inside a group by rows or columns and store complete matrix inside the hdf5 data file.
//' 
//' @param filename, character array with the name of the file where datasets are stored
//' @param group, character array with the input group name where the data sets are stored 
//' @param reducefunction, single character with function to apply, can be '+' or  '-'
//' @param outgroup, optional character array with the group name where the dataset will be stored after reduction. If `outgroup` is NULL, the resulting dataset is stored in the same input group. 
//' @param outdataset, optional character with the dataset name for the resulting dataset after reduction if `outdataset` is NULL, then the input group name is used as outdataset name 
//' @param overwrite, boolean if true, previous results in same location inside hdf5 will be overwritten.
//' @param remove, boolean if true, removes original matrices, by default bremove = false.
//' @return Full matrix with results from reduction
//' @export
// [[Rcpp::export]]
void bdReduce_hdf5_dataset( std::string filename, std::string group, 
                            std::string reducefunction, 
                            Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                            Rcpp::Nullable<std::string> outdataset = R_NilValue,
                            Rcpp::Nullable<bool> overwrite = false , 
                            Rcpp::Nullable<bool> remove = false )
{
    
    try
    {
        // BigDataStatMeth::hdf5Dataset* dsOut;
        
        std::string strOutgroup,
                    strOutdatset,
                    stroutdata;
        
        bool boverwrite, bremove;
        
        if (reducefunction.compare("+") != 0 && reducefunction.compare("-") != 0 ) {
            throw std::range_error( "Function to apply must be \"+\" or \"-\" other values are not allowed" );
            return void();
        }
        
        if(outgroup.isNull()) {  strOutgroup = group ;
        } else {   strOutgroup = Rcpp::as<std::string>(outgroup);}
        
        if(remove.isNull()) { bremove = false ;
        } else {   bremove = Rcpp::as<bool>(remove);}
        
        if(overwrite.isNull()) { boverwrite = false; } 
        else {   boverwrite = Rcpp::as<bool>(overwrite); }
        
        if(outdataset.isNull()){  strOutdatset = group ;
        } else {   strOutdatset = Rcpp::as<std::string>(outdataset);}
        
        
        BigDataStatMeth::RcppReduce_dataset_hdf5( filename, group, strOutgroup, strOutdatset, reducefunction, boverwrite, bremove, false);
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception bdReduce_hdf5_dataset (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdReduce_hdf5_dataset (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdReduce_hdf5_dataset (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdReduce_hdf5_dataset" << ex.what();
        return void();
    }
    
    return void();
}
