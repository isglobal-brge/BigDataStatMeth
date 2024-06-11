#include <BigDataStatMeth.hpp>

//' Gets all dataset names inside a group
//'
//' Gets a list of all dataset names inside a group or all the datasets names 
//' starting with a prefix under a group
//' 
//' @param filename, character array with the name of the file to be accessed
//' @param group, character array with the input group name where the data sets are stored 
//' @param prefix, character array optional, indicates the prefix with which the dataset
//' names begin, if null, then the function returns all datasets inside the group
//' @return character array with the name of all datasets inside the group
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdgetDatasetsList_hdf5(std::string filename, std::string group, Rcpp::Nullable<std::string> prefix = R_NilValue)
{
    
    // H5File* file = nullptr;
    Rcpp::StringVector groupDatasets;
    
    try
    {
        std::string strprefix;
        
        if(prefix.isNull()){  strprefix = "" ;
        } else {   strprefix = Rcpp::as<std::string>(prefix);}
        
        BigDataStatMeth::hdf5File* fQuery = new BigDataStatMeth::hdf5File(filename, false);
        fQuery->openFile("r");
        
        // Get dataset names without prefix, all datasets inside the group
        groupDatasets =  fQuery->getDatasetNames(group, strprefix);
        
        delete fQuery;
        
    }catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception bdgetDatasetsList_hdf5 (File IException)" );
        return(R_NilValue);
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception bdgetDatasetsList_hdf5 (DataSet IException)" );
        return(R_NilValue);
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return(R_NilValue);
    }
    
    return(groupDatasets);
    
}




/*** R

*/

