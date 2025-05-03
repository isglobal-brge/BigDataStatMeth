#include <BigDataStatMeth.hpp>

//' Gets all dataset names inside a group
//'
//' Gets a list of all dataset names inside a group or all the datasets names 
//' starting with a prefix under a group
//' 
//' @param filename, character array with the name of the file to be accessed
//' @param group, character array with the input group name where the data sets 
//' are stored 
//' @param prefix, character array optional, indicates the prefix with which the 
//' dataset names begin, if prefix and sufix are null, the function returns all 
//' datasets inside the group
//' @param sufix, character array optional, indicates the sufix with which the 
//' dataset names ends, if prefix and sufix are null, the function returns all 
//' datasets inside the group
//' @return character array with the name of all datasets inside the group
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdgetDatasetsList_hdf5(std::string filename, std::string group, Rcpp::Nullable<std::string> prefix = R_NilValue)
{
    
    // H5File* file = nullptr;
    Rcpp::StringVector groupDatasets;
    BigDataStatMeth::hdf5File* fQuery = nullptr;
    
    try
    {
        std::string strprefix;
        
        if(prefix.isNull()){  strprefix = "" ;
        } else {   strprefix = Rcpp::as<std::string>(prefix);}
        
        fQuery = new BigDataStatMeth::hdf5File(filename, false);
        fQuery->openFile("r");
        
        // Get dataset names without prefix, all datasets inside the group
        groupDatasets =  fQuery->getDatasetNames(group, strprefix, "");
        
        delete fQuery; fQuery = nullptr;
        
    } catch( H5::FileIException& error ) { 
        fQuery = nullptr;
        Rcpp::Rcerr<<"\nc++ exception bdgetDatasetsList_hdf5 (File IException)\n";
        return(R_NilValue);
    } catch( H5::DataSetIException& error ) { 
        delete fQuery; fQuery = nullptr;
        Rcpp::Rcerr<<"\nc++ exception bdgetDatasetsList_hdf5 (DataSet IException)\n";
        return(R_NilValue);
    } catch(std::exception &ex) {
        delete fQuery; fQuery = nullptr;
        Rcpp::Rcout<< ex.what();
        return(R_NilValue);
    }
    
    return(groupDatasets);
    
}




/*** R

*/

