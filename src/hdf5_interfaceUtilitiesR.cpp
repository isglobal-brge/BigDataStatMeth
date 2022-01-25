#include "include/hdf5_interfaceUtilitiesR.h"
using namespace Rcpp;



//' Gets all dataset names inside a group
//'
//' Gets a list of all dataset names inside a group or all the datasets names 
//' starting with a prefix under a group
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data sets are stored 
//' @param prefix, character array optional, indicates the prefix with which the dataset
//' names begin, if null, then the function returns all datasets inside the group
//' @return Full matrix with results from reduction
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdgetDatasetsList_hdf5(std::string filename, std::string group, Rcpp::Nullable<std::string> prefix = R_NilValue)
{
    
    H5File* file;
    StringVector groupDatasets;
    
    try
    {
        std::string strprefix;
        
        if(prefix.isNull()){  strprefix = "" ;
        } else {   strprefix = Rcpp::as<std::string>(prefix);}
        
        
        if (exist_FileGroupDataset (filename, group, "")!= 0 ) {
            file = new H5File( filename, H5F_ACC_RDWR );
        }
        
        
        // Get dataset names without prefix, all datasets inside the group
        groupDatasets =  get_dataset_names_from_group(file, group, strprefix);
        
        
    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        file->close();
        ::Rf_error( "c++ exception (File IException)" );
        return(wrap(-1));
    }
    
    file->close();
    return(groupDatasets);

    
}




/*** R

*/
