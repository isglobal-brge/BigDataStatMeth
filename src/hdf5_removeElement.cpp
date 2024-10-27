#include <BigDataStatMeth.hpp>
#include "hdf5Utilities/hdf5RemoveElements.hpp"



//' Remove element group or dataset from  hdf5 file
//'
//' Remove group or dataset from  hdf5 file
//' 
//' @param filename, character array indicating the name of the file to create
//' @param element string vector with one or multiple elements to be removed, 
//' each element in the string vectur must be a complete route to the element to be removed.
//' @return none
//' @export
//' 
//' @examples
//' 
//' matA <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), nrow = 3, byrow = TRUE)
//' matB <- matrix(c(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,5,3,4,5,2,6,2,3,4,
//'                    42, 23, 23, 423,1,2), nrow = 3, byrow = TRUE)
//'                    
//' bdCreate_hdf5_matrix("BasicMatVect.hdf5", matA, "INPUT", "matA")
//' bdCreate_hdf5_matrix("BasicMatVect.hdf5", matB, "INPUT", "matB")
//' 
//' bdRemove_hdf5_element("BasicMatVect.hdf5", c("INPUT/matA", "INPUT/matB"))
//' 
//' 
//' # Remove file (used as example)
//'   if (file.exists("BasicMatVect.hdf5")) {
//'     # Delete file if it exist
//'     file.remove("BasicMatVect.hdf5")
//'   }
//' 
// [[Rcpp::export]]
void bdRemove_hdf5_element(std::string filename, std::vector<std::string> elements)
{
    
    bool allOk;
    try
    {
        
        BigDataStatMeth::hdf5File* objFile = new BigDataStatMeth::hdf5File(filename, false);
        objFile->openFile("rw");
        
        BigDataStatMeth::RcppRemove_hdf5_elements(objFile, elements);
        
        delete objFile;
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception bdBind_hdf5_datasets (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdBind_hdf5_datasets (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdBind_hdf5_datasets (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdBind_hdf5_datasets" << ex.what();
        return void();
    }
    
    // if(allOk == false)
    //     std::string strmessage = "Not all elements had been removed";
    // 
    // Rcpp::message(Rcpp::wrap(strmessage));    
    return void();
    
}
