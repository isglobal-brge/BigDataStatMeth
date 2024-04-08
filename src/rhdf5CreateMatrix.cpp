#include <BigDataStatMeth.hpp>
// #include "hdf5Utilities/hdf5Utilities.hpp"


//' Create hdf5 data file and write data to it
//'
//' Creates a hdf5 file with numerical data matrix,
//' 
//' @param filename, character array indicating the name of the file to create
//' @param object numerical data matrix
//' @param group, character array indicating folder name to put the matrix in hdf5 file
//' @param dataset, character array indicating the dataset name to store the matix data
//' @param transp boolean, if trans=true matrix is stored transposed in hdf5 file
//' @param overwriteFile, optional boolean by default overwriteFile = false, if true and file exists, removes old file and creates a new file with de dataset data.
//' @param overwriteDataset, optional boolean by default overwriteDataset = false,  if true and dataset exists, removes old dataset and creates a new dataset.
//' @param unlimited, optional boolean by default unlimited = false, if true creates a dataset that can growth.
//' @return none
//' 
//' @examples
//' 
//' matA <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), nrow = 3, byrow = TRUE)
//' bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
//'                     object = a, group = "datasets", 
//'                     dataset = "datasetA", transp = FALSE, 
//'                     overwriteFile = TRUE, 
//'                     overwriteDataset = TRUE,
//'                     unlimited = FALSE,)
//' 
//' # Remove file (used as example)
//'   if (file.exists("test_temp.hdf5")) {
//'     # Delete file if it exist
//'     file.remove("test_temp.hdf5")
//'   }
//' 
//' @export
// [[Rcpp::export]]
void bdCreate_hdf5_matrix(std::string filename, 
                          Rcpp::RObject object, 
                          Rcpp::Nullable<std::string> group = R_NilValue, 
                          Rcpp::Nullable<std::string> dataset = R_NilValue,
                          Rcpp::Nullable<bool> transp = R_NilValue, 
                          Rcpp::Nullable<bool> overwriteFile = R_NilValue,
                          Rcpp::Nullable<bool> overwriteDataset = R_NilValue,
                          Rcpp::Nullable<bool> unlimited = R_NilValue)
{
    
    
    try
    {
        
        Rcpp::IntegerVector dims(2);
        
        std::string strsubgroup, 
                    strdataset,
                    strdatatype;
        
        bool bforceFile, 
             bforceDataset, 
             bunlimited;
        
        int iRes;
        
        if(group.isNull())  strsubgroup = "INPUT" ;
        else    strsubgroup = Rcpp::as<std::string>(group);
        
        if(dataset.isNull())  strdataset = "A" ;
        else    strdataset = Rcpp::as<std::string>(dataset);
        
        if(unlimited.isNull())  bunlimited = false ;
        else    bunlimited = Rcpp::as<bool>(unlimited);
        
        if(overwriteDataset.isNull())  bforceDataset = false ;
        else    bforceDataset = Rcpp::as<bool>(overwriteDataset);
        
        if(overwriteFile.isNull())  bforceFile = false ;
        else    bforceFile = Rcpp::as<bool>(overwriteFile);
        
        
        strdatatype = BigDataStatMeth::getObjecDataType(object);
        
        if ( object.sexp_type()==0   )
            throw std::range_error("Unknown data type");
        
        if ( object.sexp_type()==0   )
            throw std::range_error("Data matrix must exsits and mustn't be null");
        
        dims = BigDataStatMeth::getObjectDims(object, strdatatype);
        
        BigDataStatMeth::hdf5File* objFile = new BigDataStatMeth::hdf5File(filename, bforceFile);
        iRes = objFile->createFile();
        
        if(iRes == EXEC_WARNING) {
            objFile->openFile("rw");
        }
        
        BigDataStatMeth::hdf5Dataset* objDataset = new BigDataStatMeth::hdf5Dataset(objFile, strsubgroup, strdataset, bforceDataset );
        
        if( bunlimited == false){
            objDataset->createDataset(dims[0], dims[1], strdatatype);
        } else{
            objDataset->createUnlimitedDataset(dims[0], dims[1], strdatatype);
        }
        objDataset->writeDataset(object); 
        
        
        delete objDataset;
        delete objFile;
        
    }  catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception getObjecDataType: "<<ex.what()<< " \n";
    }
        
        

    return void();
    
}

