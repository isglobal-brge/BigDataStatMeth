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
        
        // bool transposed, 
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
        
        // if(transp.isNull())  transposed = false ;
        // else    transposed = Rcpp::as<bool>(transp);
        
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
        
        // Rcpp::IntegerVector offset = {0, 0 };
        // Rcpp::IntegerVector stride = {1, 1 };
        // Rcpp::IntegerVector block = {1, 1 };
        
        BigDataStatMeth::hdf5Dataset* objDataset = new BigDataStatMeth::hdf5Dataset(objFile, strsubgroup, strdataset, bforceDataset );
        
        if( bunlimited == false){
            objDataset->createDataset(dims[0], dims[1], strdatatype);
        } else{
            objDataset->createUnlimitedDataset(dims[0], dims[1], strdatatype);
        }
        objDataset->writeDataset(object); 
        
        // /*** 
        // if ( object.sexp_type()==0   )
        //     throw std::range_error("Data matrix must exsits and mustn't be null");
        // 
        // std::string strdatasetUnlim = strdataset + "_unlimited";
        // 
        // Rcpp::Rcout<<"Creem l'objecte\n";
        // BigDataStatMeth::hdf5File* objFile = new BigDataStatMeth::hdf5File(filename, bforce);
        // Rcpp::Rcout<<"Objecte creat, anem a crear el fitxer\n";
        // objFile->createFile();
        // Rcpp::Rcout<<"Fitxer creat\n";
        // 
        // Rcpp::IntegerVector offset = {0, 0 };
        // Rcpp::IntegerVector stride = {1, 1 };
        // Rcpp::IntegerVector block = {1, 1 };
        // 
        // 
        // Rcpp::Rcout<<"Anem a crear el dataset\n";
        // BigDataStatMeth::hdf5Dataset* objDataset = new BigDataStatMeth::hdf5Dataset(objFile, strsubgroup, strdataset, bforce );
        // 
        // ***/
        //     
        //     
        // //..// objDataset->createDataset(nrows, ncols, "real");
        // //..// objDataset_unlim->createUnlimitedDataset(nrows, ncols, "int");
        // 
        // Rcpp::Rcout<<"Write Dataset\n";
        // // objDataset->writeDatasetBlock(object, offset, stride, block );
        // objDataset->writeDataset(object);
        // // 
        // // 
        
        // 
        // objDataset->writeDatasetBlock(object, offset, stride, block );
        // Rcpp::Rcout<<"END Write Dataset Block\n";
        // 
        // BigDataStatMeth::hdf5Dataset* objDataset_unlim = new BigDataStatMeth::hdf5Dataset(objFile, strsubgroup, strdatasetUnlim, bforce );
        // objDataset_unlim->createUnlimitedDataset(nrows, ncols, "int");
        // objDataset_unlim->extendUnlimitedDataset(50, 50);
        // 
        // 
        // // Vectors: 
        // int nrowsv = 1, 
        //     ncolsv = 5;
        // offset[0] = 0;
        // offset[1] = 0;
        // 
        // std::string strdatasetVect = strdataset + "_vector";
        // BigDataStatMeth::hdf5Dataset* objDataset_vect = new BigDataStatMeth::hdf5Dataset(objFile, strsubgroup, strdatasetVect, bforce );
        // objDataset_vect->createDataset(nrowsv, ncolsv, "int");
        // Rcpp::Rcout<<"Write Dataset VECTOR\n";
        // // objDataset_vect->writeDataset(object);
        // // objDataset_vect->writeDatasetBlock(object, offset, stride, block );
        // Rcpp::Rcout<<"END Write DatasetVECTOR\n";
        
        
        // BigDataStatMeth::hdf5Dataset* objDataset = new BigDataStatMeth::hdf5Dataset(objFile, strsubgroup, strdataset, bforce );
        // objDataset->createDataset(nrows, ncols, "string");
        // Rcpp::Rcout<<"Write ALL\n";
        // objDataset->writeDataset(object);
        // Rcpp::Rcout<<"Write BLOCK\n";
        // objDataset->writeDatasetBlock(object, offset, stride, block );
        
        
        delete objDataset;
        delete objFile;
        
        // objDataset_unlim->createUnlimitedDataset(ncols, nrows, "real");
        // Rcpp::Rcout<<"Anem a extendre l'unlimited\n";
        // objDataset_unlim->extendUnlimitedDataset(ncols, nrows);
        
        
        // if ( TYPEOF(object) == INTSXP ) {
        //     write_HDF5_matrix_from_R_ptr(file, strsubgroup + "/" + strdataset, Rcpp::as<IntegerMatrix>(object), transposed);
        // } else{
        //     write_HDF5_matrix_from_R_ptr(file, strsubgroup + "/" + strdataset, Rcpp::as<NumericMatrix>(object), transposed);
        // }
        // 
        // 
        // H5::H5File* file = nullptr;
        // int nrows, ncols;
        // Rcpp::CharacterVector svrows, svrcols;
        // Rcpp::List dimnames;
        // 
        // nrows = object.attr("dim")(1);
        // ncols = object.attr("dim")(2);
        // 
        // // Create new file
        // BigDataStatMeth::hdf5File h5file(filename, bforce);
        // h5file.createFile();
        // 
        // BigDataStatMeth::hdf5Dataset h5dataset(h5file.getFileptr(), strsubgroup, strdataset, bforce);
        // h5dataset.createDataset(nrows, ncols, TYPEOF(object));
        
    }  catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception getObjecDataType: "<<ex.what()<< " \n";
    }
        
        

    return void();
    
}


/***R

setwd("/Users/mailos/PhD/dummy")


a <- matrix(seq(1:50), nrow = 10, ncol = 5)

bdCreate_hdf5_matrix("testhdf5File.hdf5",a, "grup", "dataset", unlimited = TRUE)


*/
