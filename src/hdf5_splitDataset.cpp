#include <BigDataStatMeth.hpp>
// #include "hdf5Utilities/hdf5SplitDataset.hpp"

//' Split hdf5 dataset
//'
//' Split hdf5 dataset in small datasets by rows or columns and store splitted 
//' submatrices inside an hdf5 file.
//' 
//' @param filename, character array indicating the name of the file where 
//' dataset to split is stored
//' @param group, character array indicating the input group where the data set 
//' to be splitted is. 
//' @param dataset, character array indicating the input dataset to be splitted
//' @param outgroup, optional character array indicating group where the data set 
//' will be saved after split process if `outgroup` is NULL, output dataset is 
//' stored in the same input group. 
//' @param outdataset, optional character array indicating basename for the 
//' splitted dataset if `outdataset` is NULL, input dataset name is used adding .x, 
//' where x is the splitted block number. 
//' @param nblocks, integer number of blocks in which we want to split the data
//' @param blocksize, integer, number of elements in each block
//' @param bycols, boolean by default = true, true indicates that the imputation 
//' will be done by columns, otherwise, the imputation will be done by rows
//' @param overwrite, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
//' @return Splitted datasets inside an hdf5 data file
//' @export
// [[Rcpp::export]]
void bdSplit_matrix_hdf5( std::string filename, std::string group, std::string dataset, 
                          Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                          Rcpp::Nullable<std::string> outdataset = R_NilValue, 
                          Rcpp::Nullable<int> nblocks = R_NilValue,  
                          Rcpp::Nullable<int> blocksize = R_NilValue,
                          Rcpp::Nullable<bool> bycols = true, 
                          Rcpp::Nullable<bool> overwrite = false  )
{
    
    BigDataStatMeth::hdf5Dataset* dsIn;
    std::string stroutgroup, stroutdataset, stroutdata;
    
    try
    {
        std::string strdataset = group + "/" + dataset;
        std::string strdatasetout;
        int iblocksize = 0; //, iwholesize = 0;
        bool bcols; //, bforce;
        
        if(bycols.isNull()) { bcols = true ;
        } else {   bcols = Rcpp::as<bool>(bycols);}
        
        // if(overwrite.isNull()) { bforce = false ;
        // } else {   bforce = Rcpp::as<bool>(overwrite);}
        
        if(outgroup.isNull()) {  stroutgroup = group ;
        } else {   stroutgroup = Rcpp::as<std::string>(outgroup);}
        
        if(outdataset.isNull()){  stroutdataset = dataset ;
        } else {   stroutdataset = Rcpp::as<std::string>(outdataset);}
        
        dsIn = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
        dsIn->openDataset();
        
        if( dsIn->getDatasetptr() != nullptr ) { 
            
            hsize_t nrows = dsIn->nrows(),
                ncols = dsIn->ncols();
            
            
            if( nblocks.isNull() && blocksize.isNull()){
                delete dsIn;
                Rcpp::Rcerr<<"\n Block size or number of blocks are needed to proceed with matrix split. Please, review parameters";
                return void();
                
            } else if (!nblocks.isNull() && !blocksize.isNull()) {
                delete dsIn;
                Rcpp::Rcerr<<"\nBlock size and number of blocks are defined, please define only one option, split by number of blocks or by block size";
                return void();
                
            } else if(!nblocks.isNull()) {
                
                if ( Rcpp::as<int>(nblocks) == 1) {
                    delete dsIn;
                    Rcpp::Rcerr<<"\nNumbers of blocks = 1, no data to split";
                    return void();
                    
                } else {
                    
                    double module;
                    
                    if(bcols == true) {
                        iblocksize = nrows / Rcpp::as<int>(nblocks);
                        module = nrows % iblocksize;
                    } else {
                        iblocksize = ncols / Rcpp::as<int>(nblocks);
                        module = ncols % iblocksize;
                    }
                    if (module > 0) { iblocksize = iblocksize + 1; }
                }
                
            } else {
                iblocksize = Rcpp::as<int>(blocksize);
                
                if( bcols == true ) {
                    if( iblocksize == nrows) {  
                        delete dsIn;
                        throw std::range_error( "c++ exception in bdSplit_matrix_hdf5 ( No data to split)" ); 
                    }
                } else {
                    if( iblocksize == ncols) {  
                        delete dsIn;
                        throw std::range_error( "c++ exception in bdSplit_matrix_hdf5 ( No data to split)" ); 
                    }
                }
            }
            
            RcppSplit_matrix_hdf5 ( dsIn, bcols, stroutgroup, stroutdataset, iblocksize, ncols, nrows );    
            
        }
    
        delete dsIn; dsIn = nullptr;
        Rcpp::Rcout<<"Dataset has been splitted, results can be found in "<< stroutgroup + "/" + stroutdataset <<"\n";
        
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        checkClose_file(dsIn);
        Rcpp::Rcerr<<"\nc++ exception bdSplit_matrix_hdf5 (File IException)\n";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        checkClose_file(dsIn);
        Rcpp::Rcerr<<"\nc++ exception bdSplit_matrix_hdf5 (DataSet IException)";
        return void();
    } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        checkClose_file(dsIn);
        Rcpp::Rcerr<<"\nc++ exception bdSplit_matrix_hdf5 (DataSpace IException)";
        return void();
    } catch( H5::DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        checkClose_file(dsIn);
        Rcpp::Rcerr<<"\nc++ exception bdSplit_matrix_hdf5 (DataType IException)";
        return void();
    } catch(std::exception &ex) {
        checkClose_file(dsIn);
        Rcpp::Rcerr<<"\nC++ exception bdSplit_matrix_hdf5 : "<< ex.what();
    } catch (...) {
        checkClose_file(dsIn);
        Rcpp::Rcerr<<"\nC++ exception bdSplit_matrix_hdf5 (unknown reason)";
        return void();
    }
    
    return void();
}

