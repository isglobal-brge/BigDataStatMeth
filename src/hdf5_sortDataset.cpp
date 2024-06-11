#include <BigDataStatMeth.hpp>
#include "hdf5Utilities/hdf5SortDataset.hpp"

//' Sort existing dataset 
//'
//' Sort an existing dataset taking in to account a list with sorted positions
//' 
//' @param filename, character array indicating the name of the file to be sorted
//' @param group, character array indicating the input group where the data set 
//' to be sorted is stored.
//' @param dataset, character array indicating the input dataset to be sorted
//' @param outdataset, character array indicating the name for the new sorted 
//' dataset. This dataset 
//' @param blockedSortlist, a list with blocks with sorted positions, see example
//' $`1`
//'                       chr order newOrder Diagonal
//' TCGA-OR-A5J1 TCGA-OR-A5J1     1        1        1
//' TCGA-OR-A5J2 TCGA-OR-A5J2     2        2        1
//' TCGA-OR-A5J3 TCGA-OR-A5J3     3        3        1
//' TCGA-OR-A5J4 TCGA-OR-A5J4     4        4        1
//' 
//' $`2`
//'                       chr order newOrder
//' TCGA-OR-A5J5 TCGA-OR-A5JA    10        5        1
//' TCGA-OR-A5J6 TCGA-OR-A5JB    11        6        1
//' TCGA-OR-A5J7 TCGA-OR-A5JC    12        7        0
//' TCGA-OR-A5J8 TCGA-OR-A5JD    13        8        1
//' 
//' $`3`
//'                       chr order newOrder
//' TCGA-OR-A5J9 TCGA-OR-A5J5     5        9        1
//' TCGA-OR-A5JA TCGA-OR-A5J6     6       10        1
//' TCGA-OR-A5JB TCGA-OR-A5J7     7       11        1
//' TCGA-OR-A5JC TCGA-OR-A5J8     8       12        1
//' TCGA-OR-A5JD TCGA-OR-A5J9     9       13        0
//' 
//' where rowname is the current rowname, chr is the new rowname, order is the
//' current position and newOrder is the new position
//' @param func, character array function to be applyed
//' \describe{
//'     \item{sortRows}{sort datasets rows}
//'     \item{sortCols}{sort datasets columns}
//' }
//' @param outgroup, optional, character array indicating group where the data 
//' set will be saved after imputation if `outgroup` is NULL, output dataset is 
//' stored in the same input group. 
//' @param overwrite, boolean if true, previous results in same location inside hdf5
//' will be overwritten.
//' @return Original hdf5 data file with sorted dataset
//' @examples
//' 
//' print("See vignette")
//' @export
// [[Rcpp::export]]
void bdSort_hdf5_dataset( std::string filename, std::string group, 
                          std::string dataset, std::string outdataset, 
                          Rcpp::List blockedSortlist, std::string func, 
                          Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                          Rcpp::Nullable<bool> overwrite = false )
{
    
    try
    {
        
        BigDataStatMeth::hdf5Dataset* dsIn;
        BigDataStatMeth::hdf5Dataset* dsOut;
        
        std::string strOutgroup;
        bool boverwrite;
        hsize_t ncols = 0,
                nrows = 0;
        
        if( blockedSortlist.length()<=0 ) {
            Rcpp::Rcout<<"\nList is empty, please create a list with the new sort";
            return void();
        }
        
        if( overwrite.isNull() ) { boverwrite = false; } 
        else { boverwrite = Rcpp::as<bool>(overwrite); }
        
        if( outgroup.isNull() ) {  strOutgroup = group;  } 
        else {  strOutgroup = Rcpp::as<std::string>(outgroup); }
        
        dsIn = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
        dsIn->openDataset();
        
        ncols = dsIn->ncols();

        // Get the nomber of rows in dataframes inside the list
        for(int i=0; i<blockedSortlist.size(); i++) {     
            Rcpp::DataFrame df(blockedSortlist[i]);
            nrows = nrows + df.nrow();
        } 
        
        dsOut = new BigDataStatMeth::hdf5Dataset(filename, strOutgroup, outdataset, boverwrite);
        dsOut->createDataset( ncols, nrows, "real");
        
        RcppSort_dataset_hdf5(dsIn, dsOut, blockedSortlist, func);
        
        delete dsIn;
        delete dsOut;
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception bdSort_hdf5_dataset (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdSort_hdf5_dataset (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdSort_hdf5_dataset (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdSort_hdf5_dataset" << ex.what();
        return void();
    }
    
    Rcpp::Rcout<<outdataset<<" dataset has been sorted \n";
    return void();
    
}
