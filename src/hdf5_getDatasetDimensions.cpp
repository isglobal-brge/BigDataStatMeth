#include <BigDataStatMeth.hpp>

//' Get hdf5 dataset dimensions
//'
//' Get hdf5 dataset dimensions
//' 
//' @param filename, character array with the name of the file where datasets are stored
//' @param dataset, character array with the input dataset name, with complete route inside the hdf5 file
//' @return Dataset dimensions
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject  bdgetDim_hdf5( std::string filename, std::string dataset)
{
    
    BigDataStatMeth::hdf5Dataset* ds;
    Rcpp::IntegerVector dims(4);
    
    try
    {
        
        ds = new BigDataStatMeth::hdf5Dataset(filename, dataset, false);
        ds->openDataset();
        
        if( ds->getDatasetptr() != nullptr ) { 
            dims[0] = ds->nrows_r();  
            dims[1] = ds->ncols_r();
        }
        
        delete ds;
        
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        if( ds->isOpen()) ds->close_file();
        Rcpp::Rcerr<<"\nc++ exception bdReduce_hdf5_dataset (File IException)";
        return(dims);
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        if( ds->isOpen()) ds->close_file();
        Rcpp::Rcerr<<"\nc++ exception bdReduce_hdf5_dataset (Group IException)";
        return(dims);
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        if( ds->isOpen()) ds->close_file();
        Rcpp::Rcerr<<"\nc++ exception bdReduce_hdf5_dataset (DataSet IException)";
        return(dims);
    } catch(std::exception& ex) {
        if( ds->isOpen()) ds->close_file();
        Rcpp::Rcerr<<"\nc++ exception bdReduce_hdf5_dataset" << ex.what();
        return(dims);
    }
    
    return(dims);
}
