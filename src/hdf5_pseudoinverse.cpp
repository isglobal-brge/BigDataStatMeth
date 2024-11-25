#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixPseudoinverse.hpp"
#include "Utilities/Utilities.hpp"



// 
//' Pseudo-Inverse
//' 
//' Compute the pseudo-inverse of a singular matrix
//' 
//' @param X Singular matrix (m x n)
//' @param threads (optional) only if bparal = true, number of concurrent 
//' threads in parallelization if threads is null then threads =  maximum number of threads available
//' @return Pseudo-inverse matrix of A
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdpseudoinv( Rcpp::RObject X,
                           Rcpp::Nullable<int> threads = R_NilValue)
{
    
    try {
        
        Eigen::MatrixXd A;
        
        try{
            A = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
        }catch(std::exception &ex) {
            A = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(X);
        }
        
        Eigen::MatrixXd pinv = BigDataStatMeth::RcppPseudoinv(&A, threads);
        
        return(Rcpp::wrap(pinv));
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return Rcpp::wrap(-1);
    }
    
}



// 
//' Pseudo-Inverse
//' 
//' Compute the pseudo-inverse of a singular matrix in hdf5 data files
//' 
//' @param filename, character array with the name of an existin hdf5 data file containing the dataset
//' @param group, character array indicating the input group where the dataset is stored
//' @param dataset, character array indicating the input dataset that contains an m x n singular matrix
//' @param outgroup optional, character array with output group name where we want to 
//' store results if not provided then results are stored in the same group as original dataset
//' @param outdataset character array with output dataset name where we want to store results
//' @param overwrite, optional boolean if true, previous results in same location inside 
//' hdf5 will be overwritten, by default force = false, data was not overwritten.
//' @param threads optional parameter. Integer with numbers of threads to be used
//' @return Pseudo-inverse matrix of A
//' @export
// [[Rcpp::export]]
void bdpseudoinv_hdf5(std::string filename, std::string group, std::string dataset,
                               Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                               Rcpp::Nullable<std::string> outdataset = R_NilValue, 
                               Rcpp::Nullable<bool> overwrite = R_NilValue,
                               Rcpp::Nullable<int> threads = R_NilValue)
{
     
    try {
         
        Eigen::MatrixXd A;
        std::string strOutgroup, strOutdataset;
        bool bforce;
        
        if(outgroup.isNull()) { strOutgroup = "PseudoInverse"; } 
        else {   strOutgroup = Rcpp::as<std::string>(outgroup); }
        
        if(outdataset.isNull()) { strOutdataset = dataset; } 
        else { strOutdataset = Rcpp::as<std::string>(outdataset); }
        
        if(overwrite.isNull()) { bforce = false ; }
        else { bforce = Rcpp::as<bool>(overwrite); }
        
        BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
        dsA->openDataset();
        
        BigDataStatMeth::hdf5Dataset* dsRes = new BigDataStatMeth::hdf5Dataset(filename, strOutgroup, strOutdataset, bforce);
        
        RcppPseudoinvHdf5(dsA, dsRes, threads);
        
        delete dsA;
        delete dsRes;
         
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception bdCholesky_hdf5 (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdCholesky_hdf5 (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdCholesky_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdCholesky_hdf5" << ex.what();
        return void();
    }
    
    return void();
}