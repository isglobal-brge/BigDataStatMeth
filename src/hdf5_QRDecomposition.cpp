#include <BigDataStatMeth.hpp>
#include "Utilities/Utilities.hpp"
#include "hdf5Algebra/matrixQR.hpp"



//' QR Decomposition 
//' 
//' This function compute QR decomposition (also called a QR factorization) 
//' of a matrix \code{A} into a product \code{A = QR} of an 
//' orthogonal matrix Q and an upper triangular matrix R.
//' 
//' @param X a real square matrix 
//' @param thin boolean, (optional) if thin = true returns Q thin decomposition else 
//' returns Q full decomposition, default is thin = false
//' @param block_size (optional) block size to perform computation
//' @param threads (optional) number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @return List with orthogonal matrix \code{Q} and upper triangular matrix \code{R} 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdQR( const Rcpp::RObject & X, 
                    Rcpp::Nullable<bool> thin = R_NilValue, 
                    Rcpp::Nullable<int> block_size = R_NilValue,
                    Rcpp::Nullable<int> threads = R_NilValue)
{
    
    
    Eigen::MatrixXd A;
    bool bthin;
    BigDataStatMeth::strQR decQR;
    
    try{
        A = Rcpp::as<Eigen::MatrixXd >(X);
    }catch(std::exception &ex) {
        A = Rcpp::as<Eigen::VectorXd >(X);
    }
    
    
    if( thin.isNull()) {
        bthin = false;
    } else {
        bthin = Rcpp::as<bool> (thin);
    }
    
    decQR = BigDataStatMeth::RcppQR(A, bthin);
    
    return Rcpp::List::create(Rcpp::Named("Q") = decQR.Q,
                              Rcpp::Named("R") = decQR.R
    );
}



//' QR Decomposition hdf5
//' 
//' This function compute QR decomposition (also called a QR factorization) 
//' of a matrix \code{A} into a product \code{A = QR} of an 
//' orthogonal matrix Q and an upper triangular matrix R.
//' @param filename, character array with the name of an existin hdf5 data file 
//' containing the dataset to be modified
//' @param group, character array indicating the input group where the data set 
//' to be modified. 
//' @param dataset, character array indicating the input dataset to be modified
//' @param outgroup optional, character array with output group name where we want to 
//' store results if not provided then results are stored in the same group as 
//' original dataset
//' @param outdataset character array with output dataset name where we want to 
//' store results, results are stored as Q.<outdataset> and R.<outdataset>
//' @param thin boolean, if thin = true returns Q thin decomposition else 
//' returns Q full decomposition, default is thin = false
//' @param block_size (optional) block size to perform computation
//' @param overwrite boolean(optional) if overwrite = TRUE, if exis 
//' @param threads (optional) number of concurrent threads in parallelization if 
//' threads is null then threads =  maximum number of threads available
//' @return List with orthogonal matrix \code{Q}  and upper triangular matrix \code{R}
//' @export
// [[Rcpp::export]]
void bdQR_hdf5( std::string filename, std::string group, std::string dataset,
                Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                Rcpp::Nullable<std::string> outdataset = R_NilValue,
                Rcpp::Nullable<bool> thin = R_NilValue, 
                Rcpp::Nullable<int> block_size = R_NilValue,
                Rcpp::Nullable<bool> overwrite = R_NilValue,
                Rcpp::Nullable<int> threads = R_NilValue)
{
    
    
    bool bthin, bforce;
    std::string strOutgroup, strOutdataset_Q, strOutdataset_R;
    int iblock_size;
    
    try {
        
        if(outgroup.isNull()) { strOutgroup = group + "/QRDec"; } 
        else {   strOutgroup = Rcpp::as<std::string>(outgroup); }
        
        if(outdataset.isNull()) { 
            strOutdataset_Q = "Q." + dataset;
            strOutdataset_R = "R." + dataset; 
        } 
        else {   
            strOutdataset_Q = "Q." + Rcpp::as<std::string>(outdataset);
            strOutdataset_R = "R." + Rcpp::as<std::string>(outdataset);
        }
        
        if( thin.isNull()) {
            bthin = false;
        } else {
            bthin = Rcpp::as<bool> (thin);
        }
        
        if(overwrite.isNull())  bforce = false ;
        else    bforce = Rcpp::as<bool>(overwrite);
        
        BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, bforce);
        dsA->openDataset();
        
        BigDataStatMeth::hdf5Dataset* dsQ = new BigDataStatMeth::hdf5Dataset(filename, strOutgroup, strOutdataset_Q, bforce);
        // dsQ->openDataset();
        
        BigDataStatMeth::hdf5Dataset* dsR = new BigDataStatMeth::hdf5Dataset(filename, strOutgroup, strOutdataset_R, bforce);
        // dsR->openDataset();
        
        RcppQRHdf5(dsA, dsQ, dsR, bthin, block_size, threads);
        
        delete dsA;
        delete dsQ;
        delete dsR;
        
        
    }catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception bdQR_hdf5 (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdQR_hdf5 (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdQR_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdQR_hdf5" << ex.what();
        return void();
    }
    
    return void();
    
    
}