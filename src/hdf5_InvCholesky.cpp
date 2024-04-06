#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixInvCholesky.hpp"
#include "hdf5Algebra/matrixTriangular.hpp"



//' Compute inverse cholesky with hdf5 data files
//'
//' Compute inverse cholesky with datasets stored in hdf5 data files
//'
//' @param filename, character array with the name of an existin hdf5 data file containing the dataset to be modified
//' @param group, character array indicating the input group where the data set to be modified. 
//' @param dataset, character array indicating the input dataset to be modified
//' @param outdataset character array with output dataset name where we want to store results
//' @param outgroup optional, character array with output group name where we want to 
//' store results if not provided then results are stored in the same group as original dataset
//' @param fullMatrix boolean, optional parameter, by default false. 
//' If fullMatrix = true, in the hdf5 file the complete matrix is stored. 
//' If false, only the lower triangular matrix is saved
//' @param force, optional boolean if true, previous results in same location inside 
//' hdf5 will be overwritten, by default force = false, data was not overwritten.
//' @param threads optional parameter. Integer with numbers of threads to be used
//' @param elementsBlock, optional integer defines de maximum number of elements 
//' to read from hdf5 data file in each block. By default this value is set 
//' to 10000. If matrix is bigger thant 5000x5000 then block is set to number 
//' of rows or columns x 2
//' @return Original hdf5 data file with Inverse of Cholesky
//' @examples
//' 
//' library(BigDataStatMeth)
//' library(rhdf5)
//' 
//' set.seed(1234)
//'     Y <- matrix(sample.int(10, 100, replace = TRUE), ncol = 10)
//'     
//' # devtools::reload(pkgload::inst("BigDataStatMeth"))
//' Ycp <- crossprod(Y)
//'         
//' # devtools::reload(pkgload::inst("BigDataStatMeth"))
//' bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
//'                         object = Ycp, group = "data", dataset = "matrix",
//'                         transp = FALSE,
//'                         overwriteFile = TRUE, overwriteDataset = TRUE, 
//'                         unlimited = FALSE)
//'        
//' # Get Inverse Cholesky
//'        res <- bdInvCholesky_hdf5(filename = "test_temp.hdf5", group = "data", 
//'        dataset = "matrix", outdataset = "invmatrix", outgroup = "InvCholesky", 
//'        fullMatrix = FALSE, force = TRUE)
//' 
//' @export
// [[Rcpp::export]]
void bdInvCholesky_hdf5( std::string filename, std::string group, std::string dataset,
                         std::string  outdataset,
                         Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                         Rcpp::Nullable<bool> fullMatrix = R_NilValue, 
                         Rcpp::Nullable<bool> force = R_NilValue,
                         Rcpp::Nullable<int> threads = 2,
                         Rcpp::Nullable<long> elementsBlock = 1000000)
{
    
    long dElementsBlock;
    // std::string strOutgroup, strIndataset, 
    //             strOutdataset, strOutdataset_tmp;
    
    std::string strOutgroup, strOutdataset;
    // int ithreads;
    // bool bforce, bfull;
    bool bfull;
    int nrows = 0, ncols = 0;
    
    try
    {
        
        // Get default values for Nullable variables
        // if(force.isNull()) { bforce = false; } 
        // else {  bforce = Rcpp::as<bool>(force); }
        
        if(fullMatrix.isNull()) { bfull = false; } 
        else {  bfull = Rcpp::as<bool>(fullMatrix); }
        
        if(outgroup.isNull()) { strOutgroup = group; } 
        else {   strOutgroup = Rcpp::as<std::string>(outgroup); }
        
        // if(threads.isNull()) { ithreads = 2; } 
        // else { ithreads = Rcpp::as<int>(threads); }
        
        if(elementsBlock.isNull()) { dElementsBlock = MAXELEMSINBLOCK; } 
        else { dElementsBlock = Rcpp::as<long>(elementsBlock); }
        
        
        // strIndataset = group + "/" + dataset;
        strOutdataset = strOutgroup + "/" + outdataset;
        // strOutdataset_tmp = "tmp/tmp_L";

        BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
        dsA->openDataset();
        
        nrows = dsA->nrows();
        ncols = dsA->ncols();
        
        
        if(nrows == ncols) {

            BigDataStatMeth::hdf5DatasetInternal* dstmp = new BigDataStatMeth::hdf5DatasetInternal(filename, strOutdataset, true);
            dstmp->createDataset(nrows, ncols, "real");
            
            BigDataStatMeth::Rcpp_InvCholesky_hdf5( dsA, dstmp, bfull, dElementsBlock, threads);
            
            // int res = Cholesky_decomposition_hdf5(dsA, dstmp, nrows, ncols, dElementsBlock, threads);
            // delete dsA;
            // 
            // if(res == 0)
            // {
            //     Inverse_of_Cholesky_decomposition_hdf5( dstmp, nrows, ncols, dElementsBlock, threads);
            //     Inverse_Matrix_Cholesky_parallel( dstmp, nrows, ncols, dElementsBlock, threads); 
            //     
            //     if( bfull==true ) {
            //         setUpperTriangularMatrix( dstmp, dElementsBlock);
            //     }
            // }
            
            delete dsA;
            delete dstmp;
            
        } else {
            delete dsA;
            Rcpp::Rcout<<"\n Can't get inverse matrix for "<< group + "/" + dataset <<" using Cholesky decomposition \n";
            return void();
        }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception bdInvCholesky_hdf5 (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdInvCholesky_hdf5 (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdInvCholesky_hdf5 (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdInvCholesky_hdf5" << ex.what();
        return void();
    }
    
    return void();
}
