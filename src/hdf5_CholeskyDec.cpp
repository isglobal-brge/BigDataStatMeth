#include <BigDataStatMeth.hpp>
// #include "hdf5Algebra/matrixInvCholesky.hpp"
// #include "hdf5Algebra/matrixTriangular.hpp"



//' Compute Cholesky decomposition with hdf5 data files
//'
//' Compute cholesky decomposition with datasets stored in hdf5 data files. Function returns the upper triangular matrix.
//'
//' @inheritParams bdblockmult_hdf5
//' @inheritParams bdNormalize_hdf5
//' @param fullMatrix boolean, optional parameter, by default false. 
//' If fullMatrix = true, in the hdf5 file the complete matrix is stored. 
//' If false, only the lower triangular matrix is saved
//' @param elementsBlock integer (optional), an integer that specifies the 
//' maximum number of elements to read from the HDF5 data file in each block. 
//' By default, this value is set to 100,000. If the matrix size exceeds 5000x5000, 
//' the block size is automatically adjusted to number of rows or columns * 2
//' @return Original hdf5 data file with Cholesky decomposition
//' 
//' @details 
//' The **Cholesky decomposition** is a factorization of a **symmetric positive-definite matrix** \eqn{A} 
//' into the product of a **lower triangular matrix** \eqn{L} and its transpose.
//' \deqn{A = L L^\top}
//' where:
//'   * \eqn{A} is a symmetric positive-definite matrix of size \eqn{n \times n},
//'   * \eqn{L} is a lower triangular matrix of size \eqn{n \times n},
//'   * \eqn{L^\top} is the transpose of \eqn{L}.
//'   
//' ### Key Properties
//' 1. **Positive-Definiteness**: The matrix \eqn{A} must be positive-definite
//' 2. **Positive Diagonal**: The diagonal elements of \eqn{L}, \eqn{l_{ii}}, are strictly positive.
//' 
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
//' bdCholesky_hdf5(filename = "test_temp.hdf5", group = "data", 
//'    dataset = "matrix", outdataset = "matrixDec", outgroup = "Cholesky_Dec", 
//'    fullMatrix = FALSE, overwrite = TRUE)
//'        
//' res <-  h5read("test_temp.hdf5", "Cholesky_Dec/matrixDec")
//' 
//' @export
// [[Rcpp::export]]
 void bdCholesky_hdf5( std::string filename, std::string group, std::string dataset,
                          std::string  outdataset,
                          Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                          Rcpp::Nullable<bool> fullMatrix = R_NilValue, 
                          Rcpp::Nullable<bool> overwrite = R_NilValue,
                          Rcpp::Nullable<int> threads = R_NilValue,
                          Rcpp::Nullable<long> elementsBlock = 1000000)
 {
     
     long dElementsBlock;
     std::string strOutgroup, strIndataset, 
     strOutdataset, strOutdataset_tmp;
     
     BigDataStatMeth::hdf5Dataset* dsA = nullptr;
     BigDataStatMeth::hdf5DatasetInternal* dstmp = nullptr;
     
     int nrows = 0, ncols = 0;
     
     try
     {
         
         if(outgroup.isNull()) { strOutgroup = group; } 
         else {   strOutgroup = Rcpp::as<std::string>(outgroup); }
         
         if(elementsBlock.isNull()) { dElementsBlock = MAXELEMSINBLOCK; } 
         else { dElementsBlock = Rcpp::as<long>(elementsBlock); }
         
         
         strIndataset = group + "/" + dataset;
         strOutdataset = strOutgroup + "/" + outdataset;
         strOutdataset_tmp = "tmp/tmp_L";
         
         dsA = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
         dsA->openDataset();
         
         if( dsA->getDatasetptr() != nullptr) { 
             nrows = dsA->nrows();
             ncols = dsA->ncols();
             
             if(nrows == ncols) {
                 
                 dstmp = new BigDataStatMeth::hdf5DatasetInternal(filename, strOutdataset, true);
                 dstmp->createDataset(nrows, ncols, "real");
                 
                 int res = Cholesky_decomposition_hdf5(dsA, dstmp, nrows, ncols, dElementsBlock, threads);
                 
                 if(res != 0) {
                     Rcpp::Rcout<<"\n Can't get Cholesky decomposition \n";
                 }
                 
                 delete dstmp; dstmp = nullptr;
                 
             } else {
                 Rcpp::Rcout<<"\n Can't get Cholesky decomposition \n";
             }    
         }
         
         delete dsA; dsA = nullptr;
         
         
     } catch( H5::FileIException& error ) { 
         checkClose_file(dsA, dstmp);
         Rcpp::Rcerr<<"c++ exception bdCholesky_hdf5 (File IException)";
     } catch( H5::GroupIException & error ) { 
         checkClose_file(dsA, dstmp);
         Rcpp::Rcerr << "c++ exception bdCholesky_hdf5 (Group IException)";
     } catch( H5::DataSetIException& error ) { 
         checkClose_file(dsA, dstmp);
         Rcpp::Rcerr << "c++ exception bdCholesky_hdf5 (DataSet IException)";
     } catch(std::exception& ex) {
         checkClose_file(dsA, dstmp);
         Rcpp::Rcerr << "c++ exception bdCholesky_hdf5" << ex.what();
     } catch (...) {
         checkClose_file(dsA, dstmp);
         Rcpp::Rcerr<<"\nC++ exception bdCholesky_hdf5 (unknown reason)";
     }
     
     return void();
 }


/***
 //' @param filename, character array with the name of an existin hdf5 data file containing the dataset to be modified
 //' @param group, character array indicating the input group where the data set to be modified. 
 //' @param dataset, character array indicating the input dataset to be modified
 //' @param outdataset character array with output dataset name where we want to store results
 //' @param outgroup optional, character array with output group name where we want to 
 //' store results if not provided then results are stored in the same group as original dataset
 //' @param overwrite, optional boolean if true, previous results in same location inside 
 //' hdf5 will be overwritten, by default overwrite = false, data was not overwritten.
 //' @param threads optional parameter. Integer with numbers of threads to be used
 */