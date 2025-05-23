#include <BigDataStatMeth.hpp>
// #include "hdf5Algebra/multiplication.hpp"

//' Hdf5 datasets multiplication
//'
//' The bdblockmult_hdf5 function performs block-wise matrix multiplication 
//' between two matrices stored in an HDF5 file. This approach is also efficient 
//' for large matrices that cannot be fully loaded into memory.
//' 
//' @param filename string specifying the path to the HDF5 file
//' @param group string specifying the group within the HDF5 file containing 
//' matrix A.
//' @param A string specifying the dataset name for matrix A.
//' the data matrix to be used in calculus
//' @param B string specifying the dataset name for matrix B.
//' @param groupB string, (optional), An optional string specifying the group 
//' for matrix B. Defaults to the value of `group` if not provided.
//' @param block_size integer (optional), an optional parameter specifying the 
//' block size for processing the matrices. If not provided, a default block 
//' size is used. The block size should be chosen based on the available memory 
//' and the size of the matrices
//' @param paral boolean (optional), an optional parameter to enable parallel 
//' computation. Defaults to FALSE. Set `paral = true` to force parallel execution
//' @param threads integer (optional), an optional parameter specifying the 
//' number of threads to use if paral = TRUE. Ignored if paral = FALSE.
//' @param outgroup string (optional), An optional parameger specifying the group 
//' where the output matrix will be stored. If NULL, the output will be stored 
//' in the default group "OUTPUT".
//' @param outdataset string (optional), An optional parameter specifying the 
//' dataset name for the output matrix. If NULL, the default name will be 
//' constructed as the name of dataset A concatenated with _x_ and the 
//' name of dataset B.
//' @param overwrite logical (optional), An optional parameter to indicate whether 
//' existing results in the HDF5 file should be overwritten. Defaults to FALSE. 
//' If FALSE and the dataset already exists, an error will be displayed, and 
//' no calculations will be performed. If TRUE and a dataset with the same 
//' name as specified in outdataset already exists, it will be overwritten.
//' @details
//' * The function `bdblockmult_hdf5()` is efficient for both matrices that cannot 
//' fit into memory (by processing in blocks) and matrices that can be fully 
//' loaded into memory, as it optimizes computations based on available resources.
//' * Ensure that the dimensions of `A` and `B` matrices are compatible for 
//' matrix multiplication.
//' * The `block size` should be chosen based on the available memory and 
//' the size of the matrices.
//' * If `bparal = true`, number of concurrent threads in parallelization. If 
//' `paral = TRUE` and `threads = NULL` then `threads` is set to a half of a 
//' maximum number of available threads 
//' @return a dataset inside the hdf5 data file with A*B 
//' @examples
//' library("BigDataStatMeth")
//' 
//' N = 1500; M = 1500
//' 
//' set.seed(555)
//' a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
//' b <- matrix( rnorm( N*M, mean=0, sd=1), M, N) 
//' 
//' bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
//'                      object = a, group = "groupA", 
//'                      dataset = "datasetA",
//'                      transp = FALSE,
//'                      overwriteFile = TRUE, 
//'                      overwriteDataset = FALSE, 
//'                      unlimited = FALSE)
//'                      
//' bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
//'                      object = t(b), 
//'                      group = "groupA", 
//'                      dataset = "datasetB",
//'                      transp = FALSE,
//'                      overwriteFile = FALSE, 
//'                      overwriteDataset = TRUE, 
//'                      unlimited = FALSE)
//'                      
//' # Multiply two matrix
//' bdblockmult_hdf5(filename = "test_temp.hdf5", group = "groupA", 
//'     A = "datasetA", B = "datasetB", outgroup = "results", 
//'     outdataset = "res", overwrite = TRUE ) 
//'     
//' bdblockmult_hdf5(filename = "test_temp.hdf5", group = "groupA", 
//'     A = "datasetA", B = "datasetB", outgroup = "results", 
//'     outdataset = "res", block_size = 1024, overwrite = TRUE ) 
//' 
//' @export
// [[Rcpp::export]]
void bdblockmult_hdf5(std::string filename, 
                    std::string group, 
                    std::string A, 
                    std::string B,
                    Rcpp::Nullable<std::string> groupB = R_NilValue, 
                    Rcpp::Nullable<int> block_size = R_NilValue,
                    Rcpp::Nullable<bool> paral = R_NilValue,
                    Rcpp::Nullable<int> threads = R_NilValue,
                    Rcpp::Nullable<std::string> outgroup = R_NilValue,
                    Rcpp::Nullable<std::string> outdataset = R_NilValue,
                    Rcpp::Nullable<bool> overwrite = R_NilValue)
{
    
    
    BigDataStatMeth::hdf5Dataset* dsA = nullptr;
    BigDataStatMeth::hdf5Dataset* dsB = nullptr;
    BigDataStatMeth::hdf5Dataset* dsC = nullptr;
    
    try {
        
        H5::Exception::dontPrint();  
        
        bool bforce;
        
        std::string strsubgroupOut, 
        strdatasetOut, 
        strsubgroupIn,
        strsubgroupInB;
        
        strsubgroupIn = group;
        
        if( outgroup.isNull()) { strsubgroupOut = "OUTPUT";
        } else { strsubgroupOut = Rcpp::as<std::string> (outgroup); }
        
        if(groupB.isNotNull()){ strsubgroupInB =  Rcpp::as<std::string> (groupB) ; } 
        else { strsubgroupInB =  group; }
        
        if( outdataset.isNotNull()) { strdatasetOut =  Rcpp::as<std::string> (outdataset); } 
        else { strdatasetOut =  A + "_x_" + B; }
        
        if (overwrite.isNull()) { bforce = false; } 
        else { bforce = Rcpp::as<bool> (overwrite); }
        
        dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupIn, A, false);
        dsA->openDataset();
        
        dsB = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupInB, B, false);
        dsB->openDataset();
        
        if( dsA->getDatasetptr() != nullptr && dsB->getDatasetptr() != nullptr) {
            dsC = new BigDataStatMeth::hdf5Dataset(filename, strsubgroupOut, strdatasetOut, bforce);
            BigDataStatMeth::multiplication(dsA, dsB, dsC, paral, block_size, threads);    
            
            delete dsC; dsC = nullptr;
        }
        
        delete dsA; dsA = nullptr;
        delete dsB; dsB = nullptr;
        
    }  catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        checkClose_file(dsA, dsB, dsC);
        Rcpp::Rcerr<<"\nc++ c++ exception blockmult_hdf5 (File IException)\n";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        checkClose_file(dsA, dsB, dsC);
        Rcpp::Rcerr<<"\nc++ exception blockmult_hdf5 (DataSet IException)\n";
        return void();
    } catch(std::exception &ex) {
        checkClose_file(dsA, dsB, dsC);
        Rcpp::Rcerr << "c++ exception blockmult_hdf5: " << ex.what();
        return void();
    }  catch (...) {
        checkClose_file(dsA, dsB, dsC);
        Rcpp::Rcerr<<"\nC++ exception blockmult_hdf5 (unknown reason)";
        return void();
    }
    
    // return List::create(Named("filename") = filename,
    //                     Named("dataset") = strsubgroupOut + "/" + strdatasetOut);
    return void();
}

