#include "include/hdf5_blockSum.h"




// IMPORTANT : R data stored as Row-major in hdf5 file (stored transposed data from R)
// 
// This function differs from "hdf5_block_matrix_sum_hdf5_transposed" in parameters, 
// 
// Working directly with C-matrix in hdf5 file
// browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
//             else, indicates that R data is stored in hdf5 as column major
void hdf5_block_matrix_sum_hdf5_indatasets_transposed( std::string matA, std::string matB, 
                                                      IntegerVector sizeA, IntegerVector sizeB, int hdf5_block, 
                                                      std::string filename, std::string strsubgroupIN, 
                                                      std::string strsubgroupINB, std::string strsubgroupOUT, 
                                                      std::string strdatasetOUT, 
                                                      int mem_block_size, bool bparal, bool browmajor, 
                                                      Rcpp::Nullable<int> threads  = R_NilValue)
{
    

    int K = sizeA[0];
    int N = sizeA[1];
    
    int M = sizeB[0];
    int L = sizeB[1];
    
    
    //..// Rcpp::Rcout<<"Estem tractant unes matriu de tamany : \n\tMatriu A : "<<sizeA<<"\n\tMatriu B : "<<sizeB<<"\n\tMatriu C : "<<M<<" "<<N;
    
    IntegerVector stride = {1,1};
    IntegerVector block = {1,1};
    H5File* file  = nullptr;
    DataSet* datasetA = nullptr;
    DataSet* datasetB = nullptr;
    DataSet* datasetC = nullptr;
    
    try {
        
        if( K == M && N == L)
        {
            
            int isize = hdf5_block + 1;
            
            IntegerVector stride = IntegerVector::create(1, 1);
            IntegerVector block = IntegerVector::create(1, 1); 
            
            create_HDF5_dataset( filename, strsubgroupOUT + strdatasetOUT, K, N, "real");
            
            // Open file and get dataset
            file = new H5File( filename, H5F_ACC_RDWR );
            
            datasetA = new DataSet(file->openDataSet(strsubgroupIN + matA));
            datasetB = new DataSet(file->openDataSet(strsubgroupINB + matB));
            datasetC = new DataSet(file->openDataSet(strsubgroupOUT + strdatasetOUT));
            
            if( K<N ) {
                for (int ii = 0; ii < N; ii += hdf5_block)
                {
                    if( ii + hdf5_block > N ) isize = N - ii;
                    // Get blocks from hdf5 data file
                    Eigen::MatrixXd A = GetCurrentBlock_hdf5( file, datasetA, 0, ii, K, std::min(hdf5_block,isize));
                    Eigen::MatrixXd B = GetCurrentBlock_hdf5( file, datasetB, 0, ii, K, std::min(hdf5_block,isize));
                    Eigen::MatrixXd C = A + B;
                    
                    IntegerVector count = {K, std::min(hdf5_block,isize)};
                    IntegerVector offset = {0,ii};
                    
                    write_HDF5_matrix_subset_v2( file, datasetC, offset, count, stride, block, Rcpp::wrap(C));
                    
                    if( ii + hdf5_block > N ) isize = hdf5_block + 1;
                }
                
            } else {
                for (int ii = 0; ii < K; ii += hdf5_block)
                {
                    if( ii + hdf5_block > K ) isize = K - ii;
                    // Get blocks from hdf5 data file
                    Eigen::MatrixXd A = GetCurrentBlock_hdf5( file, datasetA, ii, 0, std::min(hdf5_block,isize), N);
                    Eigen::MatrixXd B = GetCurrentBlock_hdf5( file, datasetB, ii, 0, std::min(hdf5_block,isize), N);
                    Eigen::MatrixXd C = A + B;
                    
                    IntegerVector count = {std::min(hdf5_block,isize), N};
                    IntegerVector offset = {ii,0};
                    
                    write_HDF5_matrix_subset_v2( file, datasetC, offset, count, stride, block, Rcpp::wrap(C));
                    
                    if( ii + hdf5_block > K ) isize = hdf5_block + 1;
                }
                
            }
            
            
            
        }else {
            Rcpp::Rcout<<"non-conformable arguments\n";
            return void();
        }
        
        
    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        datasetA->close();
        datasetB->close();
        datasetC->close();
        file->close();
        ::Rf_error( "c++ exception hdf5_block_matrix_mul_hdf5_indatasets_transposed (File IException)" );
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        datasetA->close();
        datasetB->close();
        datasetC->close();
        file->close();
        ::Rf_error( "c++ exception hdf5_block_matrix_mul_hdf5_indatasets_transposed (DataSet IException)" );
        return void();
    } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        datasetA->close();
        datasetB->close();
        datasetC->close();
        file->close();
        ::Rf_error( "c++ exception hdf5_block_matrix_mul_hdf5_indatasets_transposed (DataType IException)" );
        return void();
    }catch(std::exception &ex) {
        datasetA->close();
        datasetB->close();
        datasetC->close();
        file->close();
        Rcpp::Rcout<< ex.what();
        return void();
    }
    
    datasetA->close();
    datasetB->close();
    datasetC->close();
    file->close();
    
    return void();
    
    

}



// [[Rcpp::export(.blockSum_hdf5)]]
Rcpp::RObject blockSum_hdf5(std::string filename, 
                             const std::string group, 
                             std::string A, 
                             std::string B,
                             Rcpp::Nullable<std::string> groupB = R_NilValue, 
                             Rcpp::Nullable<int> block_size = R_NilValue, 
                             Rcpp::Nullable<bool> paral = R_NilValue,
                             Rcpp::Nullable<int> threads = R_NilValue,
                             Rcpp::Nullable<std::string> outgroup = R_NilValue,
                             Rcpp::Nullable<std::string> outdataset = R_NilValue)
{
  
  int iblock_size, res;
  bool bparal, bexistgroup;// = Rcpp::as<double>;
  Eigen::MatrixXd C;
  
  H5File* file = nullptr;
  
  std::string strsubgroupOut, strdatasetOut;

  IntegerVector dsizeA, dsizeB;
  
  // hdf5 parameters
  try{
    
    H5::Exception::dontPrint();  
    
    if( outgroup.isNull()) {
      strsubgroupOut = "OUTPUT";
    } else {
      strsubgroupOut = Rcpp::as<std::string> (outgroup);
    }
    
    //..// std::string strsubgroup = "Base.matrices/";
    std::string strsubgroupIn = group + "/";
    std::string strsubgroupInB;
    std::string strGroupB;
    
    if(groupB.isNotNull()){
        strsubgroupInB =  Rcpp::as<std::string> (groupB) + "/";
        strGroupB = Rcpp::as<std::string> (groupB);
    } else {
        strsubgroupInB =  group + "/";
        strGroupB = group;
    }
    
    if( outdataset.isNotNull()) {
        strdatasetOut =  Rcpp::as<std::string> (outdataset);    
    } else {
        strdatasetOut =  A + "_x_" + B;
    }

    if(!ResFileExist(filename)) {
        Rcpp::Rcout<<"\nFile not exits, create file and datasets\n";  
        return wrap(-1);
    }
    // Open file and get dataset
    file = new H5File( filename, H5F_ACC_RDWR );
    
    
    if(exists_HDF5_element_ptr(file, group)==0) {
        Rcpp::Rcout<<"\nGroup not exits"<<group<<" \n";
        file->close();
        return Rcpp::wrap(-1);
    }  else{
        if(!exists_HDF5_element_ptr(file, group + "/" + A)) {
            Rcpp::Rcout<<"\n Dataset"<< A<<" not exits \n";
            file->close();
            return Rcpp::wrap(-1);
        }
    }
    
    if(exists_HDF5_element_ptr(file, strGroupB)==0) {
        Rcpp::Rcout<<"\nGroup not exits"<<strGroupB<<" \n";
        file->close();
        return Rcpp::wrap(-1);
    }  else{
        if(!exists_HDF5_element_ptr(file, strGroupB + "/" + B)) {
            Rcpp::Rcout<<"\n Dataset"<< B<<" not exits \n";
            file->close();
            return Rcpp::wrap(-1);
        }
    }
    
    
    DataSet dsA = file->openDataSet(strsubgroupIn + A);
    IntegerVector dsizeA = get_HDF5_dataset_size(dsA);
    DataSet dsB = file->openDataSet(strsubgroupInB + B);
    IntegerVector dsizeB = get_HDF5_dataset_size(dsB);

    bexistgroup = exists_HDF5_element_ptr(file,strsubgroupOut+ "/" );

    if(bexistgroup) {
        std::string strdataset = strsubgroupOut + "/" + strdatasetOut;
        if(exists_HDF5_element_ptr(file, strdataset )) {
            remove_HDF5_element_ptr(file, strdataset);
        }
    } 
    
    file->close();
    dsA.close();
    dsB.close();

    if(block_size.isNotNull())
    {
      iblock_size = Rcpp::as<int> (block_size);
      
    } else {
      iblock_size = std::min(  std::min(dsizeA[0],dsizeA[1]),  std::min(dsizeB[0],dsizeB[1]));
      if (iblock_size>512)
        iblock_size = 512;
    }

    if( paral.isNull()) {
      bparal = false;
    } else {
      bparal = Rcpp::as<bool> (paral);
    }
    

    if(!bexistgroup) {
      res = create_HDF5_group(filename, strsubgroupOut + "/" );
    } 
    

        
        // Test mix versión read block from file and calculate multiplication in memory (with paral·lel algorithm)
        hdf5_block_matrix_sum_hdf5_indatasets_transposed(A, B, dsizeA, dsizeB, iblock_size, filename, 
                                                         strsubgroupIn, strsubgroupInB, 
                                                         strsubgroupOut + "/", strdatasetOut, 
                                                         bparal, true, threads);


  } catch( FileIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception blockmult_hdf5 (File IException)" );
    return wrap(-1);
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    file->close();
    ::Rf_error( "c++ exception blockmult_hdf5 (DataSet IException)" );
    return wrap(-1);
  } catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }
  
  
  //..// return wrap(wrap(C));
  
  //..// return(C);
  return List::create(Named("filename") = filename,
                      Named("dataset") = strsubgroupOut + "/" + strdatasetOut,
                      Named("result") = wrap(0));
  
}


/***R
library(BigDataStatMeth)
library(rhdf5)

setwd("C:/tmp_test/")

# devtools::reload(pkgload::inst("BigDataStatMeth"))

# Prepare data and functions
X <- matrix(rnorm(200*10), nrow = 10, ncol = 200)
diag(X) <- 0.5

# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("test_file3.hdf5", X, "data", "X", force = T)
bdAdd_hdf5_matrix(X, "test_file3.hdf5", "data", "Y", force = T)

# Update diagonal
diagonal <- bdblockSum_hdf5( filename = "test_file3.hdf5", group = "data", a = "X", b = "Y", groupB = "data", block_size = 128, outgroup = "tmp", outdataset = "Z")

*/
