#include "include/hdf5_blocktCrossprod.h"



// Working directly with C-matrix in hdf5 file
// If option paral·lel is enabled, this function loads hdf5 read blocks
// of medium size into memory and calculates the multiplication of blocks
// by applying the parallel algorithm Bblock_matrix_mul_parallel (in-memory process)
// browmajor : if = true, indicates that R data is stored in hdf5 as row major (default in hdf5)
//             else, indicates that R data is stored in hdf5 as column major
int hdf5_block_matrix_tcrossprod_hdf5( std::string matA, IntegerVector sizeA, 
                                       std::string matB, IntegerVector sizeB, 
                                       int hdf5_block,
                                       std::string filename, 
                                       std::string strsubgroupIN,std::string strsubgroupINB, 
                                       std::string strsubgroupOUT, 
                                       int mem_block_size, bool bparal, bool browmajor, 
                                       Rcpp::Nullable<int> threads  = R_NilValue)
{
  
  
  int N = sizeA[1];
  int K = sizeA[0];
  int M = sizeB[1];
  int L = sizeB[0];
  
  IntegerVector stride = {1,1};
  IntegerVector block = {1,1};
  
  H5File* file = nullptr;
  DataSet* datasetA = nullptr;
  DataSet* datasetB = nullptr;
  DataSet* datasetC = nullptr;
  
  if( K == L)
  {
      
    int isize = hdf5_block + 1;
    int ksize = hdf5_block + 1;
    int jsize = hdf5_block + 1;
    
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1); 

    // Create an empty dataset for C-matrix into hdf5 file
    //.Works OK.// create_HDF5_dataset( filename, strsubgroupOUT + "tCrossProd_" + matA + "x" + matB, N, M, "real");
    create_HDF5_dataset( filename, strsubgroupOUT + "tCrossProd_" + matA + "x" + matB, M, N, "real");

    // Open file and get dataset
    file = new H5File( filename, H5F_ACC_RDWR );
    
    datasetA = new DataSet(file->openDataSet(strsubgroupIN + matA));
    datasetB = new DataSet(file->openDataSet(strsubgroupINB + matB));
    datasetC = new DataSet(file->openDataSet(strsubgroupOUT + "tCrossProd_" + matA + "x" + matB));
    
    
    for (int ii = 0; ii < N; ii += hdf5_block)
    {
      
      if( ii + hdf5_block > N ) isize = N - ii;
      for (int jj = 0; jj < M; jj += hdf5_block)
      {
        
        if( jj + hdf5_block > M) jsize = M - jj;
        for(int kk = 0; kk < K; kk += hdf5_block)
        {
          if( kk + hdf5_block > K ) ksize = K - kk;
          
          // Get blocks from hdf5 file
          Eigen::MatrixXd A = GetCurrentBlock_hdf5_Original( file, datasetA, kk, ii, 
                                                    std::min(hdf5_block,ksize),std::min(hdf5_block,isize));
          
          Eigen::MatrixXd B = GetCurrentBlock_hdf5( file, datasetB, kk, jj, 
                                                    std::min(hdf5_block,ksize),std::min(hdf5_block,jsize));

          //.Works OK.// Eigen::MatrixXd C = GetCurrentBlock_hdf5( file, datasetC, ii, jj, 
          //.Works OK.//                                           std::min(hdf5_block,isize),std::min(hdf5_block,jsize));
          
          Eigen::MatrixXd C = GetCurrentBlock_hdf5( file, datasetC, jj, ii, 
                                                    std::min(hdf5_block,jsize),std::min(hdf5_block,isize));
          C.transposeInPlace();

          if( bparal == false)
            C = C + A*B;
          else
            C = C + Bblock_matrix_mul_parallel(A, B, mem_block_size, threads);

          //.Works OK.// IntegerVector count = {std::min(hdf5_block,isize), std::min(hdf5_block,jsize)};
          //.Works OK.// IntegerVector offset = {ii,jj};
          
          IntegerVector count = {std::min(hdf5_block,jsize), std::min(hdf5_block,isize)};
          IntegerVector offset = {jj,ii};
          
          write_HDF5_matrix_subset_v2( file, datasetC, offset, count, stride, block, Rcpp::wrap(C.transpose()));
          
          if( kk + hdf5_block > K ) ksize = hdf5_block + 1;
        }
        
        if( jj + hdf5_block > M ) jsize = hdf5_block + 1;
      }
      
      if( ii + hdf5_block > N ) isize = hdf5_block + 1;
    }
    
    
  }else {
    throw std::range_error("non-conformable arguments");
  }
  
  datasetA->close();
  datasetB->close();
  datasetC->close();
  file->close();
  
  delete(datasetA);
  delete(datasetC);
  delete(datasetB);
  delete(file);
  
  return(0);
}







//' Transposed Crossprod with hdf5 matrix
//' 
//' This function performs the transposed crossprod from a matrix inside and hdf5 data file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group, string, group name where dataset A is stored
//' @param A string name inside HDF5 file
//' @param groupB, string, group name where dataset b is stored
//' @param B string, dataset name for matrix B inside HDF5 file
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @param mixblock_size (optional) only for debug pourpose
//' @param outgroup (optional) group name to store results from Crossprod inside hdf5 data file
//' @examples
//' 
//' library(BigDataStatMeth)
//' library(rhdf5)
//' 
//' matA <- matA <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), nrow = 3, byrow = TRUE)
//' matB <- matrix(c(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,5,3,4,5,2,6,2,3,4,
//'                  42, 23, 23, 423,1,2), ncol = 5, byrow = TRUE)
//' 
//' bdCreate_hdf5_matrix_file("BasicMatVect.hdf5", matA, "INPUT", "matA")
//' bdAdd_hdf5_matrix( matB, "BasicMatVect.hdf5", "INPUT", "matB")
//' 
//' res <- bdtCrossprod_hdf5("BasicMatVect.hdf5", "INPUT",
//'                                         "matA", block_size = 2)
//' res2 <- bdtCrossprod_hdf5("BasicMatVect.hdf5", "INPUT", 
//'                                         "matA", "INPUT","matB", block_size = 2)
//' 
//' # Open file
//' h5fdelay = H5Fopen("BasicMatVect.hdf5")
//' 
//' res <- h5fdelay$OUTPUT$tCrossProd_matAxmatA
//' res2 <- h5fdelay$OUTPUT$tCrossProd_matAxmatB
//' 
//' all.equal(tcrossprod(matA), res)
//' all.equal(tcrossprod(matA, matB), res2)
//' 
//' # Close delayed.hdf5 file
//' H5Fclose(h5fdelay)
//' 
//' # Remove file (used as example)
//' if (file.exists("BasicMatVect.hdf5")) {
//'   # Delete file if it exist
//'   file.remove("BasicMatVect.hdf5")
//' }
//' 
//' @return If all process is ok, returns a list with : 
//' \itemize{
//'   \item{"filename"}{ File name where results are storesd }
//'   \item{"dataset"}{ route to results inside hdf5 data file }
//' }
//' 
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdtCrossprod_hdf5(std::string filename, const std::string group, 
                                std::string A,
                                Rcpp::Nullable<std::string> groupB = R_NilValue, 
                                Rcpp::Nullable<std::string>  B = R_NilValue, 
                                Rcpp::Nullable<int> block_size = R_NilValue, 
                                Rcpp::Nullable<bool> paral = R_NilValue,
                                Rcpp::Nullable<int> threads = R_NilValue,
                                Rcpp::Nullable<double> mixblock_size = R_NilValue,
                                Rcpp::Nullable<std::string> outgroup = R_NilValue)
{


  int iblock_size, res;
  bool bparal, bexistgroup;// = Rcpp::as<double>;
  Eigen::MatrixXd C;
  std::string matB;
  std::string strsubgroupIn, strsubgroupInB;
  
  
  H5File* file = nullptr;
  
  std::string strsubgroupOut;

  IntegerVector dsizeA;
  
  // hdf5 parameters
  try {
    
    H5::Exception::dontPrint();  
    
    if( outgroup.isNull()) {
      strsubgroupOut = "OUTPUT/";
    } else {
      strsubgroupOut = Rcpp::as<std::string> (outgroup) + "/";
    }
    
    //..// std::string strsubgroup = "Base.matrices/";
    strsubgroupIn = group + "/";
    
    

    // Open file and get dataset
    file = new H5File( filename, H5F_ACC_RDONLY );
    
    DataSet dsA = file->openDataSet(strsubgroupIn + A);
    IntegerVector dsizeA = get_HDF5_dataset_size(dsA);
    bexistgroup = exists_HDF5_element_ptr(file,strsubgroupOut );
    
    dsA.close();
    file->close();
    delete(file);
    
    DataSet dsB;
    IntegerVector dsizeB;
    
    
    if( groupB.isNull()) {
      strsubgroupInB = group + "/";
      dsizeB = dsizeA;
      matB =  A;
      
    } else {
      
      strsubgroupInB = as<std::string>(groupB) + "/";
      
      // std::string strgroup = as<std::string>(groupB) + "/";
      matB = as<std::string>(B);
      
      // Open file and get dataset
      file = new H5File( filename, H5F_ACC_RDONLY );
      
      DataSet dsB = file->openDataSet(strsubgroupInB + matB);
      dsizeB = get_HDF5_dataset_size(dsB);
      
      dsB.close();
      file->close();
      delete(file);
    }
    
    if(block_size.isNotNull())
    {
      iblock_size = Rcpp::as<int> (block_size);
    } else {
      iblock_size =  std::min(std::min(dsizeA[0], dsizeA[1]), std::min(dsizeB[0], dsizeB[1]));
      if (iblock_size>1024) {
        iblock_size = 1024;
      }
    }
    
    if( paral.isNull()) {
      bparal = false;
    } else {
      bparal = Rcpp::as<bool> (paral);
    }
    
    if(!bexistgroup) {
      res = create_HDF5_group(filename, strsubgroupOut );
    }

    if(bparal == true)
    {
      //.. TODO : Work with parallel hdf5 access
      //..// C = hdf5_block_matrix_mul_parallel( dsizeA, dsizeB, iblock_size, filename, strsubgroup, threads );
      
      int memory_block; // Block size to apply to read hdf5 data to paralelize calculus
      
      if(mixblock_size.isNotNull())
        memory_block = Rcpp::as<int> (mixblock_size);
      else 
        memory_block = 128;
      
      // Test mix versión read block from file and calculate multiplication in memory (with paral·lel algorithm)
      hdf5_block_matrix_tcrossprod_hdf5(A, dsizeA, matB, dsizeB, iblock_size, filename, strsubgroupIn, strsubgroupInB, strsubgroupOut, 
                                             memory_block, bparal,true, threads);
      
    }else if (bparal == false)
    {
      
      // Not parallel
      hdf5_block_matrix_tcrossprod_hdf5(A, dsizeA, matB, dsizeB, iblock_size, filename, strsubgroupIn, strsubgroupInB, strsubgroupOut, 
                                       0, bparal,true, threads);
      
    }
    

  } catch( FileIException& error ) { // catch failure caused by the H5File operations
      file->close();
      delete(file); 
      ::Rf_error( "c++ exception tCrossprod_hdf5 (File IException)" );
      return wrap(-1);
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
      file->close();
      delete(file);
      ::Rf_error( "c++ exception tCrossprod_hdf5 (DataSet IException)" );
      return wrap(-1);
  } catch(std::exception &ex) {
      file->close();
      delete(file);
      Rcpp::Rcout<< ex.what();
      return wrap(-1);
  }
  

  std::string strname;
  
  if(B.isNull()) {
    strname = "tCrossProd_" + A;
  } else {
    strname = "tCrossProd_" + A + "x" + as<std::string>(B);
  }
  
  return List::create(Named("filename") = filename,
                      Named("dataset") = strsubgroupOut + "/" +  strname);
  
}

