#include "include/hdf5_blockmult.h"


//' Multiply hdf5 matrix
//' 
//' This function multiply matrix stored in hdf5 data file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string or Delayed Array Matrix
//' @param dataset string or Delayed Array Matrix
//' @param bcenter logical (default = TRUE) if TRUE, centering is done by subtracting the column means
//' @param bscale logical (default = TRUE) if TRUE, centering is done by subtracting the column means
//' @param wsize integer (default = 1000), file block size to read to perform normalization
//' @return file with scaled, centered or scaled and centered dataset
//' @examples
//'   a = "See vignette"
//' @export
// [[Rcpp::export]]
Rcpp::RObject blockmult_hdf5(std::string filename, const std::string group, 
                             std::string A, std::string B,
                             Rcpp::Nullable<int> block_size = R_NilValue, 
                             Rcpp::Nullable<bool> paral = R_NilValue,
                             Rcpp::Nullable<int> threads = R_NilValue,
                             Rcpp::Nullable<double> mixblock_size = R_NilValue,
                             Rcpp::Nullable<std::string> outgroup = R_NilValue)
{
  
  int iblock_size, res;
  bool bparal, bexistgroup;// = Rcpp::as<double>;
  Eigen::MatrixXd C;
  
  
  H5File* file;
  
  
  std::string strsubgroupOut;

  
  IntegerVector dsizeA, dsizeB;
  
  // hdf5 parameters
  try{
    
    H5::Exception::dontPrint();  
    
    if( outgroup.isNull()) {
      strsubgroupOut = "OUTPUT/";
    } else {
      strsubgroupOut = Rcpp::as<std::string> (outgroup) + "/";
    }
    
    //..// std::string strsubgroup = "Base.matrices/";
    std::string strsubgroupIn = group + "/";

    // Open file and get dataset
    file = new H5File( filename, H5F_ACC_RDONLY );
    
    DataSet dsA = file->openDataSet(strsubgroupIn + A);
    IntegerVector dsizeA = get_HDF5_dataset_size(dsA);
    DataSet dsB = file->openDataSet(strsubgroupIn + B);
    IntegerVector dsizeB = get_HDF5_dataset_size(dsB);
    
    bexistgroup = exists_HDF5_element_ptr(file,strsubgroupOut );

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
      hdf5_block_matrix_mul_hdf5_indatasets_transposed(A, B, dsizeA, dsizeB, iblock_size, filename, strsubgroupIn, strsubgroupOut, 
                                             memory_block, bparal,true, threads);
      
      
      //..// C = Eigen::MatrixXd::Zero(2,2);
      
    }else if (bparal == false)
    {
      
      // Not parallel
      hdf5_block_matrix_mul_hdf5_indatasets_transposed(A, B, dsizeA, dsizeB, iblock_size, filename, strsubgroupIn, strsubgroupOut, 
                                             0, bparal,true, threads);
      //..// C = Eigen::MatrixXd::Zero(2,2);
      
    }

  } catch( FileIException error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception (File IException)" );
    return wrap(-1);
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    file->close();
    ::Rf_error( "c++ exception (DataSet IException)" );
    return wrap(-1);
  } catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }
  
  
  //..// return wrap(wrap(C));
  
  //..// return(C);
  return List::create(Named("filename") = filename,
                      Named("dataset") = strsubgroupOut + "/" + A + "_x_" + B,
                      Named("result") = wrap(0));
  
}


/***R
library(BigDataStatMeth)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/BigDataStatMeth")
devtools::document()


library(BigDataStatMeth)
library(DelayedArray)
library(rhdf5)

setwd( "/Users/mailos/BitDataStatMeth_test/BasicMatVect" )

n <- 7
m <- 5
p <- 5
q <- 12

# R Object

set.seed(222)
matA <- matrix(runif(n*m), nrow = n, ncol = m)
set.seed(333)
matB <- matrix(runif(p*q), nrow = p, ncol = q)


matAD <- DelayedArray(matA)
matBD <- DelayedArray(matB)


Create_HDF5_matrix_file("BasicMatVect.hdf5", matA, "INPUT", "matAA")
#..# Create_HDF5_matrix_file("BasicMatVect.hdf5", matAD, "INPUT", "matA")
Create_HDF5_matrix(matBD, "BasicMatVect.hdf5", "INPUT", "matB")
Create_HDF5_matrix(matAD, "BasicMatVect.hdf5", "INPUT", "matA")
#..# Create_HDF5_matrix(matA, "BasicMatVect.hdf5", "INPUT", "matAA")
Create_HDF5_matrix(matB, "BasicMatVect.hdf5", "INPUT", "matBB")

BigDataStatMeth::blockmult_hdf5("BasicMatVect.hdf5", "INPUT","matA","matB")
BigDataStatMeth::blockmult_hdf5("BasicMatVect.hdf5", "INPUT","matAA","matBB")



# Examine hierarchy before open file
h5ls("BasicMatVect.hdf5")

# Open file
h5fdelay = H5Fopen("BasicMatVect.hdf5")
# Show hdf5 hierarchy (groups)
h5fdelay

mDA = h5fdelay$INPUT$matA
mDB = h5fdelay$INPUT$matB

mA = h5fdelay$INPUT$matAA
mB = h5fdelay$INPUT$matBB


matA
mA
mDA

matB
mB
mDB

# Close delayed.hdf5 file
H5Fclose(h5fdelay)





matAxmatB


matC <- matA[1:5,1:5]%*%matB[1:5,1:5]
matC

*/