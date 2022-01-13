#include "include/hdf5_blockmultSparse.h"


//' Block matrix multiplication
//' 
//' This function performs a block matrix-matrix multiplication with numeric matrix
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string Matrix
//' @param A, string with dataset name where matrix is stored
//' @param B, string with dataset name where matrix is stored
//' @param outgroup string with de group name under the matrix will be stored
//' @return list with filename and the group and dataset name under the results are stored
//' @examples
//' 
//' library(Matrix)
//' library(BigDataStatMeth)
//' 
//' k <- 1e3
//' set.seed(1)
//' x_sparse <- sparseMatrix(
//'     i = sample(x = k, size = k),
//'     j = sample(x = k, size = k),
//'     x = rnorm(n = k)
//' )
//' set.seed(2)
//' y_sparse <- sparseMatrix(
//'     i = sample(x = k, size = k),
//'     j = sample(x = k, size = k),
//'     x = rnorm(n = k)
//' )
//' 
//' if( isTRUE(file.exists('BasicMatVect.hdf5'))) {
//'      file.remove('BasicMatVect.hdf5')
//' }
//' bdCreate_hdf5_matrix_file("BasicMatVect.hdf5", as.matrix(x_sparse), "SPARSE", "x_sparse")
//' bdAdd_hdf5_matrix(as.matrix(y_sparse), "BasicMatVect.hdf5", "SPARSE", "y_sparse")
//' 
//' d <- bdblockmult_sparse_hdf5("BasicMatVect.hdf5", "SPARSE", "x_sparse", "y_sparse")
//' 
//' # Remove file (used as example)
//' if (file.exists("BasicMatVect.hdf5")) {
//'   # Delete file if it exist
//'   file.remove("BasicMatVect.hdf5")
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdblockmult_sparse_hdf5(std::string filename, const std::string group, 
                             std::string A, std::string B,
                             Rcpp::Nullable<std::string> outgroup = R_NilValue )
{
   
   
   std::string strsubgroupOut, strdataset;
   std::string spMatrix("dgCMatrix");
   
   bool bexistgroup;
   bool bsparseA, bsparseB;
   int res;
   
   IntegerVector dsizeA, dsizeB;
   
   H5File* file;
   DataSet* dsA; 
   DataSet* dsB;
   
   
   
   try{
      
      H5::Exception::dontPrint();  
      
      if( outgroup.isNull()) {
         strsubgroupOut = "OUTPUT";
      } else {
         strsubgroupOut = Rcpp::as<std::string> (outgroup);
      }
      
      //..// std::string strsubgroup = "Base.matrices/";
      std::string strsubgroupIn = group + "/";

      // Open file and get dataset
      file = new H5File( filename, H5F_ACC_RDWR );

      dsA = new DataSet(file->openDataSet(strsubgroupIn + A));
      Rcpp::IntegerVector dsizeA = get_HDF5_dataset_size(*dsA);

      dsB = new DataSet(file->openDataSet(strsubgroupIn + B));
      Rcpp::IntegerVector dsizeB = get_HDF5_dataset_size(*dsB);
      
     
      bexistgroup = exists_HDF5_element_ptr(file,strsubgroupOut+ "/" );
      strdataset = strsubgroupOut+ "/" + A + "_x_" + B;
      
      if(bexistgroup) {
         if(exists_HDF5_element_ptr(file, strdataset )) {
            remove_HDF5_element_ptr(file, strdataset);
         }
      }
      
      
     

      Eigen::MatrixXd A = GetCurrentBlock_hdf5_Original(file, dsA, 0, 0, dsizeA[0], dsizeA[1] );
      Eigen::SparseMatrix<double> A_sp(dsizeA[0], dsizeA[1]);
      bsparseA = is_sparse(A);
      
      if(bsparseA) {
         A_sp = A.sparseView();
         A.resize(0,0); }
      
      Eigen::MatrixXd B = GetCurrentBlock_hdf5_Original(file, dsB, 0, 0, dsizeB[0], dsizeB[1] );
      Eigen::SparseMatrix<double> B_sp(dsizeB[0], dsizeB[1]);
      
      bsparseB = is_sparse(B);
      
      if(bsparseB) {
         B_sp = B.sparseView();
         B.resize(0,0); }
      
      dsA->close();
      dsB->close();
      
      
      Eigen::SparseMatrix<double> C_sp(dsizeA[0], dsizeB[1]);
      
      if(bsparseA || bsparseB) 
      {
         
         if(!bsparseA){
            Rcpp::Rcout<<"Matrix A isn't a sparse matrix";
         } 
         
         if(!bsparseB){
            Rcpp::Rcout<<"Matrix B isn't a sparse matrix";
         }
         
         res = create_HDF5_group_ptr(file, strsubgroupOut );
         
         C_sp = A_sp * B_sp;

         write_HDF5_matrix_transposed_ptr(file, strdataset, wrap( Eigen::MatrixXd(C_sp)) );
         
      } else {
         
         Rf_error("No sparse matrix found");
         file->close();
         return(wrap(-1));
      }
      
      
      file->close();
      return(wrap(C_sp));
      
   } catch( FileIException& error ) { // catch failure caused by the H5File operations
      file->close();
      ::Rf_error( "c++ exception bdblockmult_sparse_hdf5 (File IException)" );
      return wrap(-1);
   } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
      file->close();
      ::Rf_error( "c++ exception bdblockmult_sparse_hdf5 (DataSet IException)" );
      return wrap(-1);   
   } catch(std::exception &ex) {
      Rcpp::Rcout<< ex.what();
      return wrap(-1);
   }
   
   
}



