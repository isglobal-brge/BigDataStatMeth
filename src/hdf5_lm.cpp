#include "include/hdf5_lm.h"

using namespace RcppEigen;
using namespace Rcpp;


// Apply mlr_mr algorith to perform lm regression with big data
// Algorithm from : 
//    https://isglobal-brge.github.io/Master_Modelling/dealing-with-big-data-in-r.html#linear-regression-for-big-data
Eigen::MatrixXd Rcpp_mlr_mr_hdf5(Eigen::MatrixXd x, Eigen::MatrixXd y, int iblocks, Rcpp::Nullable<int> threads  = R_NilValue )
{
  
  int chunk = 1; // , tid;
  unsigned int ithreads;
  int irows =  x.rows(), 
    icols = x.cols();
  
  // Get number of threads
  if(threads.isNotNull()) 
  {
    if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
      ithreads = Rcpp::as<int> (threads);
    else 
      ithreads = std::thread::hardware_concurrency()/2;
  }
  else  ithreads = std::thread::hardware_concurrency()/2; //omp_get_max_threads()
  
  
  omp_set_num_threads(ithreads); 
  
  // Define variables
  Eigen::MatrixXd R1(iblocks*x.cols(), x.cols());
  Eigen::MatrixXd Q1( irows, icols );
  Eigen::VectorXd indexQ1( iblocks+1 );
  
#pragma omp parallel shared( R1, Q1, indexQ1, chunk) 
{
  // First steps --> (1) Read block
  //                 (2) Make decomposition by block and reduce result matrix
  
  double block_size = std::floor((double)irows/(double)iblocks);
  int maxsizetoread = block_size;
  int iSizetoRead; 
  
  
#pragma omp parallel for schedule(static)
  for( int i = 0; i<iblocks ; i++) 
  {
    strQR decQR;
    
    if( i+1 == iblocks && irows - maxsizetoread!=0)
      iSizetoRead = irows - (i*block_size);
    else
      iSizetoRead = maxsizetoread;
    
    // Read current block
    Eigen::MatrixXd obsBlock = GetCurrentBlock( x, i*block_size, 0, iSizetoRead, icols);
    
    // QR Decomposition
    decQR = rcpp_bdQR(obsBlock, true);
    
    // Complete R1
    R1.block( i*x.cols(), 0, decQR.R.rows(), decQR.R.cols() ) = decQR.R;
    
    if (i == 0) { indexQ1(i) = 0; } 
    indexQ1(i+1) = indexQ1(i) + decQR.Q.rows();
    
    // Complete Q1
    Q1.block(indexQ1(i), 0, decQR.Q.rows(), decQR.Q.cols() ) = decQR.Q;
    
  }
} 


// QR decomposition from R1
strQR decR1;
decR1 = rcpp_bdQR(R1, true);


Eigen::MatrixXd Q3;
Eigen::MatrixXd V(icols, iblocks);


double block_size = std::floor((double)decR1.Q.rows()/(double)iblocks);
int maxsizetoread = block_size;
int iSizetoRead; 
int startnext = 0;


#pragma omp parallel for ordered schedule(dynamic) 
for( int i = 0; i< iblocks ; i++) 
{
  
  Eigen::MatrixXd Q1BlockDiv;
  Eigen::MatrixXd Q2BlockDiv;
  
  if( i+1 == iblocks && decR1.Q.rows() - maxsizetoread!=0){
    iSizetoRead = decR1.Q.rows() - (i*block_size);
  }else{
    iSizetoRead = maxsizetoread;
  }
  
  // Get block from Q2    
  Q2BlockDiv = GetCurrentBlock( decR1.Q, i*block_size, 0, iSizetoRead, decR1.Q.cols());
  
  // Get block from Q1
  Q1BlockDiv = GetCurrentBlock( Q1, indexQ1(i), 0, ((indexQ1(i+1)) - indexQ1(i)), Q1.cols());
  
  // Get Q3 and V
  Q3 = Bblock_matrix_mul_parallel(Q1BlockDiv, Q2BlockDiv, 256, threads);
  
  Eigen::MatrixXd YBlock  = GetCurrentBlock( y, indexQ1(i), 0, Q3.rows(), 1);
  
  V.block(0, i, V.rows(), 1) = Bblock_matrix_mul_parallel( Q3.adjoint(), YBlock , 256, threads );
  
  startnext = startnext + Q3.adjoint().cols();
  
}


// Get Betas
Eigen::MatrixXd beta = Bblock_matrix_mul_parallel( decR1.R.inverse(), V.rowwise().sum(), 256, threads);

return(beta);
}



// // ' Linear regression using MLR-MR algorithm
// // '
// // ' Linear regression for Big Data using MLR-MR algorithm
// // '
// // ' @param X, numerical matrix with paired observations of the predictor variable X
// // ' @param Y, numerical matrix column with response variable
// // ' @param outgroup, character array indicating group where the data set will be saved after remove data with if `outgroup` is NULL, output dataset is stored in the same input group. 
// // ' @param outdataset, character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
// // ' @param blocks, integer with number of blocks we want to split matrix if null matrix is splited in blocks as maximum of 1000 variables per block
// // ' @param threads, threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
// // ' @return Lineal regression coefficients
// //' @export
// [[Rcpp::export(.bdMLR_MR_hdf5)]]
Rcpp::RObject bdMLR_MR_hdf5(std::string filename, 
                            const std::string group, 
                            std::string dataset,
                            const std::string betasgroup, 
                            std::string betasdataset,
                            int blocks, 
                            Rcpp::Nullable<std::string> outgroup = R_NilValue,
                            Rcpp::Nullable<std::string> outdataset = R_NilValue,
                            Rcpp::Nullable<int> threads  = R_NilValue) 
{
  
  try
  {
    
    H5::Exception::dontPrint();  
    std::string strOutgroup = "OUTPUT";
    std::string strOutdataset= "lm_" + dataset;
    
    
    // hdf5 variables
    H5File* file;
    DataSet* datasetOut;
    
    if( outgroup.isNotNull() ) {
      strOutgroup = as<std::string>(outgroup);
    }
    
    if( outdataset.isNotNull() ) {
      strOutdataset = as<std::string>(outdataset);
    }
    
    //.. Aquí faltaria definir on escriurà els resultats ..// 
    std::string stroutdata = strOutgroup +"/" + strOutdataset;
    std::string strdataset = group +"/" + dataset;
    std::string strdatasetbetas = betasgroup +"/" + betasdataset;
    
    
    if(!ResFileExist_filestream(filename)){
      throw std::range_error("File not exits, create file before access to dataset");
    }
    
    
    file = new H5File( filename, H5F_ACC_RDWR );
    
    if(exists_HDF5_element_ptr(file, strdataset)) 
    {
      
      
      
      /***
      
      
      // Aquí el que s'ha de fer es acabar de comprobar que tot estigui bé i desprès
      // passar els punters dels datasets a la funció per calcular el lm
      
      
      ***/
      
      
    }else{
      //.commented 20201120 - warning check().// pdataset->close();
      file->close();
      throw std::range_error("Dataset does not exits");  
    }
   
   
    
    
    
    //.  CODI A EXECUTAR .// Eigen::MatrixXd beta = Rcpp_mlr_mr(eX, eY, blocks, threads);
    
    //.. Això és el que hauria de ser : ..// return (wrap(beta));
    return (0);
    
  } catch( FileIException& error ){ // catch failure caused by the H5File operations
    //..// file->close();
    ::Rf_error( "c++ exception bdMLR_MR_hdf5 (File IException)" );
    return(wrap(-1));
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    //..// file->close();
    ::Rf_error( "c++ exception bdMLR_MR_hdf5 (DataSet IException)" );
    return(wrap(-1));
  } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
    //..// file->close();
    ::Rf_error( "c++ exception bdMLR_MR_hdf5 (DataSpace IException)" );
    return(wrap(-1));
  } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
    //..// file->close();
    ::Rf_error( "c++ exception bdMLR_MR_hdf5 (DataType IException)" );
    return(wrap(-1));
  }catch(std::exception &ex) {
    //..// file->close();
    Rcpp::Rcout<< ex.what();
    return(wrap(-1));
  }
  
  // other cases
  return (wrap(-1));
  
}


/*** R

library(BigDataStatMeth)
data(mtcars)

Y <- mtcars$mpg
X <- model.matrix(~ wt + cyl, data=mtcars)
m <- 7


res <- bdMLR_MR( X, Y, m, 1)
res

*/
