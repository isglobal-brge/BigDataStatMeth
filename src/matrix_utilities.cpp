#include "include/matrix_utilities.h"
#include "include/parallel_CrossProd.h"



Rcpp::NumericVector flatmatrm(Rcpp::NumericMatrix x)
{
  double xcol = x.ncol();
  double xrow = x.nrow();
  double dim = xrow * xcol;
  
  Rcpp::NumericVector fr( dim );
  
  for( double i=0; i< xrow; i++)
  {
    for( double j=0; j<xcol; j++)
    {
      fr[i*xcol+j] = x(i,j);
    }
  }
  
  return(fr);
}


Rcpp::NumericVector flatmatcm(Rcpp::NumericMatrix x)
{
  double xcol = x.ncol();
  double xrow = x.nrow();
  double dim = xrow * xcol;
  
  Rcpp::NumericVector fc( dim );
  
  for( double i=0; i< xrow; i++)
  {
    for( double j=0; j<xcol; j++)
    {
      fc[i+j*xrow] = x(i,j);
    }
  }
  return(fc);
}



Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X, bool bc, bool bs )
{
  Eigen::MatrixXd rX;
  
  if( bc==true && bs==true )  {
    
    Eigen::RowVectorXd mean = X.colwise().mean();
    Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
    rX = (X.rowwise() - mean).array().rowwise() / std.array();

  }   else if (bc == true  && bs==false)   {
    
    Eigen::RowVectorXd mean = X.colwise().mean();
    rX = (X.rowwise() - mean);
    
  }  else if ( bc == false && bs == true)   {
    
    Eigen::RowVectorXd mean = X.colwise().mean();
    Eigen::RowVectorXd std = (X.array().square().colwise().sum() / (X.rows() - 1)).sqrt();
    rX = X.array().rowwise() / std.array();
  } 

  return(rX);
}



// Get column means from Eigen::Matrix
Eigen::RowVectorXd get_cols_mean(Eigen::MatrixXd X)
{
  Eigen::VectorXd mean = X.rowwise().mean();
  return(mean);
}


// Get column std from Eigen::Matrix
Eigen::RowVectorXd get_cols_std(Eigen::MatrixXd X)
{
  Eigen::VectorXd mean = X.rowwise().mean();
  Eigen::VectorXd std = ((X.colwise() - mean).array().square().rowwise().sum() / (X.rows() - 1)).sqrt();
  return(std);
}






Rcpp::NumericMatrix RcppNormalize_Data_r ( Rcpp::NumericMatrix  x )
{
  Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd> (x);
  Eigen::RowVectorXd mean = X.colwise().mean();
  Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
  return Rcpp::wrap((X.rowwise() - mean).array().rowwise() / std.array());
}


//' Normalize Delayed Array matrix
//' 
//' This function performs a numerical or Delayed Array matrix normalization
//' 
//' @param X numerical or Delayed Array Matrix
//' @param bcenter logical (default = TRUE) if TRUE, centering is done by subtracting the column means
//' @param bscale logical (default = TRUE) if TRUE, centering is done by subtracting the column means
//' @return numerical matrix
//' @examples
//' library(DelayedArray)
//' 
//' m <- 500
//' n <- 100 
//' x <- matrix(rnorm(m*n), nrow=m, ncol=n)
//' 
//' # with numeric matrix
//' bdNormalize_Data(x)
//' 
//' # with Delaeyd Array
//' Dx <- DelayedArray(x)
//' 
//' # Center and scale
//' bdNormalize_Data(Dx)
//' 
//' # Only scale
//' bdNormalize_Data(Dx, bcenter = FALSE)
//' 
//' # Only center
//' bdNormalize_Data(Dx, bscale = FALSE)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdNormalize_Data ( Rcpp::RObject & X, 
                                 Rcpp::Nullable<bool> bcenter = R_NilValue,
                                 Rcpp::Nullable<bool> bscale  = R_NilValue)
{
  Eigen::MatrixXd mX;
  bool bc, bs;
  CharacterVector svrows, svrcols;
  
  
  List dimnames = X.attr( "dimnames" );
  
  if(dimnames.size()>0 ) {
    svrows= dimnames[0];
    svrcols = dimnames[1];
  }
  
  
  // Read DelayedArray's X and b
  if ( X.isS4() == true)    
  {
    mX = read_DelayedArray(X);
  } else {
    try{  
      if ( TYPEOF(X) == INTSXP ) {
        mX = Rcpp::as<Eigen::MatrixXi>(X).cast<double>()  ;
      } else{
        mX = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
      }
    }
    catch(std::exception &ex) { }
  }
  
  if( bcenter.isNull()) {
    bc = true;
  } else {
    bc = Rcpp::as<bool> (bcenter);
  }
  
  if( bscale.isNull()) {
    bs = true;
  } else {
    bs = Rcpp::as<bool> (bscale);
  }
  
  if( bc==true && bs==true )  {
    
    Eigen::RowVectorXd mean = mX.colwise().mean();
    Eigen::RowVectorXd std = ((mX.rowwise() - mean).array().square().colwise().sum() / (mX.rows() - 1)).sqrt();
    
    
    // Replace 0 to 0.00000000001 to avoid errors
    std = (std.array() == 0).select(0.00000000001, std);
    
    // Scale data
    mX = (mX.rowwise() - mean).array().rowwise() / std.array();
    
  }   else if (bc == true)   {
    
    Eigen::RowVectorXd mean = mX.colwise().mean();
    mX = (mX.rowwise() - mean);
    
  }  else if ( bs == true)   {
    
    Eigen::RowVectorXd std = (mX.array().square().colwise().sum() / (mX.rows() - 1)).sqrt();

    // Replace 0 to 0.00000000001 to avoid errors
    std = (std.array() == 0).select(0.00000000001, std);

    // Scale data
    mX = mX.array().rowwise() / std.array();
  } 

  NumericMatrix XR = wrap(mX);
  remove("mX");
  
  if(dimnames.size()>0 )  {
    rownames(XR) = as<CharacterVector>(dimnames[0]);
    colnames(XR) = as<CharacterVector>(dimnames[1]);  
  }
  
  
  return XR;
}


// Read Eigen matrix block size from X starting at (startrow,startcol) with size (nrows,ncols)
Eigen::MatrixXd GetCurrentBlock( Eigen::MatrixXd X, int startrow, int startcol, int nrows, int ncols )
{
  
  Eigen::MatrixXd m = X.block(startrow, startcol, nrows, ncols);;
  
  return(m);
  
}

