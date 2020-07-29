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

/*
Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X )
{
  Eigen::RowVectorXd mean = X.colwise().mean();
  Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
  return (X.rowwise() - mean).array().rowwise() / std.array();
}
 */


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



Eigen::MatrixXd RcppNormalize_Data_hdf5 ( Eigen::MatrixXd  X, bool bc, bool bs, bool btransp, Eigen::MatrixXd normdata )
{
  Eigen::MatrixXd rX;
  
  if( btransp == true)
  {
    if( bc==true && bs==true )  {
      
      rX = (X.colwise() - normdata.row(0).transpose() ).array().colwise() / normdata.row(1).transpose().array();

    }   else if (bc == true  && bs==false)   {
      
      Eigen::VectorXd mean = X.rowwise().mean();
      rX = (X.colwise() - normdata.row(0).transpose());
      
    }  else if ( bc == false && bs == true)   {
      rX = X.array().colwise() / normdata.row(1).transpose().array();
    }
    
  } else {
  
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
//' m <- 500
//' n <- 100 
//' x <- matrix(rnorm(m*n), nrow=m, ncol=n)
//' 
//' # with numeric matrix
//' Normalize_Data(x)
//' 
//' # with Delaeyd Array
//' Dx <- DelayedArray(x)
//' 
//' # Center and scale
//' Normalize_Data(Dx)
//' 
//' # Only scale
//' Normalize_Data(Dx, bcenter = FALSE)
//' 
//' # Only center
//' Normalize_Data(Dx, bscale = FALSE)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject Normalize_Data ( Rcpp::RObject & x, 
                               Rcpp::Nullable<bool> bcenter = R_NilValue,
                               Rcpp::Nullable<bool> bscale  = R_NilValue)
{
  Eigen::MatrixXd X;
  bool bc, bs;
  CharacterVector svrows, svrcols;
  
  List dimnames = x.attr( "dimnames" );
  
  if(dimnames.size()>0 ) {
    svrows= dimnames[0];
    svrcols = dimnames[1];
  }
  
  
  // Read DelayedArray's x and b
  if ( x.isS4() == true)    
  {
    X = read_DelayedArray(x);
  } else {
    try{  
      if ( TYPEOF(x) == INTSXP ) {
        X = Rcpp::as<Eigen::MatrixXi>(x).cast<double>()  ;
      } else{
        X = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(x);
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
    
    Eigen::RowVectorXd mean = X.colwise().mean();
    Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
    X = (X.rowwise() - mean).array().rowwise() / std.array();
    
  }   else if (bc == true)   {
    
    Eigen::RowVectorXd mean = X.colwise().mean();
    X = (X.rowwise() - mean);
    
  }  else if ( bs == true)   {
    
    Eigen::RowVectorXd mean = X.colwise().mean();
    Eigen::RowVectorXd std = (X.array().square().colwise().sum() / (X.rows() - 1)).sqrt();
    X = X.array().rowwise() / std.array();
  } 

  NumericMatrix XR = wrap(X);
  remove("X");
  
  if(dimnames.size()>0 )  {
    rownames(XR) = as<CharacterVector>(dimnames[0]);
    colnames(XR) = as<CharacterVector>(dimnames[1]);  
  }
  
  
  return XR;
}



/***R
library(BigDataStatMeth)
m <- 10
n <- 10 
x <- matrix(rnorm(m*n), nrow=m, ncol=n)

Normalize_Data(x, bscale = TRUE, bcenter = FALSE )

scale(x, scale = TRUE, center = FALSE)

Dx <- DelayedArray(x)
Normalize_Data(Dx)


*/

