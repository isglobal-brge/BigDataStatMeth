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


/***
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
***/

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


/***
//' Normalize dataset in hdf5 file
//' 
//' This function normalize data scaling, centering or scaling and centering in a dataset stored in hdf5 file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string or Delayed Array Matrix
//' @param dataset string or Delayed Array Matrix
//' @param bcenter logical (default = TRUE) if TRUE, centering is done by subtracting the column means
//' @param bscale logical (default = TRUE) if TRUE, centering is done by subtracting the column means
//' @param wsize integer (default = 1000), file block size to read to perform normalization
//' @return file with scaled, centered or scaled and centered dataset
//' @examples
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject Normalize_hdf5(std::string filename, const std::string group, std::string dataset,
                             Rcpp::Nullable<bool> bcenter = R_NilValue, Rcpp::Nullable<bool> bscale  = R_NilValue,
                             Rcpp::Nullable<int> wsize  = R_NilValue)
{
  
  bool bc, bs;
  int blocksize;
  std::string strgroupout;
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  
  H5File* file;
    
  DataSet* pdatasetin;
  DataSet* pdatasetout;


  try{
    if( bcenter.isNull())
      bc = true;
    else
      bc = Rcpp::as<bool> (bcenter);
    
    if( bscale.isNull())
      bs = true;
    else
      bs = Rcpp::as<bool> (bscale);
    
    if( wsize.isNull())
      blocksize = 1000;
    else
      blocksize = Rcpp::as<int> (wsize);
    

    if(!ResFileExist(filename)) {
      Rcpp::Rcout<<"\nFile not exits, create file before normalize dataset\n";  
      return wrap(-1);
    }
    file = new H5File( filename, H5F_ACC_RDWR );


    if(exists_HDF5_element_ptr(file, group)==0) {
      Rcpp::Rcout<<"\nGroup not exits, create file and dataset before normalize data\n";
      file->close();
      return wrap(-1);
    }  else{
      if(!exists_HDF5_element_ptr(file, group + "/" + dataset)) {
        Rcpp::Rcout<<"\n Dataset not exits, create file and dataset before normalize data \n";
        file->close();
        return wrap(-1);
      }
    }
    
    pdatasetin = new DataSet(file->openDataSet(group + "/" + dataset));

    IntegerVector dims_out = get_HDF5_dataset_size(*pdatasetin);
    
    Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2,dims_out[0]);

    // Get data to normalize matrix
    get_HDF5_mean_sd_by_column_ptr( file, pdatasetin, datanormal);

    
    Rcpp::Rcout<<"\n Mitjanes i sd : \n"<<datanormal<<"\n";
    
    for(hsize_t i=0; i*blocksize <= dims_out[0] ; i++)
    {
      int sizetoread = 0;
      if((i+1)*blocksize<dims_out[0])
        sizetoread = blocksize;
      else
        sizetoread = dims_out[0]-(i*blocksize);
      
      // Prepare file and dataset
      if(i==0)
      {

        strgroupout = "NORMALIZED/" + group;
        
        // Mirar si existeix el grup NORMALIZED al fitxer --> Crear-lo
        if(exists_HDF5_element_ptr(file, strgroupout))
        {
          if(exists_HDF5_element_ptr(file, strgroupout + "/" + dataset))
            remove_HDF5_element_ptr(file, strgroupout + "/" + dataset);
        } else {
          create_HDF5_groups_ptr( file, strgroupout);
        }

        // Create dataset
        create_HDF5_dataset_ptr(file, strgroupout + "/" + dataset, dims_out[0], dims_out[1], "real");
        pdatasetout = new DataSet(file->openDataSet(strgroupout + "/" + dataset));
      
      }
      
      IntegerVector offset = IntegerVector::create( i*blocksize, 0);
      IntegerVector count = IntegerVector::create(sizetoread, dims_out[1] );
      
      // Normalize and write data
      //..// Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, pdatasetin, offset[0], offset[1], count[0], count[1]);
      Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, pdatasetin, offset[0], offset[1], count[0], count[1]);

      X = RcppNormalize_Data_hdf5(X, bc, bs, true, datanormal.block(0,offset[0], 2, count[0]));
      
      write_HDF5_matrix_subset_v2(file, pdatasetout, offset, count, block, stride, wrap(X));
      
    }
  } catch( FileIException error ) { // catch failure caused by the H5File operations
    file->close();
    pdatasetin->close();
    pdatasetout->close();
    error.printErrorStack();
    return wrap(-1);
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    file->close();
    pdatasetin->close();
    pdatasetout->close();
    error.printErrorStack();
    return wrap(-1);
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    file->close();
    pdatasetin->close();
    pdatasetout->close();
    error.printErrorStack();
    return wrap(-1);
  } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
    file->close();
    pdatasetin->close();
    pdatasetout->close();
    error.printErrorStack();
    return wrap(-1);
  }catch(std::exception &ex) {
    file->close();
    pdatasetin->close();
    pdatasetout->close();
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }
  
  file->close();
  pdatasetin->close();
  pdatasetout->close();
  return wrap(0);
}
***/

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

