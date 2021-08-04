#include "include/hdf5_getAlleleFreq.h"

double calc_freq(NumericVector x)
{
  
  int len = x.size();
  
  std::vector<double> xc = Rcpp::as<std::vector<double> >(x);
  
  int n0 = std::count (xc.begin(), xc.end(), 0);
  int n1 = std::count (xc.begin(), xc.end(), 1);
  
  double maf = (double(n0)/len) + 0.5*(double(n1)/len);
  
  if( maf > 0.5 ) { 
    maf = 1 - maf;
  }
  
  return maf;

}






//' Get minor allele frequency
//' 
//' This function normalize data scaling, centering or scaling and centering in a dataset stored in hdf5 file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string or Delayed Array Matrix
//' @param dataset string or Delayed Array Matrix
//' @param byrows, boolean, default TRUE. If true, the frequency is calculated by rows, else, if byrows= FALSE, frequency is calculated by columns
//' @param bparallel, boolean, Perform calculous in parallel?, by default TRUE.
//' @param wsize integer (default = 1000), file block size to read to perform normalization
//' @return Numeric vector with allele frequencies
//' @examples
//' 
//' library(BigDataStatMeth)
//' 
//' maf_cols = resc <- bdget_allele_freq_hdf5("/Users/mailos/tmp/test/test.hdf5", 
//'                                           "test", "mat1", byrows = FALSE )
//' maf_rows = resc <- bdget_allele_freq_hdf5("/Users/mailos/tmp/test/test.hdf5", 
//'                                           "test", "mat1", byrows = TRUE )
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdget_maf_hdf5( std::string filename, const std::string group, std::string dataset,
                                  Rcpp::Nullable<bool> byrows = R_NilValue, Rcpp::Nullable<bool> bparallel  = R_NilValue,
                                  Rcpp::Nullable<int> wsize  = R_NilValue)
{
  
  H5File* file;
  
  DataSet* pdatasetin;

  Rcpp::NumericVector freqs;
  
  // int blocksize;
  bool readbyRows = true;
  
  try{
    
    
    if(!ResFileExist(filename)) {
      Rcpp::Rcout<<"\nFile not exits, create file before get allele frequencies\n";  
      return wrap(-1);
    }
    file = new H5File( filename, H5F_ACC_RDWR );
    
    
    if(exists_HDF5_element_ptr(file, group)==0) {
      Rcpp::Rcout<<"\nGroup not exits, create file and dataset before get allele frequencies\n";
      file->close();
      return Rcpp::wrap(-1);
    }  else{
      if(!exists_HDF5_element_ptr(file, group + "/" + dataset)) {
        Rcpp::Rcout<<"\n Dataset not exits, create file and dataset before get allele frequencies \n";
        file->close();
        return Rcpp::wrap(-1);
      }
    }
    
    
    pdatasetin = new DataSet(file->openDataSet(group + "/" + dataset));
    
    IntegerVector dims_out = get_HDF5_dataset_size(*pdatasetin);
    
    // // Define blocksize atending number of elements in rows and cols
    // if( wsize.isNull()) {
    //   int maxsize = std::max( dims_out[0], dims_out[1]);
    //   blocksize = std::ceil( maxElemBlock / maxsize);
    // } else {
    //   blocksize = Rcpp::as<int> (wsize);
    // }
    
    int readsize = 0;
    if( !byrows.isNull()) {
      readbyRows = Rcpp::as<bool> (byrows);
    }
    
    if(readbyRows == true) {
      readsize = dims_out[1];
    } else {
      readsize = dims_out[0];
    } 
      
    
    IntegerVector stride = IntegerVector::create(1, 1) ;
    IntegerVector block = IntegerVector::create(1, 1) ;
    IntegerVector offset = IntegerVector::create(0, 0) ;
    IntegerVector count = IntegerVector::create(0, 0) ;
    
    for( int i = 0; i< readsize; i++)
    {
      int length; 
      if(readbyRows == true) {
        offset[0] = 0; offset[1] =  i; 
        count[0] = dims_out[0]; count[1] =  1; 
        length = dims_out[0];
        
      } else {
        
        offset[0] = i; offset[1] =  0; 
        count[0] = 1; count[1] =  dims_out[1]; 
        length = dims_out[1];
      }
      
      NumericVector data( length);
      
      read_HDF5_matrix_subset(file, pdatasetin, offset, count, stride, block, REAL(data));
      freqs.push_back(calc_freq(data)); 
      
    }
    
    pdatasetin->close();
    file->close();
    return(freqs);
    
  } catch( FileIException error ) { // catch failure caused by the H5File operations
    pdatasetin->close();
    file->close();
    ::Rf_error( "c++ exception bdAllele_freq_hdf5 (File IException)" );
    return wrap(-1);
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    pdatasetin->close();
    file->close();
    ::Rf_error( "c++ exception bdAllele_freq_hdf5 (DataSet IException)" );
    return wrap(-1);
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    pdatasetin->close();
    file->close();
    ::Rf_error( "c++ exception bdAllele_freq_hdf5 (DataSpace IException)" );
    return wrap(-1);
  } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
    pdatasetin->close();
    file->close();
    ::Rf_error( "c++ exception bdAllele_freq_hdf5 (DataType IException)" );
    return wrap(-1);
  }catch(std::exception &ex) {
    pdatasetin->close();
    file->close();
    Rcpp::Rcout<< ex.what();
    return wrap(-1);
  }
  
  
}




/***R

library(BigDataStatMeth)

library(devtools)
reload(pkgload::inst("BigDataStatMeth"))

resc <- bdget_allele_freq_hdf5("test.hdf5", "test", "mat1", byrows = FALSE )

resr <- bdget_allele_freq_hdf5("test.hdf5", "test", "mat1", byrows = TRUE )



*/