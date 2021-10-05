#include "include/hdf5_splitMatrix.h"


// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;


// Internal call 
int RcppSplit_matrix_hdf5 ( H5File* file, DataSet* dataset, bool bycols, std::string stroutgroup, std::string stroutdataset, int blocksize, int irows, int icols )
{
  
  int blocks;
  int inrows = irows, 
      incols = icols,
      ii = 0,
      kk = 0;
  
  std::string newDatasetName = "";
  try {
    
    if( bycols == true ) {
      blocks = (icols + blocksize - 1) / blocksize;
      incols = blocksize;
    } else {
      blocks = (irows + blocksize - 1) / blocksize;
      inrows = blocksize;
    }
    
    for ( int i=0; i<blocks; i++)
    {
      newDatasetName = stroutgroup + "/" + stroutdataset + "." + std::to_string(i), incols;
      
      if( bycols == true) { 
        kk = i * blocksize;
        if( kk + blocksize > icols) { incols = icols - kk; }
      } else 
      {
        ii = i * blocksize;
        if( ii + blocksize > irows) { inrows = irows - ii; }
      }

      // Get block from complete matrix
      Eigen::MatrixXd Block = GetCurrentBlock_hdf5( file, dataset, kk, ii, incols, inrows);
      // Transform to rowmajor
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mapBlock(Block.transpose().data(), inrows, incols);
      write_HDF5_matrix_from_R_ptr(file, newDatasetName, Rcpp::wrap(mapBlock), false);
      
    }
    
  }catch( FileIException& error ) {
    file->close();
    dataset->close();
    ::Rf_error( "c++ exception (File IException )" );
    return -1;
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    file->close();
    dataset->close();
    ::Rf_error( "c++ exception (DataSet IException )" );
    return -1;
  } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
    file->close();
    dataset->close();
    ::Rf_error( "c++ exception (DataSpace IException )" );
    return -1;
  } 
  
  return(0);

  
}





//' Split hdf5 dataset
//'
//' Split hdf5 dataset by rows or columns and store splitted submatrices inside hdf5 file.
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param dataset, character array indicating the input dataset to be imputed
//' @param outgroup, optional character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, optional character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
//' @param nblocks, integer number of blocks in which we want to split the data
//' @param blocksize, integer, number of elements in each block
//' @param bycols, boolean by default = true, true indicates that the imputation will be done by columns, otherwise, the imputation will be done by rows
//' @param bforce, boolean if true, previous results in same location inside hdf5 will be overwritten.
//' @return Original hdf5 data file with imputed data
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSplit_matrix_hdf5( std::string filename, std::string group, std::string dataset, 
                                 Rcpp::Nullable<std::string> outgroup = R_NilValue, Rcpp::Nullable<std::string> outdataset = R_NilValue, 
                                 Rcpp::Nullable<int> nblocks = R_NilValue,  Rcpp::Nullable<int> blocksize = R_NilValue,
                                 Rcpp::Nullable<bool> bycols = true, Rcpp::Nullable<bool> force = false  )
{
  
  H5File* file;
  DataSet* pdataset = nullptr;
  
  try
  {
    std::string strdataset = group + "/" + dataset;
    std::string stroutgroup, stroutdataset, stroutdata;
    std::string strdatasetout;
    int iblocksize = 0, iwholesize = 0;
    bool bcols, bforce;
    
    
    if(bycols.isNull()) { bcols = true ;
    } else {   bcols = Rcpp::as<bool>(bycols);}
    
    if(force.isNull()) { bforce = false ;
    } else {   bforce = Rcpp::as<bool>(force);}
    
    if(outgroup.isNull()) {  stroutgroup = group ;
    } else {   stroutgroup = Rcpp::as<std::string>(outgroup);}
    
    if(outdataset.isNull()){  stroutdataset = dataset ;
    } else {   stroutdataset = Rcpp::as<std::string>(outdataset);}
    


    if (exist_FileGroupDataset (filename, group, dataset) != 0 ) {
      
      //..// strdataset = group +"/" + dataset;
      
      file = new H5File( filename, H5F_ACC_RDWR );
      pdataset = new DataSet(file->openDataSet(strdataset));
      
    }
    
    prepare_outGroup(file, stroutgroup, bforce);

    
    // Real data set dimension
    IntegerVector dims_out = get_HDF5_dataset_size(*pdataset);
    
    if(nblocks.isNull() && blocksize.isNull()){
      pdataset->close();
      file->close();
      Rcpp::Rcout<<"\n Block size or number of blocks needed to proceed with matrix split. Please, review parameters";
      return(wrap(-1));
      
    } else if (!nblocks.isNull() && !blocksize.isNull()) {
      pdataset->close();
      file->close();
      Rcpp::Rcout<<"\nBlock size and number of blocks are defined, please define only one option, split by number of blocks or by block size";
      return(wrap(-1));
      
    } else if(!nblocks.isNull()) {
      
      if ( Rcpp::as<int>(nblocks) == 1) {
        pdataset->close();
        file->close();
        Rcpp::Rcout<<"\nNumbers of blocks = 1, no data to split";
        return(wrap(-1));

      } else {
        
        double module;
        
        if(bcols == true) {
          iblocksize = dims_out[0] / Rcpp::as<int>(nblocks);
          module = dims_out[0] % iblocksize;
        } else {
          iblocksize = dims_out[1] / Rcpp::as<int>(nblocks);
          module = dims_out[1] % iblocksize;
        }
        if (module > 0) { iblocksize = iblocksize + 1; }
      }
      
    }else {
      iblocksize = Rcpp::as<int>(blocksize);
      
      if(bcols == true) {
        if( iblocksize == dims_out[0]) {  
          pdataset->close();
          file->close();
          throw std::range_error( "c++ exception in bdSplit_matrix_hdf5 ( No data to split)" ); 
        }
      } else {
        if( iblocksize == dims_out[1]) {  
          pdataset->close();
          file->close();
          throw std::range_error( "c++ exception in bdSplit_matrix_hdf5 ( No data to split)" ); 
        }
      }
    }
    
    
    RcppSplit_matrix_hdf5 ( file, pdataset, bcols, stroutgroup, stroutdataset, iblocksize, dims_out[1], dims_out[0] );

  }
  catch( FileIException& error ) { // catch failure caused by the H5File operations
    pdataset->close();
    file->close();
    ::Rf_error( "c++ exception (File IException)" );
    return(wrap(-1));
  }
  
  pdataset->close();
  file->close();
  
  return(wrap(0));
}



/***R

library(BigDataStatMeth)

setwd("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/UAB/DOCTORAT/BitDataStatMeth - BDSM/Analysis/BigDataStatMeth_Analysis/Cholesterol/test")

bdSplit_matrix_hdf5( "cars.hdf5", "data", "X", "dataoutCols", nblocks = 3, bycols = FALSE, force = TRUE)

bdSplit_matrix_hdf5( "cars.hdf5", "data", "X", "dataoutRows", nblocks = 3, bycols = TRUE, force = TRUE)



*/