#include "include/hdf5_removeMaf.h"

// Removes row or column with high missing data percentage
int Remove_MAF_HDF5( H5File* file, DataSet* dataset, bool bycols, std::string stroutdata, double pcent, int blocksize)
{
  
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  IntegerVector offset = IntegerVector::create(0, 0);
  IntegerVector newoffset = IntegerVector::create(0, 0);
  IntegerVector count = IntegerVector::create(0, 0);
  DataSet* unlimDataset = nullptr;
  int ilimit;
  // int blocksize = 100;
  int itotrem = 0;
  bool bcreated = false;
  
  
  try{
    
    // Real data set dimension
    IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
    
    // id bycols == true : read all rows by group of columns ; else : all columns by group of rows
    if (bycols == true) {
      ilimit = dims_out[0];
      count[1] = dims_out[1];
      offset[1] = 0;
    } else {
      ilimit = dims_out[1];
      count[0] = dims_out[0];
      offset[0] = 0;
    };
    
    for( int i=0; i<=(ilimit/blocksize); i++) 
    {
      int iread;
      int iblockrem = 0;
      
      if( (i+1)*blocksize < ilimit) iread = blocksize;
      else iread = ilimit - (i*blocksize);
      
      if(bycols == true) {
        count[0] = iread; 
        offset[0] = i*blocksize;
      } else {
        count[1] = iread; 
        offset[1] = i*blocksize;
      }
      
      // read block
      Eigen::MatrixXd data = GetCurrentBlock_hdf5(file, dataset, offset[0], offset[1], count[0], count[1]);
      
      if(bycols == true) // We have to do it by rows
      {
        int readedrows = data.rows();

        for( int row = readedrows-1 ; row>=0; row--)
        {
          if( calc_freq(wrap(data.row(row))) <= pcent ) {
            removeRow(data, row);
            iblockrem = iblockrem + 1;
          } 
        }
        
      } else {
        
        int readedcols = data.cols();
        
        for( int col = readedcols-1 ; col>=0; col--)
        { 
          if( calc_freq(wrap(data.col(col))) <= pcent ) {
            removeColumn(data, col);
            iblockrem = iblockrem + 1;
          } 
          
        }
      }
      
      
      int extendcols = data.cols();
      int extendrows = data.rows();
      
      if( extendrows>0 && extendcols>0)
      {
        
        if(bcreated == false) {
          create_HDF5_unlimited_matrix_dataset_ptr(file, stroutdata, extendrows, extendcols, "numeric");
          unlimDataset = new DataSet(file->openDataSet(stroutdata));
          bcreated = true;
        }else {
          if(bycols == true){
            extend_HDF5_matrix_subset_ptr(file, unlimDataset, extendrows, 0);
          }else{
            extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, extendcols);
          }
        }
        
        IntegerVector countblock = IntegerVector::create(extendrows, extendcols);
        write_HDF5_matrix_subset_v2(file, unlimDataset, newoffset, countblock, stride, block, wrap(data) );
        
        if(bycols == true)
          newoffset[0] =  newoffset[0] + extendrows;
        else
          newoffset[1] =  newoffset[1] + extendcols;
      }
      
      
      itotrem = itotrem - iblockrem;
      
    }
    
    if (bcreated == true) {
      unlimDataset->close();
    }
    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    unlimDataset->close();
    ::Rf_error( "c++ exception Remove_MAF_HDF5 (File IException)" );
    return -1;
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    unlimDataset->close();
    ::Rf_error( "c++ exception Remove_MAF_HDF5 (DataSet IException)" );
    return -1;
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    unlimDataset->close();
    ::Rf_error( "c++ exception Remove_MAF_HDF5 (Group IException)" );
    return -1;
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    unlimDataset->close();
    ::Rf_error( "c++ exception Remove_MAF_HDF5 (DataSpace IException)" );
    return -1;
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    unlimDataset->close();
    ::Rf_error( "c++ exception Remove_MAF_HDF5 (Data TypeIException)" );
    return -1;
  }
  
  
  return(itotrem);
}






//' Remove SNPs in hdf5 omic dataset with low data
//'
//' Remove SNPs in hdf5 omic dataset with low data
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param dataset, character array indicating the input dataset to be imputed
//' @param outgroup, character array indicating group where the data set will be saved after remove data with if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
//' @param maf, by default maf = 0.05. Numeric indicating the percentage to be considered to remove SNPs, SNPS with higest MAF will be removed from data
//' @param bycols, boolean by default = true, if true, indicates that SNPs are in cols, if SNPincols = false indicates that SNPs are in rows.
//' @param blocksize, integer, block size dataset to read/write and calculate MAF, by default this operations is made in with 100 rows if byrows = true or 100 cols if byrows = false.
//' @return Original hdf5 data file with imputed data
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdremove_maf_hdf5( std::string filename, std::string group, std::string dataset, std::string outgroup, std::string outdataset, 
                               Rcpp::Nullable<double> maf, Rcpp::Nullable<bool> bycols, Rcpp::Nullable<int> blocksize )
{
  
  H5File* file = nullptr;
  int iremoved = 0;
  int iblocksize = 100;
  
  try
  {
    bool bcols;
    double dpcent;
    
    std::string stroutdata = outgroup +"/" + outdataset;
    std::string strdataset = group +"/" + dataset;
    
    if(bycols.isNull()){  
      bcols = false ;
    }else{    
      bcols = Rcpp::as<bool>(bycols);
    }
    
    if(maf.isNull()){  
      dpcent = 0.05 ;
    }else{    
      dpcent = Rcpp::as<double>(maf);
    }
    
    if(!blocksize.isNull()){  
      iblocksize = Rcpp::as<int>(blocksize);
    }
    
    

    if(!ResFileExist_filestream(filename)){
      throw std::range_error("File not exits, create file before access to dataset");
    }
    
    
    file = new H5File( filename, H5F_ACC_RDWR );
    
    if(exists_HDF5_element_ptr(file, strdataset)) 
    {
      
      DataSet* pdataset = nullptr;
      
      pdataset = new DataSet(file->openDataSet(strdataset));
      
      if( strdataset.compare(stroutdata)!= 0)
      {
        
        // If output is different from imput --> Remve possible existing dataset and create new
        if(exists_HDF5_element_ptr(file, stroutdata))
          remove_HDF5_element_ptr(file, stroutdata);
        
        // Create group if not exists
        if(!exists_HDF5_element_ptr(file, outgroup))
          file->createGroup(outgroup);
        
      } else {
        throw std::range_error("Input and output dataset must be different");  
      }
      
      iremoved = Remove_MAF_HDF5( file, pdataset, bcols, stroutdata, dpcent, iblocksize);
      
      Function warning("warning");
      if (!bycols )
        warning( std::to_string(iremoved) + " Rows have been removed");
      else
        warning( std::to_string(iremoved) + " Columns have been removed");
      
      pdataset->close();
      
    } else{
      //.commented 20201120 - warning check().// pdataset->close();
      file->close();
      throw std::range_error("Dataset does not exits");  
    }
    
    
  }catch( FileIException& error ){ // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception bdremove_maf_hdf5 (File IException)" );
    return(wrap(-1));
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    file->close();
    ::Rf_error( "c++ exception bdremove_maf_hdf5 (DataSet IException)" );
    return(wrap(-1));
  } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
    file->close();
    ::Rf_error( "c++ exception bdremove_maf_hdf5 (DataSpace IException)" );
    return(wrap(-1));
  } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
    file->close();
    ::Rf_error( "c++ exception bdremove_maf_hdf5 (DataType IException)" );
    return(wrap(-1));
  }catch(std::exception &ex) {
    file->close();
    Rcpp::Rcout<< ex.what();
    return(wrap(-1));
  }
  
  file->close();
  return(wrap(iremoved));
  
}







/***R

*/
