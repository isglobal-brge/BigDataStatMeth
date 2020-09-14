#include "include/hdf5_QC.h"

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();
  
  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;
  
  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}



// Removes row or column with high missing data percentage
int Remove_snp_low_data_HDF5( H5File* file, DataSet* dataset, bool bycols, std::string stroutdata, double pcent)
{
  
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  IntegerVector offset = IntegerVector::create(0, 0);
  IntegerVector newoffset = IntegerVector::create(0, 0);
  IntegerVector count = IntegerVector::create(0, 0);
  DataSet* unlimDataset;
  int ilimit;
  int blocksize = 1000;
  int itotrem = 0;
  
  
  
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
      int actualrow = 0;
      int readedrows = data.rows();
      //..// for( int row = 0; row<readedrows; row++)  // COMPLETE EXECUTION
      for( int row = readedrows-1 ; row>=0; row--)  // COMPLETE EXECUTION
      {
          if((data.row(row).array() == 3).count()/(double)count[1]>= pcent )
          {
            removeRow(data, row);
            iblockrem = iblockrem + 1;
          } 
      }
    } else {
      
      int actualcol = 0;
      int readedcols = data.cols();
      
      //..//for( int col = 0; col<data.cols(); col++) 
      for( int col = readedcols-1 ; col>=0; col--)  // COMPLETE EXECUTION
      { 
        
        if((data.col(col).array() == 3).count()/(double)count[1]> pcent )
        {
          removeRow(data, col);
          iblockrem = iblockrem + 1;
        } 
        
        /***
        std::map<double, double> mapSNP;
        mapSNP = VectortoOrderedMap_SNP_counts(data. col(col));
        auto it = mapSNP.find(3);
        
        if(!(it == mapSNP.end()))
        {
          if( ( (double)it->second /(double)count[1] ) >pcent )
          {
            removeColumn(data, i-iblockrem);
            iblockrem = iblockrem + 1;
          }
        }
         ***/
      }
    }
    
    int extendcols = data.cols();
    int extendrows = data.rows();
    

    if(i==0) {
      create_HDF5_unlimited_dataset_ptr(file, stroutdata, extendrows, extendcols, "numeric");
      unlimDataset = new DataSet(file->openDataSet(stroutdata));
    }else {
      if(bycols == true)
        extend_HDF5_matrix_subset_ptr(file, unlimDataset, extendrows, 0);
      else
        extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, extendcols);
    }
    
    IntegerVector countblock = IntegerVector::create(extendrows, extendcols);
    
    write_HDF5_matrix_subset_v2(file, unlimDataset, newoffset, countblock, stride, block, wrap(data) );
    
    if(bycols == true)
      newoffset[0] =  newoffset[0] + extendrows;
    else
      newoffset[1] =  newoffset[2] + extendcols;
      
    
    itotrem = itotrem - iblockrem;
    
  }
  
  unlimDataset->close();
  
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
//' @param pcent, by default pcent = 0.5. Numeric indicating the percentage to be considered to remove SNPs, SNPS with percentage equal or higest will be removed from data
//' @param SNPincols, boolean by default = true, if true, indicates that SNPs are in cols, if SNPincols = false indicates that SNPs are in rows.
//' @return Original hdf5 data file with imputed data
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdRemovelowdata( std::string filename, std::string group, std::string dataset, std::string outgroup, std::string outdataset, 
                         Rcpp::Nullable<double> pcent, Rcpp::Nullable<bool> SNPincols )
{
  
  H5File* file;
  int iremoved = 0;
    
  try
  {
    bool bcols;
    double dpcent;
    DataSet* pdataset;
    
    std::string stroutdata = outgroup +"/" + outdataset;
    std::string strdataset = group +"/" + dataset;
    
    
    if(SNPincols.isNull())  bcols = true ;
    else    bcols = Rcpp::as<bool>(SNPincols);
    
    if(pcent.isNull())  dpcent = 0.5 ;
    else    dpcent = Rcpp::as<double>(pcent);
    
    if(!ResFileExist(filename))
      throw std::range_error("File not exits, create file before impute dataset");  
    
    file = new H5File( filename, H5F_ACC_RDWR );

    
    if(exists_HDF5_element_ptr(file, strdataset)) 
    {
      pdataset = new DataSet(file->openDataSet(strdataset));
      
      if( strdataset.compare(stroutdata)!= 0)
      {
        Rcpp::Rcout<<"\n Input and output datasets will be different \n";
        // If output is different from imput --> Remve possible existing dataset and create new
        if(exists_HDF5_element_ptr(file, stroutdata))
          remove_HDF5_element_ptr(file, stroutdata);
        
        // Create group if not exists
        if(!exists_HDF5_element_ptr(file, outgroup))
          file->createGroup(outgroup);

      } else {
        throw std::range_error("Input and output dataset must be different");  
      }
      
      iremoved = Remove_snp_low_data_HDF5( file, pdataset, bcols, stroutdata, dpcent);
      
      Function warning("warning");
      if (SNPincols )
        warning( std::to_string(iremoved) + " Columns have been removed");
      else
        warning( std::to_string(iremoved) + " Rows have been removed");
      
      pdataset->close();
      
    } else{
      pdataset->close();
      file->close();
      throw std::range_error("Dataset does not exits");  
    }
    
  
  }catch( FileIException error ){ // catch failure caused by the H5File operations
    file->close();
    error.printErrorStack();
    return(wrap(-1));
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    file->close();
    error.printErrorStack();
    return(wrap(-1));
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    file->close();
    error.printErrorStack();
    return(wrap(-1));
  } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
    file->close();
    error.printErrorStack();
    return(wrap(-1));
  }catch(std::exception &ex) {
    file->close();
    Rcpp::Rcout<< ex.what();
    return(wrap(-1));
  }
  
  file->close();
  return(wrap(iremoved));
  
}





/*** R

*/
