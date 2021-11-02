#include "include/hdf5_imputation.h"


// Get value for imputation
int get_value_to_impute_discrete(std::map<double, double> probMap)
{
  std::vector <double> probs;

  // Get values and counts for each map element
  for( auto it = probMap.begin(); it != probMap.end(); ++it )
    probs.push_back( it->second );

  // remove last element (corresponds to 3=<NA>)
  probs.erase(probs.end() - 1);

  // Get total count
  double totalSNPS = std::accumulate(probs.begin(), probs.end(), decltype(probs)::value_type(0));
  
  // Get probabilities without <NA>
  for (std::vector<double>::iterator it = probs.begin() ; it != probs.end(); ++it)
    *it = *it/totalSNPS;
  
  // Generate value with given probabilities
  std::random_device rd;
  std::mt19937 gen(rd());
  
  std::discrete_distribution<> d(probs.begin(), probs.end());
  
  return (d(gen));

}


// Convert NumericVector to map (key:vlues - value: frequency value in vector)
std::map<double, double> VectortoOrderedMap_SNP_counts( Eigen::VectorXd  vdata)
{
  std::map<double, double> mapv;
  
  try 
  {
    int position = 0;
    //.commented 20201120 - warning check().// int vlength = vdata.size();
    std::vector<double> v(vdata.data(), vdata.data()+vdata.size());

    std::sort(v.begin(), v.end() ); // Sort vector to optimize search and count
    
    for (size_t i = 0; i <=  *std::max_element(v.begin(), v.end()) ; ++i)  
    {
      double mycount = std::count(v.begin() + position, v.end(), i);
      mapv[i] = mycount;
      position = position + mycount;
    }


  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  return mapv;
}


// Pedestrian dataset imputation .... 
// TODO : 
//    - perform better imputation
void Impute_snp_HDF5(H5File* file, DataSet* dataset, bool bycols, std::string stroutdataset)
{
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  IntegerVector offset = IntegerVector::create(0, 0);
  IntegerVector count = IntegerVector::create(0, 0);
  DataSet* outdataset;
  int ilimit;
  int blocksize = 1000;
  
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
    
    
    if( stroutdataset.compare("")!=0)
    {
      hsize_t     dimsf[2];              // dataset dimensions
      dimsf[0] = dims_out[0];
      dimsf[1] = dims_out[1];
      
      DataSpace dataspace( RANK2, dimsf );
      outdataset = new DataSet(file->createDataSet(stroutdataset, PredType::NATIVE_DOUBLE, dataspace));
      
    } else  {
      outdataset = dataset;
    }
    
    
    
    for( int i=0; i<=(ilimit/blocksize); i++) 
    {
      int iread;
      
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
        for( int row = 0; row<data.rows(); row++)  // COMPLETE EXECUTION
        {
          std::map<double, double> myMap;
          myMap = VectortoOrderedMap_SNP_counts(data.row(row));
          
          /*** ORIGINAL FUNCIONA PERFECTAMENT PERÒ ASSIGNA UN MATEIX VALOR A TOS... !!!
           data.row(row) = (data.row(row).array() == 3).select( get_value_to_impute_discrete(myMap), data.row(row)); 
           ***/
          
          //..// data.row(row) = (data.row(row).array() == 0).select( -5, data.row(row)); 
          //..// data.row(row) = (data.row(row).array() == 1).select( 0, data.row(row)); 
          //..// data.row(row) = (data.row(row).array() == 2).select( 5, data.row(row)); 
          //..// data.row(row) = (data.row(row).array() == 3).select( 99, data.row(row)); 
          
          Eigen::VectorXd ev = data.row(row);
          std::vector<double> v(ev.data(), ev.data() + ev.size());
          
          auto it = std::find_if(std::begin(v), std::end(v), [](int i){return i == 3;});
          while (it != std::end(v)) {
            //..// results.emplace_back(std::distance(std::begin(v), it));
            if(*it==3) *it = get_value_to_impute_discrete(myMap);
            it = std::find_if(std::next(it), std::end(v), [](int i){return i == 3;});
          }
          
          Eigen::VectorXd X = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
          data.row(row) = X;
          
        }
        
      } else {
        for( int col = 0; col<data.cols(); col++) 
        {
          std::map<double, double> myMap;
          myMap = VectortoOrderedMap_SNP_counts(data.col(col));
          //..// data.col(col) = (data.col(col).array() == 3).select( get_value_to_impute_discrete(myMap), data.col(col));
          Eigen::VectorXd ev = data.col(col);
          std::vector<double> v(ev.data(), ev.data() + ev.size());
          
          auto it = std::find_if(std::begin(v), std::end(v), [](int i){return i == 3;});
          while (it != std::end(v)) {
            //..// results.emplace_back(std::distance(std::begin(v), it));
            if(*it==3) *it = get_value_to_impute_discrete(myMap);
            it = std::find_if(std::next(it), std::end(v), [](int i){return i == 3;});
          }
          
          Eigen::VectorXd X = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
          data.col(col) = X;
        }
      }
      //..// write_HDF5_matrix_subset_v2(file, dataset, offset, count, stride, block, wrap(data) );
      write_HDF5_matrix_subset_v2(file, outdataset, offset, count, stride, block, wrap(data) );
    }
    
    outdataset->close();
    
  } catch(FileIException& error) { // catch failure caused by the H5File operations
    outdataset->close();
    file->close();
    ::Rf_error( "c++ exception Impute_snp_HDF5 (File IException)" );
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
    outdataset->close();
    file->close();
    ::Rf_error( "c++ exception Impute_snp_HDF5 (DataSet IException)" );
  } catch(GroupIException& error) { // catch failure caused by the Group operations
    outdataset->close();
    file->close();
    ::Rf_error( "c++ exception Impute_snp_HDF5 (Group IException)" );
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    outdataset->close();
    file->close();
    ::Rf_error( "c++ exception Impute_snp_HDF5 (DataSpace IException)" );
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
    outdataset->close();
    file->close();
    ::Rf_error( "c++ exception Impute_snp_HDF5 (Data TypeIException)" );
  }
  
}




//' Impute SNPs in hdf5 omic dataset 
//'
//' Impute SNPs in hdf5 omic dataset 
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param dataset, character array indicating the input dataset to be imputed
//' @param bycols, boolean by default = true, true indicates that the imputation will be done by columns, otherwise, the imputation will be done by rows
//' @param outgroup, optional character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, optional character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
//' @return Original hdf5 data file with imputed data
//' @export
// [[Rcpp::export]]
void bdImpute_snps_hdf5(std::string filename, std::string group, std::string dataset, 
                              Rcpp::Nullable<std::string> outgroup = R_NilValue, Rcpp::Nullable<std::string> outdataset = R_NilValue, 
                              Rcpp::Nullable<bool> bycols = true )
{
  
  H5File* file;
  DataSet* pdataset = nullptr;
  
  try
  {
    std::string strdataset = group +"/" + dataset;
    std::string stroutgroup, stroutdataset, stroutdata;
    std::string strdatasetout;
    
    //.commented 20201120 - warning check().// int res;
    bool bcols;
    
    
    if(bycols.isNull())  bcols = true ;
    else    bcols = Rcpp::as<bool>(bycols);
    
    if(outgroup.isNull())  stroutgroup = group ;
    else    stroutgroup = Rcpp::as<std::string>(outgroup);
    
    if(outdataset.isNull())  stroutdataset = dataset ;
    else    stroutdataset = Rcpp::as<std::string>(outdataset);
    
    stroutdata = stroutgroup +"/" + stroutdataset;
    

    if(!ResFileExist_filestream(filename)){
      throw std::range_error("File not exits, create file before access to dataset");
    }
    
    file = new H5File( filename, H5F_ACC_RDWR );
    

    if(exists_HDF5_element_ptr(file, strdataset)) 
    {
      try
      {
        pdataset = new DataSet(file->openDataSet(strdataset));
        
        if( strdataset.compare(stroutdata)!= 0)
        {
          // If output is different from imput --> Remve possible existing dataset and create new
          if(exists_HDF5_element_ptr(file, stroutdata))
            remove_HDF5_element_ptr(file, stroutdata);
          
          // Create group if not exists
          if(!exists_HDF5_element_ptr(file, stroutgroup))
            file->createGroup(stroutgroup);
          
        } else {
          stroutdata = "";
        }
        
        Impute_snp_HDF5( file, pdataset, bcols, stroutdata);
        
      }catch(FileIException& error) {
        pdataset->close(); //.created 20201120 - warning check().//
        file->close();
      }
      
      
    } else{
      //.commented 20201120 - warning check().// pdataset->close();
      file->close();
      throw std::range_error("Dataset not exits");  
    }
    
    
  }
  catch( FileIException& error ) { // catch failure caused by the H5File operations
    //.commented 20201120 - warning check().// pdataset->close();
    file->close();
    ::Rf_error( "c++ exception (File IException)" );
    return void();
  }
  
  pdataset->close();
  file->close();
  Rcpp::Rcout<<"SNPs with missing values has been imputed\n";
  return void();
  
}




/***R

*/