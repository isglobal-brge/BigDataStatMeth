#include "include/svdBlockDecomposition_hdf5.h"


// Reads big matrix from hdf5 file in blocks and perform a svd descomposition from
// each block,results are saved in hdf5 datasets under temporal group to be processed
// if necessary
int First_level_SvdBlock_decomposition_hdf5(H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
                                            int irows, int icols, Rcpp::Nullable<int> threads = R_NilValue)
{
  
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  Eigen::MatrixXd nX;
  svdeig retsvd;
  int  M, p, n;
  int maxsizetoread;
  bool transp = false;
  std::string strGroupName  = "tmpgroup";
  Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2,icols);
  int cummoffset;
  

  if( exists_HDF5_element_ptr(file,strGroupName))
    remove_HDF5_element_ptr(file,strGroupName);

  int ret = create_HDF5_group_ptr(file, strGroupName);
  
  try{
    
    if(irows > icols) {
      // Work with transposed matrix
      n = icols;
      p = irows;
      //..// int i2 = 2;
      transp = true;

      // Get data to normalize matrix
      get_HDF5_mean_sd_by_column_ptr( file, dataset, datanormal);

    } else {
      n = irows;
      p = icols;
    }

    
    M = pow(k, q);
    if(M>p)
      throw std::runtime_error("k^q must not be greater than the number of columns of the matrix");
    
    double block_size = std::floor((double)p/(double)M); 
    int realsizeread;


    // Get data from M blocks in initial matrix
    for( int i = 0; i< M ; i++) 
    {
      //..// Rcpp::Rcout<<"\nObtaining k = "<< i+1 <<" of "<< M;
     
       
      // 1.- Get SVD from all blocks
      //    a) Get all blocks from initial matrix 
      
      IntegerVector offset = getInitialPosition( transp, (unsigned long long)(i*block_size) ); // Initial read position
      IntegerVector count; // Blocks size
      maxsizetoread = block_size;

      // Get max block size to read - for blocks smaller than default block size 
      if(transp == true)
      {
        if( ((i+1)*block_size) > irows)
          maxsizetoread = irows - (i*block_size);
        
        if( i+1 == M && irows - maxsizetoread!=0)
          realsizeread = irows - (i*block_size);
        else
          realsizeread = maxsizetoread;
          
      } else {
        
        if( ((i+1)*block_size) > icols)
          maxsizetoread = icols - (i*block_size);
        
        if( i+1 == M && icols - maxsizetoread!=0)
          realsizeread = icols - (i*block_size);
        else
          realsizeread = maxsizetoread;
      }

      count = getSizetoRead(transp, (unsigned long long)(realsizeread), icols, irows );

      Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, dataset, offset[0], offset[1], count[0], count[1]);
      

      if(transp==false)
        X.transposeInPlace();


      // Normalize data
      if (bcenter==true || bscale==true)
        X = RcppNormalize_Data_hdf5(X, bcenter, bscale, transp, datanormal);
  
      //    b) SVD for each block
      retsvd = RcppbdSVD_lapack(X, false, false);

      //    c)  U*d
      // Create diagonal matrix from svd decomposition d
      int isize = (retsvd.d).size();
      Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
      d.diagonal() = retsvd.d;
      

      std::string strDatasetName = strGroupName + "/A" + std::to_string(i/(M/k));
      
      Eigen::MatrixXd restmp = Bblock_matrix_mul_parallel(retsvd.u, d, 128, threads);
      

      //    d) Write results to hdf5 file
      offset[0] = 0; offset[1] = 0;
      if(transp == true )
      {
        count[0] = restmp.rows();
        count[1] = restmp.cols();
      }else {
        count[0] = restmp.cols();
        count[1] = restmp.rows();
      }

      
      if(i%(M/k) == 0) {
        // If dataset exists --> remove dataset
        if( exists_HDF5_element_ptr(file,strDatasetName))
          remove_HDF5_element_ptr(file,strDatasetName);
        
        
        // Create unlimited dataset in hdf5 file
        if( transp == true ){
          create_HDF5_unlimited_dataset_ptr(file, strDatasetName, count[0], count[1], "numeric");
        }
        else{
          create_HDF5_unlimited_dataset_ptr(file, strDatasetName, restmp.cols(), restmp.rows(), "numeric");
        }
        cummoffset = 0;
        
      } //else {
        
        
       
        
        // Get write position
          if(transp == true){
            offset[1] = cummoffset;
            cummoffset = cummoffset + restmp.cols();
          }
          else{
            offset[0] = cummoffset;
            cummoffset = cummoffset + restmp.rows();
          }

      //}
 
    
  /*** else {
        // Get write position
        if(maxsizetoread == block_size)
        {
          if(transp == true)
            offset[1] = (i%(M/k))*block_size;
          else
            offset[0] = (i%(M/k))*block_size;
        }
        else{
          if(transp==true)
            offset[1] = ( (i%(M/k))-1 )*block_size + maxsizetoread ;
          else
            offset[0] = ( (i%(M/k))-1 )*block_size + maxsizetoread ;
        }
      }***/
      
      DataSet* unlimDataset = new DataSet(file->openDataSet(strDatasetName));
      
      // Extend dataset before put data
      if((i%(M/k)) != 0)
      {
        if(transp == true)
          extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
        else
          extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[0]);
        
      }
      
      if(transp == true)
        write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp)  );  
      else
        write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp.transpose())  );  
      unlimDataset->close();
      
    }
    
    Rcpp::Rcout<<"\n";
    
    
  } catch(FileIException error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception (File IException)" );
    return -1;
  } catch(DataSetIException error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception (DataSet IException)" );
    return -1;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception (Group IException)" );
    return -1;
  } catch(DataSpaceIException error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception (DataSpace IException)" );
    return -1;
  } catch(DataTypeIException error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception (Data TypeIException)" );
    return -1;
  }catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    file->close();
    return -1;
  }
  
  return 0;
  
  
}



// Reads small datasets from hdf5 and perform a svd descomposition from each block,
// results are saved in hdf5 datasets under temporal group to be processed if necessary
int Next_level_SvdBlock_decomposition_hdf5(H5File* file, std::string strGroupName, int k, int q, 
                                           bool bcenter, bool bscale, Rcpp::Nullable<int> threads = R_NilValue)
{
  
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  IntegerVector count = IntegerVector::create(0, 0);
  IntegerVector offset = IntegerVector::create(0, 0);
  int cbefore = 0; //.. MODIFICAT 07/09/2020 ..//
  Eigen::MatrixXd nX;
  svdeig retsvd;
  int M;
  // int nconv, M, p, n;
  // int maxsizetoread;
  
  CharacterVector strvmatnames = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
  
  try{
    
    // Get dataset names
    StringVector joindata =  get_dataset_names_from_group(file, strGroupName, (std::string)strvmatnames[q-1]);
  
    M = joindata.size();
    
    // Get data from M blocks in initial matrix
    for( int i = 0; i< M ; i++) 
    {
      
      //    a) Get dataset
      DataSet currentdataset = file->openDataSet(strGroupName + "/" + joindata[i]);
      IntegerVector dims_out = get_HDF5_dataset_size(currentdataset);
      
      Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, &currentdataset, 0, 0, dims_out[0], dims_out[1]);
      
      //    b) Get dataset svd
      retsvd = RcppbdSVD_lapack(X, false, false);
      
      
      //    c) U*d
      int isize = (retsvd.d).size();
      Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
      d.diagonal() = retsvd.d;
      
      std::string strDatasetName = strGroupName + "/" + strvmatnames[q] + std::to_string(i/(M/k));
      
      Eigen::MatrixXd restmp = Bblock_matrix_mul_parallel(retsvd.u, d, 128, threads);
      
      //    d) Write results to dataset
      //.. MODIFICAT 07/09/2020 ..// int cbefore = restmp.cols();
      //.. MODIFICAT 07/09/2020 ..// offset[0] = 0; offset[1] = 0;
      count[0] = restmp.rows();
      count[1] = restmp.cols();
      
      if(i%(M/k) == 0) {
        // If dataset exists --> remove dataset
        if( exists_HDF5_element_ptr(file,strDatasetName))
          remove_HDF5_element_ptr(file,strDatasetName);
        // Create unlimited dataset in hdf5 file
        create_HDF5_unlimited_dataset_ptr(file, strDatasetName, restmp.rows(), restmp.cols(), "numeric");
        
      } else {
        // Get initial write position
        offset[1] = offset[1] + cbefore;
      }
      
      DataSet* unlimDataset = new DataSet(file->openDataSet(strDatasetName));
      // Extend dataset before put data
      if((i%(M/k)) != 0)
        extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
      
      cbefore = restmp.cols();
      
      write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp)  );  
      unlimDataset->close();
      
    }
    
    //..ONLY DEBUG !!!...//remove_HDF5_multiple_elements_ptr(file, strGroupName, joindata);
    
  }catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return -1;
  }
  
  return 0;
  
  
}

