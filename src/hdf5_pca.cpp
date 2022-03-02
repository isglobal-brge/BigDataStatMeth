#include "include/hdf5_pca.h"


Eigen::VectorXd cumsum_hdf5(Eigen::VectorXd x){
  // initialize an accumulator variable
  double acc = 0;
  // initialize the result vector
  Eigen::VectorXd res = Eigen::VectorXd::Zero(x.size());
  for(int i = 0; i < x.size(); i++){
    acc += x[i];
    res[i] = acc;
  }
  
  
  return res;
}


// Get's variance and cumulative variance from svd decomposition
// and write results to hdf5 file
void get_HDF5_PCA_variance_ptr(  H5File* file, std::string strdataset)
{
  
  DataSet* dataset;
  
  try
  {
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1);
    IntegerVector offset = IntegerVector::create(0, 0);
    IntegerVector count = IntegerVector::create(0, 0);
    
    Exception::dontPrint();
    
    std::string strSVDdataset = "SVD/"+strdataset+"/d";
    std::string strlocpcadataset = "PCA/" + strdataset;
    

    if( exists_HDF5_element_ptr(file, strSVDdataset ) ){
      dataset = new DataSet(file->openDataSet(strSVDdataset));
    }else{
      //.Remove note.// throw("Dataset does not exist !");
      // stop("Dataset does not exist !");
      file->close();
      ::Rf_error( "c++ exception in get_HDF5_PCA_variance_ptr (Dataset does not exist !)" );
      return void();
    }

    // Real data set dimension
    IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
    count[0] = dims_out[0];
    count[1] = dims_out[1];
    
    if(exists_HDF5_element_ptr(file, strlocpcadataset))
      remove_HDF5_element_ptr(file, strlocpcadataset);
    
    create_HDF5_groups_ptr(file, strlocpcadataset );
    
    Eigen::VectorXd data = GetCurrentBlock_hdf5(file, dataset, 0, 0, count[0], count[1]);
    Eigen::VectorXd vvar = data.array().pow(2);
    
    // Write lambda dataset
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/lambda", wrap(vvar));

    // Write variance dataset
    vvar = vvar/vvar.sum();
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/var", wrap(vvar));
    
    // Write cumulative variance dataset
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/cumvar", wrap(cumsum_hdf5(vvar)));
    

  }catch( FileIException& error ) {
    ::Rf_error( "c++ exception (File IException )" );
    return void();
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception (DataSet IException )" );
    return void();
  } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception (DataSpace IException )" );
    return void();
  } 
  
  dataset->close();
  return void();
}



// Get's variance and cumulative variance from svd decomposition
// var.contr, C, var.coord and var.cos^2 and write results to hdf5 file
void get_HDF5_PCA_variables_ptr(  H5File* file, std::string strdataset)
{
  
  DataSet* d_ds;
  DataSet* v_ds;
  DataSet* u_ds;
  
  try
  {
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1);
    IntegerVector offset = IntegerVector::create(0, 0);
    IntegerVector count = IntegerVector::create(0, 0);
    IntegerVector offset_v = IntegerVector::create(0, 0);
    IntegerVector count_v = IntegerVector::create(0, 0);
    
    Exception::dontPrint();
    
    std::string strSVDdataset_d = "SVD/"+strdataset+"/d";
    std::string strSVDdataset_u = "SVD/"+strdataset+"/u";
    std::string strSVDdataset_v = "SVD/"+strdataset+"/v";
    
    std::string strlocpcadataset = "PCA/" + strdataset;
    
    
    if( exists_HDF5_element_ptr(file, strSVDdataset_d ) ){
      d_ds = new DataSet(file->openDataSet(strSVDdataset_d));
    }else{
      //.Remove note.// throw("Dataset does not exist !");
      // stop("Dataset does not exist !");
      file->close();
      ::Rf_error( "c++ exception in get_HDF5_PCA_variables_ptr(Dataset does not exist !)" );
    }
    
    
    if( exists_HDF5_element_ptr(file, strSVDdataset_v ) ){
      v_ds = new DataSet(file->openDataSet(strSVDdataset_v));
    }else{
      //.Remove note.// throw("Dataset does not exist !");
      // stop("Dataset does not exist !");
      file->close();
      ::Rf_error( "c++ exception in get_HDF5_PCA_variables_ptr(Dataset does not exist !)" );
    }

    // Real data set dimension
    IntegerVector dims_out = get_HDF5_dataset_size(*d_ds);
    count[0] = dims_out[0];
    count[1] = dims_out[1];
    
    
    // Real data set dimension
    IntegerVector dims_out_v = get_HDF5_dataset_size(*v_ds);
    count_v[0] = dims_out_v[0];
    count_v[1] = dims_out_v[1];
    
    // // Get only first n components (if defined)
    // if(ncomponents < count[0] && ncomponents != 0) {
    //     count[0] = ncomponents;
    // }
    // 
    // if(ncomponents < count_v[0] && ncomponents != 0) {
    //     count_v[0] = ncomponents;
    // }
    
    Eigen::VectorXd d = GetCurrentBlock_hdf5(file, d_ds, 0, 0, count[0], count[1]);
    Eigen::MatrixXd v = GetCurrentBlock_hdf5(file, v_ds, 0, 0, count_v[0], count_v[1]);
    
    
    if(exists_HDF5_element_ptr(file, strlocpcadataset))
      remove_HDF5_element_ptr(file, strlocpcadataset);
    
    create_HDF5_groups_ptr(file, strlocpcadataset );
    
    Rcpp::Rcout<<"\nGetting Lambda";
    Eigen::VectorXd vvar = d.array().pow(2);
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/lambda", wrap(vvar));
    
    Rcpp::Rcout<<"\nGetting Variance";
    // {
      // Eigen::VectorXd dComplete = GetCurrentBlock_hdf5(file, d_ds, 0, 0, dims_out[0], dims_out[1]);
    vvar = d.array().pow(2)/d.array().pow(2).sum()  ;
    // }
    
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/variance", wrap(vvar));
    
    Rcpp::Rcout<<"\nGetting Cumulative Variance";
    int ielements = 0;
    if(vvar.size()>1000){  ielements = 1000;    }
    else{ ielements = vvar.size();    }
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/cumvar", wrap(cumsum_hdf5(vvar.head(ielements))));
    
    
    Rcpp::Rcout<<"\nGetting Variables - Coord";
    Eigen::MatrixXd var_coord = v.array().colwise() * d.array();
    write_HDF5_matrix_transposed_ptr(file, strlocpcadataset+"/var.coord", wrap(var_coord.transpose()));
    
    Rcpp::Rcout<<"\nGetting Variables - Cos2";
    Eigen::MatrixXd var_cos2 = var_coord.unaryExpr([](double d) {return std::pow(d, 2);});;
    write_HDF5_matrix_transposed_ptr(file, strlocpcadataset+"/var.cos2", wrap(var_cos2.transpose()));
    
    Rcpp::Rcout<<"\nGetting Components";
    
    // if( ncomponents > 0 && ncomponents < dims_out[0] ) {
    // 
    //   if( exists_HDF5_element_ptr(file, strSVDdataset_u ) ){
    //     u_ds = new DataSet(file->openDataSet(strSVDdataset_u));
    //   }else {
    //     file->close();
    //     ::Rf_error( "c++ exception in get_HDF5_PCA_variables_ptr(Dataset u does not exist !)" );
    //   }
    //   
    //   Eigen::MatrixXd u = GetCurrentBlock_hdf5(file, u_ds, 0, 0, ncomponents, dims_out_v[0]);
    //   write_HDF5_matrix_transposed_ptr(file, strlocpcadataset + "/components", wrap(u));
    //   u_ds->close();
    //   
    // } else {
      create_hardLink(file, "/" + strSVDdataset_u, strlocpcadataset+"/components");
    // }


  }catch( FileIException& error ) {
    d_ds->close();
    v_ds->close();
    u_ds->close();
    file->close();
    ::Rf_error( "c++ exception (File IException )" );
    return void();
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    d_ds->close();
    v_ds->close();
    u_ds->close();
    file->close();
    ::Rf_error( "c++ exception (DataSet IException )" );
    return void();
  } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
    d_ds->close();
    v_ds->close();
    u_ds->close();
    file->close();
    ::Rf_error( "c++ exception (DataSpace IException )" );
    return void();
  } 
  
  d_ds->close();
  v_ds->close();
  return void();
}





// // Get's variance and cumulative variance from svd decomposition
// // var.contr, C, var.coord and var.cos^2 and write results to hdf5 file
// int get_HDF5_PCA_individuals_ptr(  H5File* file, std::string strgroup, std::string strdataset, bool bcenter, bool bscale)
// {
//   
//   
//   
//   try
//   {
//     
//     if(bcenter==true){
//       // Get mean and scale used to normalize after SVT
//       if(exists_HDF5_element_ptr(file, strgroup + "/" + strdataset)){
//         dataset = new DataSet(file->openDataSet(strgroup + "/" + strdataset));
//       } else {
//         file->close();
//         throw std::range_error("Dataset not exits"); 
//       }
//       
//       IntegerVector dims_out = get_HDF5_dataset_size(dataset);
//       Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2,dims_out[1]);
//       get_HDF5_mean_sd_by_column_ptr( file, dataset, datanormal);
//     }
//     
//     Eigen::MatrixXd Y = 
//     
//     
//   }catch( FileIException error ) {
//     ::Rf_error( "c++ exception (File IException )" );
//     return -1;
//   } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
//     ::Rf_error( "c++ exception (DataSet IException )" );
//     return -1;
//   } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
//     ::Rf_error( "c++ exception (DataSpace IException )" );
//     return -1;
//   } 
//   
// 
//   return(0);
// }
// 
// 





// 
//' PCA Descomposition
//' 
//' Compute PCA
//' 
//' @param filename string, file name where dataset is stored 
//' @param group string group name  where dataset is stored in file
//' @param dataset string dataset name with data to perform PCA
//' @param ncomponents integer, number of components to be computed, by default ncomponents = 0, all components are computed
//' @param bcenter logical value if true data is centered to zero
//' @param bscale logical value, if true data is scaled
//' @param k number of local SVDs to concatenate at each level, performance parameter 
//' @param q number of levels to compute SVD for PCA, performance parameter
//' @param rankthreshold double, threshold used to determine the range of the array. The matrix rank is equal to the number of
//'  singular values different from the threshold. By default, threshold = 0 is used to get the matrix rank , but it can be
//'  changed to an approximation of 0.
//' @param force logical value, if true, the SVD is forced to be computed although the SVD exists
//' @param threads integer number of threads used to run PCA
//' @return original file with results in folder PCA/<datasetname>
//' @export
// [[Rcpp::export]]
void bdPCA_hdf5(std::string filename, std::string group, std::string dataset,
                Rcpp::Nullable<int> ncomponents = 0,
                Rcpp::Nullable<bool> bcenter = false, Rcpp::Nullable<bool> bscale = false, 
                Rcpp::Nullable<int> k=2, Rcpp::Nullable<int> q=1,
                Rcpp::Nullable<double> rankthreshold = 0.0,
                Rcpp::Nullable<bool> force = false, Rcpp::Nullable<int> threads = R_NilValue)
{
  

  H5File* file;
  int incomponents = 0;
  double dthreshold;
  
  try
  {
    
    Eigen::MatrixXd X;
    Eigen::MatrixXd U, V, Lambda;
    Eigen::VectorXd lambda;
    Eigen::MatrixXd C, D, varcoord;
    int ks, qs = 1, nvs = 0;
    
    
    if(ncomponents.isNull())  incomponents = 0 ;
    else    incomponents = Rcpp::as<int>(ncomponents);
    
    if(q.isNull())  qs = 1 ;
    else    qs = Rcpp::as<int>(q);
    
    if(k.isNull())  ks = 2 ;
    else    ks = Rcpp::as<int>(k);
    
    if(rankthreshold.isNull()) {  
        dthreshold = 0 ;
    } else {
        if( Rcpp::as<double>(rankthreshold) > 0.1 ) {
            Rcpp::Rcout<< "Threshold to big, please set threshold with value lower than 0.1";
            return void();
        } else if( Rcpp::as<double>(rankthreshold) < 0 ) {
            Rcpp::Rcout<< "Threshold must be a positive value near zero";
            return void();
        } else {
            dthreshold = Rcpp::as<double>(rankthreshold);
        }
    }
    
    bool bcent= as<bool>(bcenter), 
         bscal = as<bool>(bscale),
         bforce = as<bool>(force);
    
    pcaeig pcaRes;
    
    // Open an existing file and dataset.
    file = new H5File( filename, H5F_ACC_RDWR );

    std::string strSVDdataset = "SVD/" + dataset;

    // Look up for svd decomposition in hdf5 file or if we have to recompute again the SVD
    // With parameter true we force to write data to disk 
    if( !(exists_HDF5_element_ptr(file, strSVDdataset )) || bforce == true) {
      svdeig retsvd = RcppbdSVD_hdf5_ptr( file, group, dataset, ks, qs, incomponents, bcent, bscal, true, dthreshold, threads ); 
    }

    // ------------ Variables ----------------

    get_HDF5_PCA_variables_ptr(file, dataset);
    
    // ------------ Individuals ----------------

    
    /*** TO DO : 
     * 
     *    get_HDF5_PCA_individuals_ptr(file, group, dataset, bcenter, bscale);
     *    
     ***/


  } catch (std::exception &ex) {
    file->close();
    forward_exception_to_r(ex);
    return void();
  }catch( FileIException& error ) {
    file->close();
    ::Rf_error( "c++ exception bdPCA_hdf5 (File IException)" );
    return void();
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    file->close();
    ::Rf_error( "c++ exception bdPCA_hdf5 (DataSet IException)" );
    return void();
  } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
    file->close();
    ::Rf_error( "c++ exception bdPCA_hdf5 (DataSpace IException)" );
    return void();
  } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
    file->close();
    ::Rf_error( "c++ exception bdPCA_hdf5 (DataType IException)" );
    return void();
  } catch(std::exception &ex) {
      file->close();
      Rcpp::Rcout<< "C++ exception bdPCA_hdf5 : "<< ex.what();
  } catch (...) {
    file->close();
    ::Rf_error("C++ exception bdPCA_hdf5 (unknown reason)");
    return void();
  }
  
  file->close();
  return void();
  
}

