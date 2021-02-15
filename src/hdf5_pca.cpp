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
int get_HDF5_PCA_variance_ptr(  H5File* file, std::string strdataset)
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
      ::Rf_error( "c++ exception (Dataset does not exist !)" );
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
    

  }catch( FileIException error ) {
    ::Rf_error( "c++ exception (File IException )" );
    return -1;
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception (DataSet IException )" );
    return -1;
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception (DataSpace IException )" );
    return -1;
  } 
  
  dataset->close();
  return(0);
}



// Get's variance and cumulative variance from svd decomposition
// var.contr, C, var.coord and var.cos^2 and write results to hdf5 file
int get_HDF5_PCA_variables_ptr(  H5File* file, std::string strdataset)
{
  
  DataSet* d;
  DataSet* v;
  
  try
  {
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1);
    IntegerVector offset = IntegerVector::create(0, 0);
    IntegerVector count = IntegerVector::create(0, 0);
    
    Exception::dontPrint();
    
    std::string strSVDdataset_d = "SVD/"+strdataset+"/d";
    std::string strSVDdataset_v = "SVD/"+strdataset+"/v";
    
    std::string strlocpcadataset = "PCA/" + strdataset;
    
    if( exists_HDF5_element_ptr(file, strSVDdataset_d ) ){
      d = new DataSet(file->openDataSet(strSVDdataset_d));
    }else{
      //.Remove note.// throw("Dataset does not exist !");
      // stop("Dataset does not exist !");
      ::Rf_error( "c++ exception (Dataset does not exist !)" );
    }
    
    // Real data set dimension
    IntegerVector dims_out = get_HDF5_dataset_size(*d);
    count[0] = dims_out[0];
    count[1] = dims_out[1];
    
    if(exists_HDF5_element_ptr(file, strlocpcadataset))
      remove_HDF5_element_ptr(file, strlocpcadataset);
    
    create_HDF5_groups_ptr(file, strlocpcadataset );
    
    Eigen::VectorXd data = GetCurrentBlock_hdf5(file, d, 0, 0, count[0], count[1]);
    Eigen::VectorXd vvar = data.array().pow(2);
    
    // Write lambda dataset
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/lambda", wrap(vvar));
    
    // Write variance dataset
    vvar = vvar/vvar.sum();
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/var", wrap(vvar));
    
    // Write cumulative variance dataset
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/cumvar", wrap(cumsum_hdf5(vvar)));
    
    if( exists_HDF5_element_ptr(file, strSVDdataset_v ) ){
      v = new DataSet(file->openDataSet(strSVDdataset_v));
    }else{
      //.Remove note.// throw("Dataset does not exist !");
      // stop("Dataset does not exist !");
      ::Rf_error( "c++ exception (Dataset does not exist !)" );
    }
    
    // Real data set dimension
    dims_out = get_HDF5_dataset_size(*v);
    count[0] = dims_out[0];
    count[1] = dims_out[1];
    
    Eigen::MatrixXd V = GetCurrentBlock_hdf5(file, v, 0, 0, count[0], count[1]);

    // Variable contribution var.contr
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/var.contr", wrap(V.pow(2)));
    
    // Correlations C
    Eigen::MatrixXd tmp =Bblock_matrix_mul_parallel(Eigen::MatrixXd::Identity(count[0], count[1]), V, 1024, R_NilValue);
    
    tmp =Bblock_matrix_mul_parallel( tmp, data.asDiagonal(), 1024, R_NilValue);
    
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/C", wrap(tmp));
    
    // Cosinus Variable var.cos2
    write_HDF5_matrix_ptr(file, strlocpcadataset+"/var.cos2", wrap(tmp.pow(2) ));
    
    // Quality Variable var.qual

    
  }catch( FileIException error ) {
    ::Rf_error( "c++ exception (File IException )" );
    return -1;
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception (DataSet IException )" );
    return -1;
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception (DataSpace IException )" );
    return -1;
  } 
  
  d->close();
  v->close();

  return(0);
}








// 
//' PCA Descomposition
//' 
//' Compute PCA
//' 
//' @param filename string, file name where dataset is stored 
//' @param group string group name  where dataset is stored in file
//' @param dataset string dataset name with data to perform PCA
//' @param bcenter logical value if true data is centered to zero
//' @param bscale logical value, if true data is scaled
//' @param threads integer number of threads used to run PCA
//' @return original file with results in folder PCA/<datasetname>
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdPCA_hdf5(std::string filename, std::string group, std::string dataset, 
                         Rcpp::Nullable<bool> bcenter = false, Rcpp::Nullable<bool> bscale = false, 
                         Rcpp::Nullable<int> threads)
{
  
  H5File* file;
  
  try
  {
    
    Eigen::MatrixXd X;
    Eigen::MatrixXd U, V, Lambda;
    Eigen::VectorXd lambda;
    Eigen::MatrixXd C, D, varcoord;
    // int ithreads=1;
    int ks = 4, qs = 1, nvs = 0;
    
    bool bcent= as<bool>(bcenter), 
         bscal = as<bool>(bscale);
    
    pcaeig pcaRes;
    
    // Open an existing file and dataset.
    file = new H5File( filename, H5F_ACC_RDWR );

    std::string strSVDdataset = "SVD/" + dataset;
    
    // Look up for svd decomposition in hdf5 file
    if( !(exists_HDF5_element_ptr(file, strSVDdataset )) )
    {
      
      /*
      if(threads.isNotNull()) 
      {
        if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency())
          ithreads = Rcpp::as<int> (threads);
        else 
          ithreads = std::thread::hardware_concurrency() - 1;
      }
      else    ithreads = std::thread::hardware_concurrency() - 1; //omp_get_max_threads(); 
      */
      //..// svdeig retsvd = RcppbdSVD_hdf5( filename, group, dataset, 4, 1, 0, true, true, ithreads );
      svdeig retsvd = RcppbdSVD_hdf5( filename, group, dataset, ks, qs, nvs, bcent, bscal, threads );
    }
    
    // Gets variance related variables
    //. Works ok but obsolete.//get_HDF5_PCA_variance_ptr(file, dataset);
    
    
    Rcpp::Rcout<<"\n Ara anem a la funció ?? !! \n";
    // Gets var.contr, C, var.coord and var.cos^2 
    get_HDF5_PCA_variables_ptr(file, dataset);
      
      //..// ans <- list(varcoord=var.coord, Y=Y, var = svdX$d^2, percvar=variance, components = svdX$v, method = method)
    
    /****

      {
        Eigen::VectorXd prov = lambda.array().sqrt();
        D = prov.asDiagonal();
      }
      
      // Get correlations
      {
        Rcpp::Rcout<<"\n Calculem produc per a C....: \n";
        Eigen::MatrixXd mprov = block_matrix_mul_parallel( Eigen::MatrixXd::Identity(V.rows(), V.cols()), V, 512 );
        C = block_matrix_mul_parallel( mprov, D, 512 );
      }

      // Get coordinates
      {
        svdeig singX = RcppbdSVD(X, k, int(), false);
        if(singX.bokuv == true && singX.bokd==true)
        {
          varcoord = block_matrix_mul_parallel(X.adjoint(), singX.u, 512 );  
        } else {
          throw(Rcpp::exception("Error with svd decomposition","pca.cpp",4));
        }
      }

      Lambda = lambda.asDiagonal(); 
    
    return Rcpp::List::create(Rcpp::Named("U") = U,
                              Rcpp::Named("V") = V,
                              Rcpp::Named("D") = D,
                              Rcpp::Named("Lambda") = Lambda,
                              Rcpp::Named("lambda") = lambda,
                              Rcpp::Named("pvac") = cumsum(lambda/lambda.array().sum()), // percentatge de vari`ancia acumulada pels components pvac,
                              Rcpp::Named("var.contr") = V.pow(2), // Variable contrib.
                              Rcpp::Named("C") = C,
                              Rcpp::Named("var.coord") = varcoord,
                              Rcpp::Named("var.cos2") = C.pow(2)
    );
    ***/
    

  } catch (std::exception &ex) {
    file->close();
    forward_exception_to_r(ex);
    return(wrap(-1));
  }catch( FileIException error ) {
    file->close();
    ::Rf_error( "c++ exception (File IException)" );
    return(wrap(-1));
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    file->close();
    ::Rf_error( "c++ exception (DataSet IException)" );
    return(wrap(-1));
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    file->close();
    ::Rf_error( "c++ exception (DataSpace IException)" );
    return(wrap(-1));
  } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
    file->close();
    ::Rf_error( "c++ exception (DataType IException)" );
    return(wrap(-1));
  } catch (...) {
    file->close();
    ::Rf_error("C++ exception (unknown reason)");
    return(wrap(-1));
  }
  
  file->close();
  return Rcpp::List::create(Rcpp::Named("filename") = filename);
  
}

