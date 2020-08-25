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
    

    if( exists_HDF5_element_ptr(file, strSVDdataset ) )
      dataset = new DataSet(file->openDataSet(strSVDdataset));
    else
      throw ("Dataset does not exist !");
    
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
    error.printErrorStack();
    return -1;
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    error.printErrorStack();
    return -1;
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    error.printErrorStack();
    return -1;
  } 
  
  dataset->close();
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
//' @return original file with results in folder PCA/<datasetname>
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdPCA_hdf5(std::string filename, std::string group, std::string dataset)
{
  
  H5File* file;
  
  try
  {
    
    Eigen::MatrixXd X;
    Eigen::MatrixXd U, V, Lambda;
    Eigen::VectorXd lambda;
    Eigen::MatrixXd C, D, varcoord;
    
    pcaeig pcaRes;
    
    // Open an existing file and dataset.
    file = new H5File( filename, H5F_ACC_RDWR );

    std::string strSVDdataset = "SVD/" + dataset;

    // Look up for svd decomposition in hdf5 file
    if( !(exists_HDF5_element_ptr(file, strSVDdataset )) )
    {
      
      Rcpp::Rcout<<"\n Doncs ara ens toca calcular .... Fatal !!!!\n";
      svdeig retsvd = RcppbdSVD_hdf5( filename, group, dataset, 4, 1, 0, true, true );
    }
    
    get_HDF5_PCA_variance_ptr(file, dataset);
      
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
  } catch (...) {
    file->close();
    ::Rf_error("C++ exception (unknown reason)");
    return(wrap(-1));
  }
  
  file->close();
  return Rcpp::List::create(Rcpp::Named("filename") = filename);
  
}

