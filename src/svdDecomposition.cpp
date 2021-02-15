#include "include/svdDecomposition.h"
#include "include/ReadDelayedData.h"

// SVD decomposition 
svdeig RcppbdSVD( Eigen::MatrixXd& X, int k, int ncv, bool bcenter, bool bscale )
{
  
  svdeig retsvd;
  Eigen::MatrixXd nX;
  int nconv;
  
  if( k==0 )    k = (std::min(X.rows(), X.cols()))-1;
  else if (k > (std::min(X.rows(), X.cols()))-1 ) k = (std::min(X.rows(), X.cols()))-1;
  
  if(ncv == 0)  ncv = k + 1 ;
  if(ncv<k) ncv = k + 1;
  
  {
    Eigen::MatrixXd Xtcp;
    //..//if(normalize ==true )  {
    if(bcenter ==true || bscale == true)  {
      // Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(RcppNormalize_Data(X))));
      nX = RcppNormalize_Data(X, bcenter, bscale);
      Xtcp =  bdtcrossproduct(nX);
    }else {
      //Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(X)));
      Xtcp =  bdtcrossproduct(X);
      
    }
    
    Spectra::DenseSymMatProd<double> op(Xtcp);
    Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, ncv);
    
    // Initialize and compute
    eigs.init();
    nconv = eigs.compute();
    
    if(eigs.info() == Spectra::SUCCESSFUL)
    {
      retsvd.d = eigs.eigenvalues().cwiseSqrt();
      retsvd.u = eigs.eigenvectors();
      retsvd.bokuv = true;
    } else {
      retsvd.bokuv = false;
    }
    
  }
  
  if(retsvd.bokuv == true)
  {
    Eigen::MatrixXd Xcp;
    if(bcenter ==true || bscale==true )  {
      // Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(RcppNormalize_Data(X))));  
      Xcp =  bdcrossproduct(nX);  
    }else {
      // Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(X)));
      Xcp =  bdcrossproduct(X);
    }  

    Spectra::DenseSymMatProd<double> opv(Xcp);
    Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigsv(&opv, k, ncv);
    
    // Initialize and compute
    eigsv.init();
    nconv = eigsv.compute();
    
    // Retrieve results
    if(eigsv.info() == Spectra::SUCCESSFUL)
    {
      retsvd.v = eigsv.eigenvectors();
    } else {
      retsvd.bokd = false;
    }
    
  }
  
  return retsvd;
}




// Lapack SVD decomposition
svdeig RcppbdSVD_lapack( Eigen::MatrixXd& X, bool bcenter, bool bscale )
{
  
  svdeig retsvd;
  
  char Schar='S';
  char Achar='A';
  int info = 0;
  
  
   if(bcenter ==true || bscale == true)
     X = RcppNormalize_Data(X, bcenter, bscale);
  
  int m = X.rows();
  int n = X.cols();
  int lda = std::max(1,m);
  int ldu = std::max(1,m);
  int ldvt = std::min(m, n);
  int k = std::min(m,n);
  int lwork;
  
  /*if(n>5*m)
    lwork = std::max(1,5*std::min(m,n));
  else
    lwork = std::max( 1, 4*std::min(m,n)* std::min(m,n) + 7*std::min(m, n) );
  */
  lwork = std::max( 5*std::min(m,n)+ std::max(m,n), 9*std::min(m, n) ); //.. ORIGINAL ..//
  //..// lwork = std::max( 3*std::min(m,n)+ std::max(m,n), 5*std::min(m, n) );
  

  Eigen::VectorXd s = Eigen::VectorXd::Zero(k);
  Eigen::VectorXd work = Eigen::VectorXd::Zero(lwork);
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(ldu,k);
  Eigen::MatrixXd vt = Eigen::MatrixXd::Zero(ldvt,n);

  dgesvd_( &Schar, &Schar, &m, &n, X.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork, &info);

  //..// Rcpp::Rcout<<"\nDescomposició - d : \n"<<s;
  //..// Rcpp::Rcout<<"\nDescomposició - u : \n"<<u;
  
  retsvd.d = s;
  retsvd.u = u;
  retsvd.v = vt.transpose();
  
  return retsvd;
}





// ##' @param k number of local SVDs to concatenate at each level 
// ##' @param q number of levels
svdeig RcppbdSVD_hdf5_Block( H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
                             int irows, int icols, Rcpp::Nullable<int> threads )
{
  
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  svdeig retsvd;
  Eigen::MatrixXd nX;
  bool transp = false;
  std::string strGroupName  = "tmpgroup";
  std::string strPrefix;
  

  CharacterVector strvmatnames = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
  strPrefix = strvmatnames[q-1];
  
  try{
    
    if(irows >= icols)
      transp = true;
    
    First_level_SvdBlock_decomposition_hdf5( file, dataset, k, q, nev, bcenter, bscale, irows, icols, threads);
    

    for(int j = 1; j < q; j++) // For each decomposition level : 
    {
      Next_level_SvdBlock_decomposition_hdf5(file, strGroupName, k, j, bcenter, bscale, threads);
    }
    

    // Get dataset names
    StringVector joindata =  get_dataset_names_from_group(file, strGroupName, strPrefix);

    // 1.- Join matrix and remove parts from file
    std::string strnewdataset = std::string((joindata[0])).substr(0,1);
    join_datasets(file, strGroupName, joindata, strnewdataset);
    //..ONLY DEBUG !!!...//remove_HDF5_multiple_elements_ptr(file, strGroupName, joindata);

    // 2.- Get SVD from Blocks full mattrix

    DataSet datasetlast = file->openDataSet(strGroupName + "/" + strnewdataset);
    IntegerVector dims_out = get_HDF5_dataset_size(datasetlast);

    
    Eigen::MatrixXd matlast;
    if(transp==true){
      matlast = GetCurrentBlock_hdf5(file, &datasetlast, 0, 0, dims_out[0],dims_out[1]);}
    else
      matlast = GetCurrentBlock_hdf5_Original(file, &datasetlast, 0, 0, dims_out[0],dims_out[1]);

    retsvd = RcppbdSVD_lapack(matlast, false, false);
    
    // Write results to hdf5 file : in folder "SVD" and dataset "SVD".<name input dataset>
    // Create structure and write d 
    StringVector name = get_dataset_names_from_dataset_ptr(dataset);
    create_HDF5_groups_ptr(file,"SVD/"+ name[0]);
    write_HDF5_matrix_ptr(file, "SVD/"+ name[0]+"/d", wrap(retsvd.d));
    
    // 3.- crossprod initial matrix and svdA$u
    IntegerVector dims_out_first = get_HDF5_dataset_size(*dataset);
    // Get initial matrix
    Eigen::MatrixXd A = GetCurrentBlock_hdf5(file, dataset, 0, 0,dims_out_first[0], dims_out_first[1] );
    
    //..// Rcpp::Rcout<< "\n Dimensions A : \n\t"<<A.cols()<<" x "<<A.rows()<<"\n Dimensions u : "<<(retsvd.u).rows()<<" x "<< (retsvd.v).cols()<<"\n";
    
    Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2,A.cols());

    Eigen::MatrixXd v;    
    if(transp==true){
      // Get data to normalize matrix
      get_HDF5_mean_sd_by_column_ptr( file, dataset, datanormal);
      
      if(bcenter == true || bscale==true)    
        A = RcppNormalize_Data_hdf5(A, bcenter, bscale, transp, datanormal );
      
      v = Bblock_matrix_mul_parallel(A.transpose(),retsvd.u,128, threads); // ABANS no PARALLEL
    } else{
      if(bcenter == true || bscale==true)    
        A = RcppNormalize_Data_hdf5(A.transpose(), bcenter, bscale, transp, datanormal );
      
      v = Bblock_matrix_mul_parallel(A.transpose(),retsvd.u,128, threads); // ABANS no PARALLEL
    }
    // 4.- resuls / svdA$d
    v = v.array().rowwise()/(retsvd.d).transpose().array();    
    
 
    // 5.- Get u and v : 
    //        --> If trans : u=v i v=u  
    
    if (transp == true)  {
      retsvd.v = retsvd.u;
      retsvd.u = v;
    } else {
      retsvd.v = v;
    }


    write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/u", wrap(retsvd.u));
    write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/v", wrap(retsvd.v));
  
  
  } catch(FileIException error) { // catch failure caused by the H5File operations
    ::Rf_error( "c++ exception (File IException)" );
    return retsvd;
  } catch(DataSetIException error) { // catch failure caused by the DataSet operations
    ::Rf_error( "c++ exception (DataSet IException)" );
    return retsvd;
  } catch(GroupIException error) { // catch failure caused by the Group operations
    ::Rf_error( "c++ exception (Group IException)" );
    return retsvd;
  } catch(DataSpaceIException error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception (DataSpace IException)" );
    return retsvd;
  } catch(DataTypeIException error) { // catch failure caused by the DataSpace operations
    ::Rf_error( "c++ exception (Data TypeIException)" );
    return retsvd;
  }catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
  }
  
  return retsvd;
}



// SVD decomposition with hdf5 file
//    input data : hdf5 file (object from crossproduct matrix) datagroup = 'strsubgroupIN'
//    output data : hdf5 file svd data in datagroup svd 
//                        svd/d 
//                        svd/u 
//                        svd/v 
//                        
//  https://github.com/isglobal-brge/svdParallel/blob/8b072f79c4b7c44a3f1ca5bb5cba4d0fceb93d5b/R/generalBlockSVD.R
//  @param k number of local SVDs to concatenate at each level 
//  @param q number of levels
//  
svdeig RcppbdSVD_hdf5( std::string filename, std::string strsubgroup, std::string strdataset,  
                       int k, int q, int nev, bool bcenter, bool bscale, Rcpp::Nullable<int> ithreads = R_NilValue )
{
  
  svdeig retsvd;
  Eigen::MatrixXd X;
  // int nconv;

  // Open an existing file and dataset.
  H5File file(filename, H5F_ACC_RDWR);
  DataSet dataset;
  
  if(exists_HDF5_element_ptr(&file, strsubgroup + "/" + strdataset)){
    dataset = file.openDataSet(strsubgroup + "/" + strdataset);
  }
  else {
    file.close();
    throw std::range_error("Dataset not exits"); 
  }

  // Get dataset dims
  IntegerVector dims_out = get_HDF5_dataset_size(dataset);

  hsize_t offset[2] = {0,0};
  //..// hsize_t count[2] = {as<hsize_t>(dims_out[0]), as<hsize_t>(dims_out[1])};
  hsize_t count[2] = { (unsigned long long)dims_out[0], (unsigned long long)dims_out[1]};

  // In memory computation for small matrices (rows or columns<5000)
  // Block decomposition for big mattrix
  if( std::max(dims_out[0], dims_out[1]) < MAXSVDBLOCK )
  {
    
    X = GetCurrentBlock_hdf5( &file, &dataset, offset[0], offset[1], count[0], count[1]);
    X.transposeInPlace();

    retsvd = RcppbdSVD(X, k, nev, bcenter, bscale);

  }
  else{

    // data stored transposed in hdf5
    int xdim = (unsigned long long)dims_out[1];
    int ydim = (unsigned long long)dims_out[0];
    
    // Remove previous results
    if(exists_HDF5_element_ptr(&file,"SVD/"+strdataset))
    {
      //..// Rcpp::Rcout<<"\n Old DATASET have been removed \n";
      remove_HDF5_element_ptr(&file,"SVD/"+strdataset);
    }
      
    retsvd = RcppbdSVD_hdf5_Block( &file, &dataset, k, q, nev, bcenter, bscale, xdim, ydim, wrap(ithreads));
  }
  

  return retsvd;
   
}





svdeig RcppCholDec(const Eigen::MatrixXd& X)
{
  Eigen::MatrixXd mX = X;
  svdeig decomp;
  
  Eigen::LDLT<Eigen::MatrixXd> cholSolv = mX.ldlt();

  // if symetric + positive definite -> info = Success
  // else no Cholesky Decomposition --> svd decomposition with Spectra
  if(cholSolv.info()==Eigen::Success)
  {
    size_t n = cholSolv.cols();
    decomp.d = cholSolv.vectorD();
    Eigen::MatrixXd preinv = Eigen::MatrixXd::Identity(n, n);
    decomp.v = cholSolv.solve(preinv);
  } else {
    Rcpp::Rcout<<"No symetric positive matrix, Cholesky decomposition not viable.";
    decomp = RcppbdSVD(mX, int(), int(), false);
  }
  return(decomp);
}




//' Inverse Cholesky of Delayed Array
//' 
//' This function get the inverse of a numerical or Delayed Array matrix. If x is hermitian and positive-definite matrix then 
//' performs get the inverse using Cholesky decomposition
//' 
//' 
//' @param x numerical or Delayed Array matrix. If x is Hermitian and positive-definite performs
//' @return inverse matrix of d 
//' @examples
//' 
//' library(DelayedArray)
//' 
//' A <- matrix(c(3,4,3,4,8,6,3,6,9), byrow = TRUE, ncol = 3)
//' bdInvCholesky(A)
//' 
//' # with Delayed Array
//' DA <- DelayedArray(A)
//' bdInvCholesky(DA)
//' 
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd bdInvCholesky (const Rcpp::RObject & x )
{
  
  svdeig result;
  Eigen::MatrixXd X;
  
  if ( x.isS4() == true)    
  {
    X = read_DelayedArray(x);
  } else {
    try{  
      X = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(x);
    }
    catch(std::exception &ex) { }
  }
  
  result = RcppCholDec(X);
  
  return (result.v);
  
}



//' k first SVD components for DelayedArray 
//' 
//' This function gets k first components from svd decomposition of numerical or Delayed Array 
//' 
//' @param x numerical or Delayed Array matrix
//' @param k number of eigen values , this should satisfy k = min(n, m) - 1
//' @param nev (optional, default nev = n-1) Number of eigenvalues requested. This should satisfy 1≤ nev ≤ n, where n is the size of matrix. 
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
//' @return u eigenvectors of AA^t, mxn and column orthogonal matrix
//' @return v eigenvectors of A^tA, nxn orthogonal matrix
//' @return d singular values, nxn diagonal matrix (non-negative real values)
//' @examples
//' 
//' library(DelayedArray)
//' 
//' n <- 500
//' A <- matrix(rnorm(n*n), nrow=n, ncol=n)
//' AD <- DelayedArray(A)
//' 
//' # svd without normalization
//' decsvd <- bdSVD( A, bscale = FALSE, bcenter = FALSE ) # No matrix normalization
//' decsvd$d
//' decsvd$u
//' 
//' # svd with normalization
//' decvsd <- bdSVD( A, bscale = TRUE, bcenter = TRUE) # Matrix normalization
//' 
//' decsvd$d
//' decsvd$u
//' 
//' # svd with scaled matrix (sd)
//' decvsd <- bdSVD( A, bscale = TRUE, bcenter = FALSE) # Scaled matrix
//' 
//' decsvd$d
//' decsvd$u
//' # svd with centered matrix (sd)
//' decvsd <- bdSVD( A, bscale = FALSE, bcenter = TRUE) # Centered matrix
//' decsvd$d
//' decsvd$u
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSVD (const Rcpp::RObject & x, Rcpp::Nullable<int> k=0, Rcpp::Nullable<int> nev=0,
                     Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true)
{
  
  auto dmtype = beachmat::find_sexp_type(x);
  int ks, nvs;
  bool bcent, bscal;
  
  if(k.isNull())  ks = 0 ;
  else    ks = Rcpp::as<int>(k);
  
  if(nev.isNull())  nvs = 0 ;
  else    nvs = Rcpp::as<int>(nev);
  
  
  if(bcenter.isNull())  bcent = true ;
  else    bcent = Rcpp::as<bool>(bcenter);
  
  if(bscale.isNull())  bscal = true ;
  else    bscal = Rcpp::as<bool>(bscale);
  
  // size_t ncols = 0, nrows=0;
  Eigen::MatrixXd X;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( x.isS4() == true){
      X = read_DelayedArray(x);
    }else {
      try{
        X = Rcpp::as<Eigen::MatrixXd >(x);
      }catch(std::exception &ex) {
        X = Rcpp::as<Eigen::VectorXd >(x);
      }
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  svdeig retsvd;
  retsvd = RcppbdSVD(X,ks,nvs, bcent, bscal);
  
  ret["u"] = retsvd.u;
  ret["v"] = retsvd.v;
  ret["d"] = retsvd.d;
  
  return Rcpp::wrap(ret);
  
}




//' Block SVD decomposition using an incremental algorithm.
//'
//' Singular values and left singular vectors of a real nxp matrix 
//' @title Block SVD decomposition using an incremental algorithm.
//' @param x a real nxp matrix in hdf5 file
//' @param group group in hdf5 data file where dataset is located
//' @param dataset matrix dataset with data to perform SVD
//' @param k number of local SVDs to concatenate at each level 
//' @param q number of levels
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
//' @param threads (optional) only used in some operations inside function. If threads is null then threads =  maximum number of threads available - 1.
//' @return a list of three components with the singular values and left and right singular vectors of the matrix
//' @return A List with : 
//' \itemize{
//'   \item{"u"}{ eigenvectors of AA^t, mxn and column orthogonal matrix }
//'   \item{"v"}{ eigenvectors of A^tA, nxn orthogonal matrix }
//'   \item{"v"}{ singular values, nxn diagonal matrix (non-negative real values) }
//' }
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSVD_hdf5 (const Rcpp::RObject & x, Rcpp::Nullable<CharacterVector> group = R_NilValue, 
                          Rcpp::Nullable<CharacterVector> dataset = R_NilValue,
                          Rcpp::Nullable<int> k=2, Rcpp::Nullable<int> q=1,
                          Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true,
                          Rcpp::Nullable<int> threads = R_NilValue)
{
  
  std::string filename;
  svdeig retsvd;
  
  try {
    
    int ks, qs, nvs = 0;
    bool bcent, bscal;
    CharacterVector strgroup, strdataset;
    // int ithreads=1;
    
    if(k.isNull())  ks = 2 ;
    else    ks = Rcpp::as<int>(k);
    
    if(q.isNull())  qs = 1 ;
    else    qs = Rcpp::as<int>(q);
    
    if(bcenter.isNull())  bcent = true ;
    else    bcent = Rcpp::as<bool>(bcenter);
    
    if(bscale.isNull())  bscal = true ;
    else    bscal = Rcpp::as<bool>(bscale);
    
    if(group.isNull())  strgroup = "" ;
    else    strgroup = Rcpp::as<std::string>(group);
    
    if(dataset.isNull())  strdataset = "";
    else    strdataset = Rcpp::as<std::string>(dataset);
   
   /* 
    if(threads.isNull())  ithreads = std::thread::hardware_concurrency() - 1;
    else    ithreads = Rcpp::as<int>(threads);*/
    
    if(is<CharacterVector>(x))
      filename = as<std::string>(x);
    else
      throw std::invalid_argument("File name must be character string");
      
      // throw std::runtime_error("unacceptable matrix type");
      
      // Rcpp::Rcout<<"Abans de cridar el procés del svd... k val : "<<ks<<"\n";
      
    retsvd = RcppbdSVD_hdf5( filename, as<std::string>(strgroup), as<std::string>(strdataset), ks, qs, nvs, bcent, bscal, threads );
    
    
  }catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return List::create(Named("d") = R_NilValue,
                        Named("u") = R_NilValue,
                        Named("v") = R_NilValue,
                        Named("file") = R_NilValue);
  }

  return List::create(Named("d") = retsvd.d,
                      Named("u") = retsvd.u,
                      Named("v") = retsvd.v,
                      Named("file") = filename);
  
}


//' Complete SVD with Lapack Functions for DelayedArray and RObjects
//' 
//' This function performs a complete svd decomposition of numerical matrix or Delayed Array with 
//' 
//' @param x numerical or Delayed Array matrix
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
//' @return u eigenvectors of AA^t, mxn and column orthogonal matrix
//' @return v eigenvectors of A^tA, nxn orthogonal matrix
//' @return d singular values, nxn diagonal matrix (non-negative real values)
//' @examples
//' 
//' library(DelayedArray)
//' 
//' n <- 500
//' A <- matrix(rnorm(n*n), nrow=n, ncol=n)
//' AD <- DelayedArray(A)
//' 
//' # svd without normalization
//' decsvd <- bdSVD_lapack( A, bscale = FALSE, bcenter = FALSE ) # No matrix normalization
//' decsvd$d
//' decsvd$u
//' 
//' # svd with normalization
//' decvsd <- bdSVD_lapack( A, bscale = TRUE, bcenter = TRUE) # Matrix normalization
//' decvsd <- bdSVD_lapack( A ) # Matrix normalization too
//' decsvd$d
//' decsvd$u
//' 
//' # svd with scaled matrix (sd)
//' decvsd <- bdSVD_lapack( A, bscale = TRUE, bcenter = FALSE) # Scaled matrix
//' 
//' decsvd$d
//' decsvd$u
//' # svd with centered matrix (sd)
//' decvsd <- bdSVD_lapack( A, bscale = FALSE, bcenter = TRUE) # Centered matrix
//' decsvd$d
//' decsvd$u
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSVD_lapack ( Rcpp::RObject x, Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true)
{
  auto dmtype = beachmat::find_sexp_type(x);
  bool bcent, bscal;
  
  if(bcenter.isNull())  bcent = true ;
  else    bcent = Rcpp::as<bool>(bcenter);
  
  if(bscale.isNull())  bscal = true ;
  else    bscal = Rcpp::as<bool>(bscale);
  

  Eigen::MatrixXd X;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( x.isS4() == true){
      X = read_DelayedArray(x);
    }else {
      try{
        X = Rcpp::as<Eigen::MatrixXd >(x);
      }catch(std::exception &ex) {
        X = Rcpp::as<Eigen::VectorXd >(x);
      }
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  svdeig retsvd =  RcppbdSVD_lapack( X, bcent, bscal);
  
  return List::create(Named("d") = retsvd.d,
                      Named("u") = retsvd.u,
                      Named("v") = retsvd.v);
}




/***R

library(microbenchmark)
library(DelayedArray)
library(BigDataStatMeth)
library(rhdf5)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/BigDataStatMeth/tmp")

# Rows --> Individuals  (small)
# Cols --> Variables (SNP's or ...) (very big)

dades <- BigDataStatMeth::bdSVD_hdf5("tmp_blockmult.hdf5", group = "INPUT", dataset = "A", k=4, q=1) 

dades$file
fprova <- H5Fopen("tmp_blockmult.hdf5")
fprova
fprova$INPUT$A[1:25,1:25]
csvd <- bdSVD_lapack(fprova$INPUT$A[1:25,1:25], bcenter = FALSE, bscale = FALSE)

csvd3 <- bdSVD_lapack(fprova$INPUT$A[1:25,1:25], bcenter = FALSE, bscale = FALSE)
csvd3$u[1:5,1:5]
csvd3$v
csvd$d
csvd$u
csvd$v

csvd2 <- svd(fprova$INPUT$A[1:25,1:25])
csvd2$v[1:5,1:5]



csvd$u %*% diag(csvd$d)

(csvd$u %*% diag(csvd$d))[1:5,1:5]
t(fprova$tmpgroup$A0)[1:5,1:5]


svd( scale(fprova$INPUT$A))$
bdSVD_lapack(fprova$INPUT$A)$d

h5closeAll()
# 



*/