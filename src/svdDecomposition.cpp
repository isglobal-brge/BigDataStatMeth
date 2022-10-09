#include "include/svdDecomposition.h"

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
        if(bcenter ==true || bscale == true)  {
            nX = RcppNormalize_Data(X, bcenter, bscale);
            Xtcp =  bdtcrossproduct(nX);
        } else {
            Xtcp =  bdtcrossproduct(X);
        }
        
        Spectra::DenseSymMatProd<double> op(Xtcp);
        Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, ncv);
        
        // Initialize and compute
        eigs.init();
        nconv = eigs.compute();
        
        if(eigs.info() == Spectra::SUCCESSFUL) {
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
        if(bcenter ==true || bscale==true ) {
            Xcp =  bdcrossproduct(nX);  
        } else {
            Xcp =  bdcrossproduct(X);
        }  
        
        Spectra::DenseSymMatProd<double> opv(Xcp);
        Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigsv(&opv, k, ncv);
        
        // Initialize and compute
        eigsv.init();
        nconv = eigsv.compute();
        
        // Retrieve results
        if(eigsv.info() == Spectra::SUCCESSFUL) {
            retsvd.v = eigsv.eigenvectors();
        } else {
            retsvd.bokd = false;
        }
    }
    
    return retsvd;
}



// Lapack SVD decomposition - Optimized algorithm with dgesdd
svdeig RcppbdSVD_lapack( Eigen::MatrixXd& X, bool bcenter, bool bscale, bool complete )
{
    
    
    svdeig retsvd;
    
    char Schar='S';
    char Achar='A';
    int info = 0;
    
    try {
    
        if(bcenter ==true || bscale == true) {
            X = RcppNormalize_Data(X, bcenter, bscale);
        }
        
        int m = X.rows();
        int n = X.cols();
        int lda = std::max(1,m);
        int ldu = std::max(1,m);
        int ldvt = std::min(m, n);
        int k = std::min(m,n);
        int lwork;
        
        if( complete == false ) {
            lwork = 4*std::min(m,n)*std::min(m,n) + 7*std::min(m,n);
        } else {
            lwork = 4*std::min(m,n)*std::min(m,n) + 6*std::min(m,n) + std::max(m,n);
        }
        
        Eigen::VectorXd s = Eigen::VectorXd::Zero(k);
        Eigen::VectorXd work = Eigen::VectorXd::Zero(lwork);
        Eigen::VectorXi iwork(8*std::min(m,n));
        Eigen::MatrixXd u;
        Eigen::MatrixXd vt = Eigen::MatrixXd::Zero(ldvt,n);
        
        if( complete == false ) {
            u = Eigen::MatrixXd::Zero(ldu,k);
            dgesdd_( &Schar, &m, &n, X.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);
        } else {
            u = Eigen::MatrixXd::Zero(ldu,m);
            dgesdd_( &Achar, &m, &n, X.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);
        }
        
        retsvd.d = s;
        retsvd.u = u;
        retsvd.v = vt.transpose();
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< "C++ exception RcppbdSVD_lapack : "<< ex.what();
        return retsvd;
    } catch (...) {
        ::Rf_error("C++ exception RcppbdSVD_lapack (unknown reason)");
        return retsvd;
    } 
    
    return retsvd;
}




// ##' @param k number of local SVDs to concatenate at each level 
// ##' @param q number of levels
svdeig RcppbdSVD_hdf5_Block( H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
                             int irows, int icols, double dthreshold, Rcpp::Nullable<int> threads = R_NilValue )
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
    
    First_level_SvdBlock_decomposition_hdf5( file, dataset, k, q, nev, bcenter, bscale, irows, icols, dthreshold, threads);
    
    for(int j = 1; j < q; j++) { // For each decomposition level :
        Next_level_SvdBlock_decomposition_hdf5(file, strGroupName, k, j, bcenter, bscale, dthreshold, threads);
    }
    
    // Get dataset names
    StringVector joindata =  get_dataset_names_from_group(file, strGroupName, strPrefix);

    // 1.- Join matrix and remove parts from file
    std::string strnewdataset = std::string((joindata[0])).substr(0,1);
    join_datasets(file, strGroupName, joindata, strnewdataset);
    remove_HDF5_multiple_elements_ptr(file, strGroupName, joindata);
    
    // 2.- Get SVD from Blocks full mattrix
    DataSet datasetlast = file->openDataSet(strGroupName + "/" + strnewdataset);
    IntegerVector dims_out = get_HDF5_dataset_size(datasetlast);

    Eigen::MatrixXd matlast;
    matlast = GetCurrentBlock_hdf5(file, &datasetlast, 0, 0, dims_out[0],dims_out[1]);
    

    if( nev < std::min( int(matlast.rows()), int(matlast.cols()) )) {
        retsvd = RcppbdSVD(matlast, nev, 0, false, false);
    } else {
        retsvd = RcppbdSVD_lapack(matlast, false, false, false);
    }
    

    // Write results to hdf5 file : in folder "SVD" and dataset "SVD".<name input dataset>
    // Create structure and write d 
    StringVector name = get_dataset_names_from_dataset_ptr(dataset);

    create_HDF5_groups_ptr(file,"SVD/"+ name[0]);
    write_HDF5_matrix_ptr(file, "SVD/"+ name[0]+"/d", wrap(retsvd.d));

    Eigen::MatrixXd v;    

    
    // 3.- crossprod initial matrix and svdA$u
/***    DataSet* normalizedData = new DataSet(file->openDataSet(strGroupName + "/normalmatrix"));
    IntegerVector dims_out_normal = get_HDF5_dataset_size(*normalizedData);
***/ 
    
    Eigen::MatrixXd A;
    

    if( bcenter == true || bscale == true) {
        DataSet* normalizedData = nullptr;
        normalizedData = new DataSet(file->openDataSet(strGroupName + "/normalmatrix"));
        IntegerVector dims_out_normal = get_HDF5_dataset_size(*normalizedData);
        A = GetCurrentBlock_hdf5_Original(file, normalizedData, 0, 0,dims_out_normal[0], dims_out_normal[1] );
        
        normalizedData->close();
        delete(normalizedData);
    } else {
        
        IntegerVector dims_out_normal = get_HDF5_dataset_size(*dataset);
        A = GetCurrentBlock_hdf5(file, dataset, 0, 0,dims_out_normal[0], dims_out_normal[1] );
    }
    
    // Rcpp::Rcout<<"\n Podria ser aquí ?? ....";
    // IntegerVector dims_out_normal = get_HDF5_dataset_size(*normalizedData);
    // 
    // Eigen::MatrixXd A = GetCurrentBlock_hdf5_Original(file, normalizedData, 0, 0,dims_out_normal[0], dims_out_normal[1] );
    // normalizedData->close();
    
    v = Bblock_matrix_mul_parallel(A, retsvd.u, 1024, threads); //  PARALLEL ==> NOT PARALLEL

    // 4.- resuls / svdA$d
    v = v.array().rowwise()/(retsvd.d).transpose().array();

    // 5.- Get u and v : 
    //        --> If trans : u=v i v=u  
    
    // Simplified in next lines (24/02/2022)
    //
    // if (transp == true)  {
    //     retsvd.v = retsvd.u;
    //     retsvd.u = v;
    // } else {
    //     retsvd.v = v;
    // }
    // write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/u", wrap(retsvd.u));
    // write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/v", wrap(retsvd.v));
    //

    
    
    if (transp == true)  {
        write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/u", wrap(retsvd.v));
        write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/v", wrap(retsvd.u));
    } else {
        write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/u", wrap(retsvd.u));
        write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/v", wrap(retsvd.v));
    }

    
    // Clean data
    if( bcenter == true || bscale == true) {
        remove_HDF5_multiple_elements_ptr(file, strGroupName, "normalmatrix");
        remove_HDF5_element_ptr(file, strGroupName);
    }
  
  
  } catch(FileIException& error) { // catch failure caused by the H5File operations
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_Block (File IException)" );
        return retsvd;
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_Block (DataSet IException)" );
        return retsvd;
  } catch(GroupIException& error) { // catch failure caused by the Group operations
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_Block (Group IException)" );
        return retsvd;
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_Block (DataSpace IException)" );
        return retsvd;
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_Block (Data TypeIException)" );
        return retsvd;  
  } catch(std::exception &ex) {
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        Rcpp::Rcout<< "C++ exception RcppbdSVD_hdf5_Block : "<< ex.what();
        return retsvd;  
  } catch (...) {
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        ::Rf_error("C++ exception RcppbdSVD_hdf5_Block (unknown reason)");
        return retsvd;
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
                       int k, int q, int nev, bool bcenter, bool bscale, double dthreshold, 
                       Rcpp::Nullable<int> ithreads = R_NilValue )
{
    
    svdeig retsvd;
    Eigen::MatrixXd X;
    
    H5File* file = nullptr;
    DataSet* dataset = nullptr;
    
    try {
    
        // Open an existing file and dataset.
        file = new H5File( filename, H5F_ACC_RDWR );
        
        if(exists_HDF5_element_ptr(file, strsubgroup + "/" + strdataset)){
            dataset = new DataSet(file->openDataSet(strsubgroup + "/" + strdataset));
        } else {
            file->close();
            delete(file);
            throw std::range_error("Dataset not exits"); 
        }
        
        // Remove previous results
        if(exists_HDF5_element_ptr(file,"SVD/"+strdataset)) {
            Rcpp::Rcout<<"SVD - old dataset have been REMOVED \n";
            remove_HDF5_element_ptr(file,"SVD/"+strdataset);
        }
        
        // Get dataset dims
        IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
        
        hsize_t offset[2] = {0,0};
        hsize_t count[2] = { (unsigned long long)dims_out[0], (unsigned long long)dims_out[1]};

        // Small matrices ==> Direct SVD
        if( dims_out[0] < MAXSVDBLOCK &&  dims_out[1] < MAXSVDBLOCK ) {
            
            X = GetCurrentBlock_hdf5_Original(file, dataset, offset[0], offset[1], count[0], count[1]); 
            
            
            // retsvd = RcppbdSVD(X, k, nev, bcenter, bscale);
            retsvd = RcppbdSVD_lapack(X, bcenter, bscale, false);
            
            create_HDF5_groups_ptr(file,"SVD/"+ strdataset);
            
            // Create filegroup
            write_HDF5_matrix_transposed_ptr(file, "SVD/"+ strdataset+"/u", wrap(retsvd.u));
            write_HDF5_matrix_transposed_ptr(file, "SVD/"+ strdataset+"/v", wrap(retsvd.v));
            write_HDF5_matrix_transposed_ptr(file, "SVD/"+ strdataset+"/d", wrap(retsvd.d));
        
        } else {
        
            // data stored transposed in hdf5
            int xdim = (unsigned long long)dims_out[1];
            int ydim = (unsigned long long)dims_out[0];
            
            retsvd = RcppbdSVD_hdf5_Block( file, dataset, k, q, nev, bcenter, bscale, xdim, ydim, dthreshold, wrap(ithreads));
            
        }
        
    }  catch( FileIException& error ) { // catch failure caused by the H5File operations
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5 (File IException)" );
        return retsvd;
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5 (DataSet IException)" );
        return retsvd;
    } catch(std::exception &ex) {
        dataset->close();
        file->close();
        delete(dataset);
        delete(file);
        Rcpp::Rcout<< ex.what();
        return retsvd;
    }
    
    dataset->close();
    file->close();
    
    delete(dataset);
    delete(file);
    
    return retsvd;
   
}



svdeig RcppbdSVD_hdf5_ptr( H5File* file, std::string strsubgroup, std::string strdataset,  
                       int k, int q, int nev, bool bcenter, bool bscale, bool bstorehdf5, 
                       double dthreshold, Rcpp::Nullable<int> ithreads = R_NilValue)
{
  
  svdeig retsvd;
  Eigen::MatrixXd X;
  std::string strdatasetFull;
  
  DataSet* dataset = nullptr;
  
  
  try
  {
      
      strdatasetFull = strsubgroup + "/" + strdataset;
    // Open an existing file and dataset.
    if(exists_HDF5_element_ptr(file, strsubgroup + "/" + strdataset)){
      //..// dataset = file->openDataSet(strsubgroup + "/" + strdataset);
        dataset = new DataSet(file->openDataSet(strsubgroup + "/" + strdataset));
    }
    else {
        file->close();
        delete(file);
        throw std::range_error("Dataset not exits"); 
    }

    // Get dataset dims
    IntegerVector dims_out = get_HDF5_dataset_size_ptr(dataset);
    
    hsize_t offset[2] = {0,0};
    hsize_t count[2] = { (unsigned long long)dims_out[0], (unsigned long long)dims_out[1]};
    
    // In memory computation for small matrices (rows or columns<5000)
    // Block decomposition for big mattrix
    if( std::max(dims_out[0], dims_out[1]) < MAXSVDBLOCK && bstorehdf5 == false )
    {

      X = GetCurrentBlock_hdf5( file, dataset, offset[0], offset[1], count[0], count[1]);
      X.transposeInPlace();
      
      retsvd = RcppbdSVD(X, k, nev, bcenter, bscale);
      
    } else {

      // data stored transposed in hdf5
      int xdim = (unsigned long long)dims_out[1];
      int ydim = (unsigned long long)dims_out[0];

      // Remove previous results
      if(exists_HDF5_element_ptr(file,"SVD/"+strdataset))  {
          remove_HDF5_element_ptr(file,"SVD/"+strdataset);
      }
      
      retsvd = RcppbdSVD_hdf5_Block( file, dataset, k, q, nev, bcenter, bscale, xdim, ydim, dthreshold, wrap(ithreads) );

    }

  } catch(FileIException& error) { // catch failure caused by the H5File operations
      dataset->close();
      delete(dataset);
        file->close();
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_ptr (File IException)" );
  } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
      dataset->close();
      delete(dataset);
        file->close();
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_ptr (DataSet IException)" );
  } catch(GroupIException& error) { // catch failure caused by the Group operations
      dataset->close();
      delete(dataset);
        file->close();
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_ptr (Group IException)" );
  } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
      dataset->close();
      delete(dataset);
        file->close();
        
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_ptr (DataSpace IException)" );
  } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
      dataset->close();
      delete(dataset);
        file->close();
        delete(file);
        ::Rf_error( "c++ exception RcppbdSVD_hdf5_ptr (Data TypeIException)" );
  } catch(std::exception &ex) {
      dataset->close();
      delete(dataset);
      file->close();
      delete(file);
        Rcpp::Rcout<<"c++ exception in RcppbdSVD_hdf5_ptr : "<< ex.what();
  } catch (...) {
      dataset->close();
      delete(dataset);
        file->close();
        delete(file);
        ::Rf_error("C++ exception RcppbdSVD_hdf5_ptr (unknown reason)");
  }
  
  
  dataset->close();
  delete(dataset);
  
  return retsvd;
  
}










svdeig RcppCholDec(const Eigen::MatrixXd& X)
{
  Eigen::MatrixXd mX = X;
  svdeig decomp;
  
  Eigen::LDLT<Eigen::MatrixXd> cholSolv = mX.ldlt();
  
  Eigen::LLT<Eigen::MatrixXd> lltOfA(X); // compute the Cholesky decomposition of A
  
  if(lltOfA.info() == Eigen::NumericalIssue)  {
    Rcpp::Rcout<<"Possibly non semi-positive definitie matrix!. Matrix returned as 0";
    decomp.v = Eigen::MatrixXd::Zero(2,2);
    decomp.u = Eigen::MatrixXd::Zero(2,2);
    return(decomp);
  //..//  throw std::runtime_error("Possibly non semi-positive definitie matrix!");
    
  } 

  // if symetric + positive definite -> info = Success
  // else no Cholesky Decomposition --> svd decomposition with Spectra
  if( cholSolv.info()==Eigen::Success )
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




//' Inverse Cholesky
//' 
//' This function get the inverse of a numerical matrix. If x is hermitian and positive-definite matrix then gets the inverse using Cholesky decomposition
//' 
//' 
//' @param X numerical matrix. If x is Hermitian and positive-definite performs
//' @return inverse matrix of d 
//' @examples
//' 
//' 
//' A <- matrix(c(3,4,3,4,8,6,3,6,9), byrow = TRUE, ncol = 3)
//' bdInvCholesky(A)
//' 
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd bdInvCholesky (const Rcpp::RObject & X )
{
  
  svdeig result;
  Eigen::MatrixXd mX;
  

    try{  
        
        if ( TYPEOF(X) == INTSXP ) {
            mX = Rcpp::as<Eigen::MatrixXi>(X).cast<double>();
        } else {
            mX = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
        }
    }
    catch(std::exception &ex) { }

  result = RcppCholDec(mX);
  
  return (result.v);
  
}



//' k first SVD components for DelayedArray 
//' 
//' This function gets k first components from svd decomposition of numerical matrix
//' 
//' @param X numerical matrix
//' @param k number of eigen values , this should satisfy k = min(n, m) - 1
//' @param nev (optional, default nev = n-1) Number of eigenvalues requested. 
//' This should satisfy 1<= nev <= n, where n is the size of matrix. 
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering 
//' is done by subtracting the column means (omitting NAs) of x from their 
//' corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is 
//' done by dividing the (centered) columns of x by their standard deviations if 
//' center is TRUE, and the root mean square otherwise. If scale is FALSE, no 
//' scaling is done.
//' @return u eigenvectors of AA^t, mxn and column orthogonal matrix
//' @return v eigenvectors of A^tA, nxn orthogonal matrix
//' @return d singular values, nxn diagonal matrix (non-negative real values)
//' @examples
//' 
//' n <- 500
//' A <- matrix(rnorm(n*n), nrow=n, ncol=n)
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
Rcpp::RObject bdSVD (const Rcpp::RObject & X, Rcpp::Nullable<int> k=0, Rcpp::Nullable<int> nev=0,
                     Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true)
{
  
  auto dmtype = beachmat::find_sexp_type(X);
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
  Eigen::MatrixXd mX;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( X.isS4() == true){
      // mX = read_DelayedArray(X);
      throw("Only numeric matrix allowd");
    }else {
      try{
        mX = Rcpp::as<Eigen::MatrixXd >(X);
      }catch(std::exception &ex) {
        mX = Rcpp::as<Eigen::VectorXd >(X);
      }
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  svdeig retsvd;
  retsvd = RcppbdSVD(mX,ks,nvs, bcent, bscal);
  
  ret["u"] = retsvd.u;
  ret["v"] = retsvd.v;
  ret["d"] = retsvd.d;
  
  return Rcpp::wrap(ret);
  
}




//' Block SVD decomposition for hdf5 files using an incremental algorithm.
//'
//' Singular values and left singular vectors of a real nxp matrix 
//' Block SVD decomposition using an incremental algorithm.
//' @param file a real nxp matrix in hdf5 file
//' @param group group in hdf5 data file where dataset is located
//' @param dataset matrix dataset with data to perform SVD
//' @param k number of local SVDs to concatenate at each level 
//' @param q number of levels
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
//' @param rankthreshold double, threshold used to determine the range of the array. The matrix rank is equal to the number of
//'  singular values different from the threshold. By default, threshold = 0 is used to get the matrix rank , but it can be
//'  changed to an approximation of 0.
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
Rcpp::RObject bdSVD_hdf5 (const Rcpp::RObject file, Rcpp::Nullable<CharacterVector> group = R_NilValue, 
                          Rcpp::Nullable<CharacterVector> dataset = R_NilValue,
                          Rcpp::Nullable<int> k=2, Rcpp::Nullable<int> q=1,
                          Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true,
                          Rcpp::Nullable<double> rankthreshold = 0.0,
                          Rcpp::Nullable<int> threads = R_NilValue)
{
  
  std::string filename;
  svdeig retsvd;
  double dthreshold;
  
  try {
    
    int ks, qs, nvs = 0;
    bool bcent, bscal;
    CharacterVector strgroup, strdataset;

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
    
    if(is<CharacterVector>(file)) {
      filename = as<std::string>(file);
    } else {
        Rcpp::Rcout<< "File name must be character string";
        return List::create(Named("file") = "");
    }
    
    if(rankthreshold.isNull()) {  
        dthreshold = 0 ;
    } else {
        if( Rcpp::as<double>(rankthreshold) > 0.1 ) {
            Rcpp::Rcout<< "Threshold to big, please set threshold with value lower than 0.1";
            return List::create(Named("file") = filename);
        } else if( Rcpp::as<double>(rankthreshold) < 0 ) {
            Rcpp::Rcout<< "Threshold must be a positive value near zero";
            return List::create(Named("file") = filename);
        } else {
            dthreshold = Rcpp::as<double>(rankthreshold);
        }
    }

      
    retsvd = RcppbdSVD_hdf5( filename, as<std::string>(strgroup), as<std::string>(strdataset), ks, qs, nvs, bcent, bscal, dthreshold, threads );
    
    
  } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return List::create(Named("file") = R_NilValue);
  }

  return List::create(Named("file") = filename);
  
}


/***

//' Complete SVD with Lapack Functions for  RObjects
//' 
//' This function performs a complete svd decomposition of numerical matrix
//' 
//' @param X numerical or Delayed Array matrix
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
//' @param complete (optional, defalut = FALSE) . If complete is TRUE svd function returns complete u and v
//' @return u eigenvectors of AA^t, mxn and column orthogonal matrix
//' @return v eigenvectors of A^tA, nxn orthogonal matrix
//' @return d singular values, nxn diagonal matrix (non-negative real values)
//' @examples
//' 
//' library(BigDataStatMeth)
//' n <- 500
//' A <- matrix(rnorm(n*n), nrow=n, ncol=n)
//' 
//' # svd without normalization
//' decsvd <- bdSVD_lapack_not_optim( A, bscale = FALSE, bcenter = FALSE ) # No matrix normalization
//' decsvd$d
//' decsvd$u
//' 
//' # svd with normalization
//' decvsd <- bdSVD_lapack_not_optim( A, bscale = TRUE, bcenter = TRUE) # Matrix normalization
//' decvsd <- bdSVD_lapack_not_optim( A ) # Matrix normalization too
//' decsvd$d
//' decsvd$u
//' 
//' # svd with scaled matrix (sd)
//' decvsd <- bdSVD_lapack_not_optim( A, bscale = TRUE, bcenter = FALSE) # Scaled matrix
//' @export
// [[Rcpp::export]]

Rcpp::RObject bdSVD_lapack_not_optim ( Rcpp::RObject X, Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true,  Rcpp::Nullable<bool> complete=false )
{
  auto dmtype = beachmat::find_sexp_type(X);
  bool bcent, bscal, bcomp;
  
  if(bcenter.isNull())  bcent = true ;
  else    bcent = Rcpp::as<bool>(bcenter);
  
  if(bscale.isNull())  bscal = true ;
  else    bscal = Rcpp::as<bool>(bscale);
  
  if(complete.isNull())  bcomp = false;
  else    bcomp = Rcpp::as<bool>(complete);
  

  Eigen::MatrixXd mX;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( X.isS4() == true){
      mX = read_DelayedArray(X);
    }else {
      try{
        mX = Rcpp::as<Eigen::MatrixXd >(X);
      }catch(std::exception &ex) {
        mX = Rcpp::as<Eigen::VectorXd >(X);
      }
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  svdeig retsvd =  RcppbdSVD_lapack( mX, bcent, bscal, bcomp);
  
  return List::create(Named("d") = retsvd.d,
                      Named("u") = retsvd.u,
                      Named("v") = retsvd.v);
}

***/





//' Complete SVD with Lapack Functions for RObjects
//' 
//' This function performs a complete svd decomposition of numerical matrix
//' 
//' @param X numerical matrix
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
//' @param complete (optional, defalut = FALSE) . If complete is TRUE svd function returns complete u and v
//' @return u eigenvectors of AA^t, mxn and column orthogonal matrix
//' @return v eigenvectors of A^tA, nxn orthogonal matrix
//' @return d singular values, nxn diagonal matrix (non-negative real values)
//' @examples
//' 
//' library(BigDataStatMeth)
//' n <- 500
//' A <- matrix(rnorm(n*n), nrow=n, ncol=n)
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
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSVD_lapack ( Rcpp::RObject X, Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true,  Rcpp::Nullable<bool> complete=false )
{
    auto dmtype = beachmat::find_sexp_type(X);
    bool bcent, bscal, bcomp;
    
    if(bcenter.isNull())  bcent = true ;
    else    bcent = Rcpp::as<bool>(bcenter);
    
    if(bscale.isNull())  bscal = true ;
    else    bscal = Rcpp::as<bool>(bscale);
    
    if(complete.isNull())  bcomp = false;
    else    bcomp = Rcpp::as<bool>(complete);
    
    
    Eigen::MatrixXd mX;
    Rcpp::List ret;
    
    if ( dmtype == INTSXP || dmtype==REALSXP ) {
        if ( X.isS4() == true){
            // mX = read_DelayedArray(X);
            throw("Only numeric matrix allowd");
        }else {
            try{
                mX = Rcpp::as<Eigen::MatrixXd >(X);
            }catch(std::exception &ex) {
                mX = Rcpp::as<Eigen::VectorXd >(X);
            }
        }
        
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    
    svdeig retsvd =  RcppbdSVD_lapack( mX, bcent, bscal, bcomp);
    
    return List::create(Named("d") = retsvd.d,
                        Named("u") = retsvd.u,
                        Named("v") = retsvd.v);
}

