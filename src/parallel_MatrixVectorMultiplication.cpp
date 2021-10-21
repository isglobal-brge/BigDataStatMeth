#include "include/parallel_MatrixVectorMultiplication.h"
#include "include/ReadDelayedData.h"




struct MatrixMult : public RcppParallel::Worker  {
  
  // input matrix to read from
  RcppParallel::RMatrix<double> mat;
  RcppParallel::RMatrix<double> tmat;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // other variables
  std::size_t numcol;
  
  
  // Constructor
  MatrixMult( Rcpp::NumericMatrix mat,  Rcpp::NumericMatrix tmat, Rcpp::NumericMatrix rmat, std::size_t numcol)
    : mat(mat), tmat(tmat), rmat(rmat), numcol(numcol)  {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    
    for (std::size_t i = begin; i < end; i++) 
    {
      RcppParallel::RMatrix<double>::Row rowmat = mat.row(i); // rows we will operate on
      
      for (std::size_t j = 0; j < numcol ; j++) 
      {
        RcppParallel::RMatrix<double>::Column coltmat = tmat.column(j); // cols we will operate on 
        size_t selem = coltmat.length();
        
        for(std::size_t k = 0; k<selem; k++) 
        {
          rmat(i,j) =  rmat(i,j) + (rowmat[k] * coltmat[k]);
          //if(i!=j)
          //  rmat(j,i) =  rmat(i,j);
        }
      }
    }
  }
};



struct MatrixMultEig : public RcppParallel::Worker  {
  
  // input matrix to read from
  const Eigen::Map<Eigen::MatrixXd> mat;
  const Eigen::Map<Eigen::MatrixXd> tmat;
  // Eigen::MatrixXd tmat;
  
  // output matrix to write to
  Eigen::Map<Eigen::MatrixXd> rmat;
  
  // other variables
  std::size_t numcol;
  
  // Constructor
  MatrixMultEig( const Eigen::Map<Eigen::MatrixXd>* mat,  const Eigen::Map<Eigen::MatrixXd>* tmat, Eigen::Map<Eigen::MatrixXd>* rmat, std::size_t numcol)
    : mat(*mat), tmat(*tmat), rmat(*rmat), numcol(numcol)  {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {

    for (std::size_t i = begin; i < end; i++) 
    {
      Eigen::VectorXd rowmat = mat.row(i);
      
      for (std::size_t j = 0; j < numcol ; j++) 
      {
        ColumnVector coltmat = tmat.col(j).transpose();
        size_t selem = tmat.col(i).size() ;
        
        for(std::size_t k = 0; k<selem; k++) 
        {
          rmat(i,j) =  rmat(i,j) + (rowmat[k] * coltmat[k]);

          //if(i!=j)
          //  rmat(j,i) =  rmat(i,j);
        }
      }
    }
  }
};



struct MatrixVectMultEig: public RcppParallel::Worker {
  
  // input matrix to read from
  const Eigen::Map<Eigen::MatrixXd> mat;
  // input vector to read from
  const Eigen::Map<Eigen::VectorXd> vect;
  
  // output matrix to write to
  Eigen::Map<Eigen::VectorXd> rmat;
  
  // other variables
  std::size_t numcol;
  
  // Constructor
  MatrixVectMultEig( const Eigen::Map<Eigen::MatrixXd>* mat,  const Eigen::Map<Eigen::VectorXd>* vect, Eigen::Map<Eigen::VectorXd>* rmat, std::size_t numcol)
    : mat(*mat), vect(*vect), rmat(*rmat), numcol(numcol)  {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    
    for (std::size_t i = begin; i < end; i++) 
    {
      //RcppParallel::RMatrix<double>::Row rowmat = mat.row(i); // rows we will operate on
      Eigen::VectorXd rowmat = mat.row(i);
      ColumnVector coltmat = vect.transpose();
      
      for (std::size_t j = 0; j < numcol ; j++) 
      {
        rmat[i] =  rmat[i] + (rowmat[j] * vect[j]);
      }
    }
  }
};


struct MatrixVectMult: public RcppParallel::Worker {
  
  // input matrix to read from
  RcppParallel::RMatrix<double> mat;
  
  // input vector to read from
  RcppParallel::RVector<double> vect;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // Constructor
  MatrixVectMult( Rcpp::NumericMatrix mat,  Rcpp::NumericVector vect, Rcpp::NumericMatrix rmat)
    : mat(mat), vect(vect), rmat(rmat)  {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    
    for (std::size_t i = begin; i < end; i++) 
    {
      RcppParallel::RMatrix<double>::Row rowmat = mat.row(i); // rows we will operate on
      size_t numcol = rowmat.length();
      
      for (std::size_t j = 0; j < numcol ; j++) 
      {
          rmat(i,0) =  rmat(i,0) + (rowmat[j] * vect[j]);

      }
    }
  }
  
};



Rcpp::NumericMatrix rcpp_parallel_Xy(Rcpp::NumericMatrix mat, Rcpp::NumericVector y) {
  
  // allocate the matrix we will return and the matrix with provisional results
  Rcpp::NumericMatrix rmat(mat.nrow(), 1);
  
  try
  {
    if(mat.ncol()==y.size() )
    {
      try{
        // Get Xy
        MatrixVectMult matvectmult( mat, y, rmat);
        RcppParallel::parallelFor(0, mat.nrow(), matvectmult);
        
      }catch(std::exception &ex) {	
        forward_exception_to_r(ex);
      }
      
    } else
    {
      throw std::range_error("non-conformable arguments");
    }
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }  
  
  
  return rmat;
}



// Eigen::VectorXd rcpp_parallel_Xy_eigen( Eigen::MatrixXd mat, Eigen::VectorXd y) {
void rcpp_parallel_Xy_eigen( Eigen::Map<Eigen::MatrixXd>* mat, const Eigen::Map<Eigen::VectorXd>* y, Eigen::Map<Eigen::VectorXd>* rmat) {
  
  // allocate the matrix we will return and the matrix with provisional results
  //..//Rcpp::NumericMatrix rmat(mat.rows(), 1);
  // Rcpp::NumericVector rmat(mat->rows());
  
  try
  {
    if(mat->cols()==y->size() )
    {
      try{
        
        // Get Xy
        MatrixVectMultEig matvectmult( mat, y, rmat, y->size());
        RcppParallel::parallelFor(0, mat->rows(), matvectmult);
        
      }catch(std::exception &ex) {	
        forward_exception_to_r(ex);
      }
      
    } else
    {
      throw std::range_error("non-conformable arguments");
    }
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }  
  
  // return (Rcpp::as<Eigen::VectorXd>(rmat));
}



Rcpp::NumericMatrix rcpp_xwxt(Rcpp::NumericMatrix mat, Rcpp::NumericVector w) {
  
  // Diagonalize weights
  Rcpp::NumericMatrix wd(Rcpp::diag(w));
  
  // Transpose de matrix
  Rcpp::NumericMatrix tmat(Rcpp::transpose(mat));
  
  // allocate the matrix we will return and the matrix with provisional results
  Rcpp::NumericMatrix rprov(mat.nrow(), mat.ncol());
  Rcpp::NumericMatrix rmat(mat.nrow(), mat.nrow());
  
  try
  {
    if(mat.ncol()==wd.nrow())
    {
      try{
        // Provisional results xw
        MatrixMult matmultp( mat, wd, rprov, wd.ncol());
        RcppParallel::parallelFor(0, mat.nrow(), matmultp);
        
        // Get xwxt
        MatrixMult matmult( rprov, tmat, rmat, tmat.ncol());
        RcppParallel::parallelFor(0, mat.nrow(), matmult);
        
        
      }catch(std::exception &ex) {	
        forward_exception_to_r(ex);
      }
      
    } else
    {
      throw std::range_error("non-conformable arguments");
    }
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }  
  
  
  return rmat;
}


// Eigen::Map<Eigen::MatrixXd> rcpp_xwxt_eig(const Eigen::Map<Eigen::MatrixXd>* mat, const Eigen::Map<Eigen::VectorXd>* w, Eigen::Map<Eigen::MatrixXd>* rmat) {
void rcpp_xwxt_eig(const Eigen::Map<Eigen::MatrixXd>* mat, 
                   const Eigen::Map<Eigen::VectorXd>* w, Eigen::Map<Eigen::MatrixXd>* rmat) {
  
  // Diagonalize weights
  
  Eigen::MatrixXd wdd = w->asDiagonal();
  Eigen::Map<Eigen::MatrixXd> wd(wdd.data(), w->size(), w->size());
  
  // Transpose de matrix
  Eigen::MatrixXd tmatd(mat->transpose());
  Eigen::Map<Eigen::MatrixXd> tmat(tmatd.data(), mat->cols(), mat->rows());
  
  // allocate the matrix we will return and the matrix with provisional results
  Eigen::MatrixXd rprovd = Eigen::MatrixXd::Zero(mat->rows(), mat->cols());
  Eigen::Map<Eigen::MatrixXd> rprov(rprovd.data(), mat->rows(), mat->cols());
  
  
  try
  {
    if(mat->cols()==wd.rows())
    {
      try{
        // Provisional results xw
        MatrixMultEig matmultp( mat, &wd, &rprov, wd.cols());
        RcppParallel::parallelFor(0, mat->rows(), matmultp);
        
        // Get xwxt
        // MatrixMult matmult( rprov, Rcpp::wrap(tmat), rmat, tmat.cols());
        MatrixMultEig matmultp2( &rprov, &tmat, rmat, tmat.cols());
        RcppParallel::parallelFor(0, rprov.rows(), matmultp2);
        // Rcpp::Rcout<<"Matriu ginv: \n"<< tmatd;
        
      }catch(std::exception &ex) {	
        forward_exception_to_r(ex);
      }
      
    } else
    {
      throw std::range_error("non-conformable arguments");
    }
    
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }  
  
  // return (Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(rmat));
}











Rcpp::NumericMatrix rcpp_xtwx(Rcpp::NumericMatrix mat, Rcpp::NumericVector w) {
  
  // Diagonalize weights
  Rcpp::NumericMatrix wd(Rcpp::diag(w));
  
  // Transpose de matrix
  Rcpp::NumericMatrix tmat(Rcpp::transpose(mat));
  
  // allocate the matrix we will return and the matrix with provisional results
  Rcpp::NumericMatrix rprov(tmat.nrow(), tmat.ncol());
  Rcpp::NumericMatrix rmat(tmat.nrow(), tmat.nrow());
  
  try
  {
    if(tmat.ncol()==wd.nrow())
    {
      try{
        // Provisional results xw
        MatrixMult matmultp( tmat, wd, rprov, wd.ncol());
        RcppParallel::parallelFor(0, tmat.nrow(), matmultp);
        
        // Get xwxt
        MatrixMult matmult( rprov, mat, rmat, mat.ncol());
        RcppParallel::parallelFor(0, tmat.nrow(), matmult);
        
        
      }catch(std::exception &ex) {	
        forward_exception_to_r(ex);
      }
      
    } else
    {
      throw std::range_error("non-conformable arguments");
    }
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
  
  return rmat;
}




/** CRIDES R **/

// [[Rcpp::export]]
Rcpp::RObject parxwxt(Rcpp::RObject X, Rcpp::RObject W)
{
  auto dmtype = beachmat::find_sexp_type(X);
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  Rcpp::List ret;
  Rcpp::NumericVector w = Rcpp::as<Rcpp::NumericVector >(W);
  
  if ( dmtype == INTSXP ) {
    if ( X.isS4() == true){
      eX = read_DelayedArray_int( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }
    
    Rcpp::NumericMatrix XX = rcpp_xwxt( Rcpp::wrap(eX), w);
    
    if ( X.isS4() == true)
    {

      ncols = XX.cols();
      nrows = XX.rows();
      
      //beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);
      beachmat::output_param oparam(X);
      auto out_dmat = beachmat::create_numeric_output(nrows, ncols, oparam);
      
      Rcpp::IntegerVector vint;
      
      for (size_t ncol=0; ncol < ncols; ++ncol) {
        vint =  XX.column(ncol) ;
        out_dmat -> set_col(ncol, vint.begin());
      }
      
      return(out_dmat->yield());
    }
    else
    {
      return(wrap(XX));
    }
  }
  else if (dmtype==REALSXP) 
  {
    if ( X.isS4() == true){
      eX = read_DelayedArray_real( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }

    
    Rcpp::NumericMatrix XX = rcpp_xwxt( Rcpp::wrap(eX), w);;
    
    
    if ( X.isS4() == true)
    {
      ncols = XX.cols();
      nrows = XX.rows();
  
      // beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);  
      beachmat::output_param oparam(X);
      auto out_dmat = beachmat::create_numeric_output(XX.rows(), XX.cols(), oparam);
      
      Rcpp::NumericVector vint;
      
      if ( XX.cols()<= XX.rows())
      {
        // ESCRIVIM PER COLUMNES !!!
        for (size_t ncol=0; ncol<ncols; ++ncol) {
          vint = XX.column(ncol);
          out_dmat ->set_col(ncol, vint.begin());
        }
      }
      else
      {
        // ESCRIVIM PER FILES !!!
        for (size_t nrow=0; nrow<nrows; ++nrow) {
          vint = XX.row(nrow);
          out_dmat ->set_row(nrow, vint.begin());
        }  
      }
      return(out_dmat->yield());
    }
    else
    {
      return(wrap(XX));
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}

// [[Rcpp::export]]
Rcpp::RObject parxtwx(Rcpp::RObject X, Rcpp::RObject W)
{
  auto dmtype = beachmat::find_sexp_type(X);
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  Rcpp::List ret;
  Rcpp::NumericVector w = Rcpp::as<Rcpp::NumericVector >(W);
  
  // Rcpp::Rcout<<"Soc una matriu, i sóc S4? "<<X.isS4()<<"\n";
  // Rcpp::Rcout<<"Soc un vector, i sóc S4? "<<W.isS4()<<"\n";
  
  
  if ( dmtype == INTSXP ) {
    if ( X.isS4() == true){
      eX = read_DelayedArray_int( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }
    
    Rcpp::NumericMatrix XX = rcpp_xtwx( Rcpp::wrap(eX), w);
    if ( X.isS4() == true)
    {
    
      ncols = XX.cols();
      nrows = XX.rows();
      
      // beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);
      beachmat::output_param oparam(X);
      auto out_dmat = beachmat::create_numeric_output(nrows, ncols, oparam);
      
      Rcpp::IntegerVector vint;
      
      for (size_t ncol=0; ncol < ncols; ++ncol) {
        vint =  XX.column(ncol) ;
        out_dmat -> set_col(ncol, vint.begin());
      }
      
      return(out_dmat->yield());
    }
    else
    {
      return(wrap(XX));
    }
  }
  else if (dmtype==REALSXP) 
  {
    
    if ( X.isS4() == true){
      eX = read_DelayedArray_real( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }
    
    Rcpp::NumericMatrix XX = rcpp_xtwx( Rcpp::wrap(eX), w);;
    
    
    if ( X.isS4() == true)
    {
      ncols = XX.cols();
      nrows = XX.rows();
      
      // beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);  
      beachmat::output_param oparam(X);
      auto out_dmat = beachmat::create_numeric_output(XX.rows(), XX.cols(), oparam);
      
      Rcpp::NumericVector vint;
      
      if ( XX.cols()<= XX.rows())
      {
        // ESCRIVIM PER COLUMNES !!!
        for (size_t ncol=0; ncol<ncols; ++ncol) {
          vint = XX.column(ncol);
          out_dmat ->set_col(ncol, vint.begin());
        }
      }
      else
      {
        // ESCRIVIM PER FILES !!!
        for (size_t nrow=0; nrow<nrows; ++nrow) {
          vint = XX.row(nrow);
          out_dmat ->set_row(nrow, vint.begin());
        }  
      }
      return(out_dmat->yield());
      
    }
    else 
    {
      return(wrap(XX));
    }
    
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}



// [[Rcpp::export]]
Rcpp::RObject parXy(Rcpp::RObject X, Rcpp::RObject Y)
{
  auto dmtype = beachmat::find_sexp_type(X);
  size_t ncols = 0, nrows=0;
  Rcpp::NumericMatrix eX;
  Rcpp::NumericVector y = Rcpp::as<Rcpp::NumericVector >(Y);
  
  if ( dmtype == INTSXP ) {
    if ( X.isS4() == true){
      eX = read_DelayedArray_int_r( X );
    }else {
      eX = Rcpp::as<Rcpp::NumericMatrix >(X);
    }
    
    Rcpp::NumericMatrix XX = rcpp_parallel_Xy( eX, y);
    
    if ( X.isS4() == true)  {
      ncols = XX.cols();
      nrows = XX.rows();
      
      // beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);
      beachmat::output_param oparam(X);
      auto out_dmat = beachmat::create_numeric_output(nrows, ncols, oparam);
      
      Rcpp::IntegerVector vint;
      
      for (size_t ncol=0; ncol < ncols; ++ncol) {
        vint =  XX.column(ncol) ;
        out_dmat -> set_col(ncol, vint.begin());
      }
      
      return(out_dmat->yield());
    } else {
      return(XX);
    }
  } else if (dmtype==REALSXP) {
    if ( X.isS4() == true){
      eX = read_DelayedArray_real_r( X );
    }else {
      eX = Rcpp::NumericMatrix(X);
    }

    Rcpp::NumericMatrix XX = rcpp_parallel_Xy( eX, y);

    if ( X.isS4() == true)
    {
      ncols = XX.cols();
      nrows = XX.rows();
      
      // beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);  
      beachmat::output_param oparam(X);
      auto out_dmat = beachmat::create_numeric_output(XX.rows(), XX.cols(), oparam);
      
      Rcpp::NumericVector vint;
      
      if ( XX.cols()<= XX.rows())
      {
        // ESCRIVIM PER COLUMNES !!!
        for (size_t ncol=0; ncol<ncols; ++ncol) {
          vint = XX.column(ncol);
          out_dmat ->set_col(ncol, vint.begin());
        }
      } else{
        // ESCRIVIM PER FILES !!!
        for (size_t nrow=0; nrow<nrows; ++nrow) {
          vint = XX.row(nrow);
          out_dmat ->set_row(nrow, vint.begin());
        }  
      }
      return(out_dmat->yield());
    }  else  {
      return(XX);
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}








/*** R

*/


