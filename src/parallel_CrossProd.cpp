/*#include <RcppEigen.h>
#include <cmath>
#include <algorithm>
#include <RcppParallel.h>
*/
#include "include/parallel_CrossProd.h"
#include "include/ReadDelayedData.h"

// [[Rcpp::depends(RcppParallel)]]

// [[Rcpp::depends(RcppEigen)]]


/** tCROSSPROD **/
struct tCrosProdd : public RcppParallel::Worker  {
  
  // input matrix to read from
  const RcppParallel::RMatrix<double> mat;
  const RcppParallel::RMatrix<double> tmat;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // other variables
  std::size_t numcol;
  
  // Constructor
  tCrosProdd(const Rcpp::NumericMatrix mat, Rcpp::NumericMatrix tmat, Rcpp::NumericMatrix rmat, std::size_t numcol)
    : mat(mat), tmat(tmat), rmat(rmat), numcol(numcol) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    // size_t ncmat = mat.ncol();
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
        }
      }
    }
  }
};




struct tCrosProddEig : public RcppParallel::Worker  {
  
  // input matrix to read from
  const Eigen::Map<Eigen::MatrixXd> mat;
  const InversMat tmat;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // other variables
  std::size_t numcol;
  
  // Constructor
  tCrosProddEig(const Eigen::Map<Eigen::MatrixXd> mat, InversMat tmat, Rcpp::NumericMatrix rmat, std::size_t numcol)
    : mat(mat), tmat(tmat), rmat(rmat), numcol(numcol) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    // size_t ncmat = mat.ncol();
    for (std::size_t i = begin; i < end; i++) 
    {
      // RcppParallel::RMatrix<double>::Row rowmat = mat.row(i); // rows we will operate on
      Eigen::VectorXd rowmat = mat.row(i);
      
      for (std::size_t j = 0; j < numcol ; j++) 
      {
        //..// RcppParallel::RMatrix<double>::Column coltmat = tmat.column(j); // cols we will operate on 
        ColumnVector coltmat = tmat.col(j).transpose();
        size_t selem = tmat.col(i).size() ;
        for(std::size_t k = 0; k<selem; k++) 
        {
          rmat(i,j) =  rmat(i,j) + (rowmat[k] * coltmat[k]);
        }
      }
    }
  }
};


struct tCrosProdd_block : public RcppParallel::Worker  {
  
  // input matrix to read from
  const RcppParallel::RVector<double> mat;
  const RcppParallel::RVector<double> tmat;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // Other Variables
  double ncola;
  
  // Constructor
  tCrosProdd_block(const Rcpp::NumericVector mat, const Rcpp::NumericVector tmat, Rcpp::NumericMatrix rmat, double ncola)
    : mat(mat), tmat(tmat), rmat(rmat), ncola(ncola) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    std::size_t sj= rmat.ncol();
    
    for (std::size_t i = begin; i < end; i++) 
    {
      for(std::size_t j=0; j<sj; j++)
      {
        for(std::size_t k=0; k<ncola ; k++)
        {          
          rmat(i,j)=rmat(i,j) + mat[begin*ncola+k]*tmat[j*ncola+k];
        }
      }
    }
  }
};


Rcpp::NumericMatrix rcpp_parallel_tCrossProd(const Rcpp::NumericMatrix& mat) {
  
  // Transpose de matrix
  Rcpp::NumericMatrix tmat(Rcpp::transpose(mat));
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(mat.nrow(), tmat.ncol());
  
  
  try
  {
    if(mat.ncol()==tmat.nrow())
    {
      try{
        
        tCrosProdd tcrossprod(mat, tmat, rmat, tmat.ncol());
        RcppParallel::parallelFor(0, mat.nrow(), tcrossprod );

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



Eigen::Map<Eigen::MatrixXd> rcpp_parallel_tCrossProd_eigen(const Eigen::Map<Eigen::MatrixXd>& mat) {
  

  InversMat tmat(mat.cols(), mat.rows());
  tmat = mat.inverse();
  Rcpp::Rcout<<"\n"<<mat<<"\n";
  Rcpp::Rcout<<tmat<<"\n";
  
  
  // Eigen::Map<Eigen::Matrix> tmat(double, mat.cols(), mat.rows()) = mat.transpose();
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(mat.rows(), mat.rows());
  
  
  try
  {
    
    try{
      
      tCrosProddEig tcrossprod(mat, tmat, rmat, mat.rows());
      RcppParallel::parallelFor(0, mat.rows(), tcrossprod );

    }catch(std::exception &ex) {	
      forward_exception_to_r(ex);
    }

  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }  
  
  return (Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(rmat));
  
}



Rcpp::NumericMatrix rcpp_parallel_tCrossProd_block(Rcpp::NumericMatrix mat) {
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(mat.nrow(), mat.nrow());
  
  // Transpose de matrix
  Rcpp::NumericMatrix tmat(Rcpp::transpose(mat));
  
  // create the worker
  tCrosProdd_block tcrossprod(flatmatrm(mat), flatmatcm(tmat), rmat, double(mat.ncol()) );
  
  size_t sloop = size_t(mat.nrow()/500);
  // call it with parallelFor
  if(sloop>0)
  {
    for( size_t i = 0; i<sloop-1; i++)  {
      RcppParallel::parallelFor(i*500, (i+1)*500, tcrossprod);
    }
  } else {
    sloop = 1;
  }
  
  // Tractament final
  
  RcppParallel::parallelFor((sloop-1)*500, mat.nrow(), tcrossprod);
  
  
  return rmat;
}



/** CROSSPROD **/
struct CrosProdd : public RcppParallel::Worker  {
  
  // input matrix to read from
  const RcppParallel::RMatrix<double> mat;
  const RcppParallel::RMatrix<double> tmat;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // Constructor
  CrosProdd(const Rcpp::NumericMatrix mat, Rcpp::NumericMatrix tmat, Rcpp::NumericMatrix rmat)
    : mat(mat), tmat(tmat), rmat(rmat) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    // size_t ncmat = mat.ncol();
    for (std::size_t i = begin; i < end; i++) 
    {
      RcppParallel::RMatrix<double>::Row rowtmat = tmat.row(i); // rows we will operate on
      for (std::size_t j = 0; j < i+1 ; j++) 
      {
        RcppParallel::RMatrix<double>::Column colmat = mat.column(j); // cols we will operate on 
        size_t selem = colmat.length();
        for(std::size_t k = 0; k<selem; k++) 
        {
          rmat(i,j) =  rmat(i,j) + (rowtmat[k] * colmat[k]);
          if(i!=j)
            rmat(j,i) =  rmat(i,j);
          
        }
      }
    }
  }
};

struct CrosProdd_block : public RcppParallel::Worker  {
  
  // input matrix to read from
  const RcppParallel::RVector<double> mat;
  const RcppParallel::RVector<double> tmat;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // Other Variables
  int ncola;
  
  // Constructor
  /** WARNIN - R CMD Check () -- NOT TESTED  -- (order tmat-mat changed)
  CrosProdd_block(const Rcpp::NumericVector tmat, const Rcpp::NumericVector mat, Rcpp::NumericMatrix rmat, int ncola)
    : tmat(tmat), mat(mat), rmat(rmat), ncola(ncola) {}
  ***/
  CrosProdd_block(const Rcpp::NumericVector mat, const Rcpp::NumericVector tmat, Rcpp::NumericMatrix rmat, int ncola)
    : mat(mat), tmat(tmat), rmat(rmat), ncola(ncola) {}
  
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    std::size_t sj= rmat.ncol();
    //.commented 20201120 - warning check().//double kp = 0;
    
    for (std::size_t i = begin; i < end; i++) 
    {
      for(std::size_t j=0; j<sj; j++)
      {
        for(std::size_t k=0; k<ncola ; k++)
        { 
          rmat(i,j)=rmat(i,j) + mat[begin*ncola+k]*tmat[j*ncola+k];
        }
      }
    }
  }
};

Rcpp::NumericMatrix rcpp_parallel_CrossProd(Rcpp::NumericMatrix mat) {
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(mat.ncol(), mat.ncol());
  
  // Transpose de matrix
  Rcpp::NumericMatrix tmat(Rcpp::transpose(mat));
  
  
  try
  {
    if(mat.ncol()==tmat.nrow())
    {
      try{
        
        // create the worker
        CrosProdd crossprod(mat, tmat, rmat);
        RcppParallel::parallelFor(0, tmat.nrow(), crossprod);
        
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

Rcpp::NumericMatrix rcpp_parallel_CrossProd_block(Rcpp::NumericMatrix mat) {
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(mat.ncol(), mat.ncol());
  
  // Transpose de matrix
  Rcpp::NumericMatrix tmat(Rcpp::transpose(mat));
  
  // create the worker
  CrosProdd_block crossprod(flatmatrm(tmat), flatmatcm(mat), rmat, tmat.ncol() );
  
  // call it with parallelFor
  RcppParallel::parallelFor(0, tmat.nrow(), crossprod);
  //RcppParallel::parallelReduce(0, mat.nrow(), crossprod);
  
  return rmat;
}


/** CRIDES R **/

// [[Rcpp::export]]
Rcpp::RObject partCrossProd(Rcpp::RObject X)
{
  auto dmtype = beachmat::find_sexp_type(X);
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP ) {
    
    if ( X.isS4() == true){
      eX = read_DelayedArray_int( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }

    Rcpp::NumericMatrix XX = rcpp_parallel_tCrossProd(Rcpp::wrap(eX));
    
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
    }else {
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
    
    Rcpp::NumericMatrix XX = rcpp_parallel_tCrossProd(Rcpp::wrap(eX));
    
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
    }else {
      return(wrap(XX));
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}

// [[Rcpp::export]]
Rcpp::RObject partCrossProd_block(Rcpp::RObject X)
{
  auto dmtype = beachmat::find_sexp_type(X);
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP ) {
    
    if ( X.isS4() == true){
      eX = read_DelayedArray_int( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }
    
    Rcpp::NumericMatrix XX = rcpp_parallel_tCrossProd_block(Rcpp::wrap(eX));
    
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
    } else {
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
    
    Rcpp::NumericMatrix XX = rcpp_parallel_tCrossProd_block(Rcpp::wrap(eX));
    
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
    } else {
      return(wrap(XX));
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
}

// [[Rcpp::export]]
Rcpp::RObject parCrossProd(Rcpp::RObject X)
{
  auto dmtype = beachmat::find_sexp_type(X);
  
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP ) {
    
    if ( X.isS4() == true){
      eX = read_DelayedArray_int( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }
    
    Rcpp::NumericMatrix XX = rcpp_parallel_CrossProd(Rcpp::wrap(eX));
    
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
    }else {
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
    
    Rcpp::NumericMatrix XX = rcpp_parallel_CrossProd(Rcpp::wrap(eX));
    
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
    }else {
      return(wrap(XX));
    }
  } else {
    Rcpp::Rcout<<"Error aquí";
    throw std::runtime_error("unacceptable matrix type");
  }
  
}

// [[Rcpp::export]]
Rcpp::RObject parCrossProd_block(Rcpp::RObject X)
{
  auto dmtype = beachmat::find_sexp_type(X);
  
  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd eX;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP ) {
    
    if ( X.isS4() == true){
      eX = read_DelayedArray_int( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }
    
    Rcpp::NumericMatrix XX = rcpp_parallel_CrossProd_block(Rcpp::wrap(eX));
    
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
    } else {
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
    
    Rcpp::NumericMatrix XX = rcpp_parallel_CrossProd_block(Rcpp::wrap(eX));
    
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
    }else {
      return(wrap(XX));
    }
    
  } else {
    Rcpp::Rcout<<"Error aquí";
    throw std::runtime_error("unacceptable matrix type");
  }
  
}



// [[Rcpp::export]]
Rcpp::RObject partCrossProdEigen(Rcpp::RObject X)
{
  
  if( X.sexp_type() == INTSXP)
  {
    const Eigen::Map<Eigen::MatrixXi> eX(Rcpp::as<Eigen::Map<Eigen::MatrixXi> >(X)); 
    Rcpp::Rcout<<"\n"<<eX<<"\n";
  } else if(X.sexp_type() == REALSXP) {
    const Eigen::Map<Eigen::MatrixXd> eX(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X));
    Rcpp::Rcout<<"\n"<<eX<<"\n";
  }
  /*
  Eigen::MatrixXd eX = Rcpp::as<Eigen::MatrixXd>(X);
  const Eigen::Map<Eigen::MatrixXd> A(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X));
  
  Rcpp::Rcout<<"\n"<<A<<"\n";
  //Eigen::Map<Eigen::MatrixXd> XX = rcpp_parallel_tCrossProd_eigen(eX);
  Eigen::Map<Eigen::MatrixXd> XX = A;
  return(Rcpp::wrap(XX));*/
  return(Rcpp::wrap('2'));

}












/***R

*/