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

struct MatrixVectMult: public RcppParallel::Worker {
  
  // input matrix to read from
  RcppParallel::RMatrix<double> mat;
  
  // input vector to read from
  RcppParallel::RVector<double> vect;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // other variables
  std::size_t numcol;
  
  
  // Constructor
  MatrixVectMult( Rcpp::NumericMatrix mat,  Rcpp::NumericVector vect, Rcpp::NumericMatrix rmat, std::size_t numcol)
    : mat(mat), vect(vect), rmat(rmat), numcol(numcol)  {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    
    for (std::size_t i = begin; i < end; i++) 
    {
      RcppParallel::RMatrix<double>::Row rowmat = mat.row(i); // rows we will operate on
      
      for (std::size_t j = 0; j < numcol ; j++) 
      {
          rmat(i,0) =  rmat(i,0) + (rowmat[j] * vect[j]);
          //if(i!=j)
          //  rmat(j,i) =  rmat(i,j);
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
        MatrixVectMult matvectmult( mat, y, rmat, y.size());
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
      
      beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);
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
  
      beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);  
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
      
      beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);
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
      
      beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);  
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
  Eigen::MatrixXd eX;
  Rcpp::List ret;
  Rcpp::NumericVector y = Rcpp::as<Rcpp::NumericVector >(Y);
  
  if ( dmtype == INTSXP ) {
    if ( X.isS4() == true){
      eX = read_DelayedArray_int( X );
    }else {
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }
    
    Rcpp::NumericMatrix XX = rcpp_parallel_Xy( Rcpp::wrap(eX), y);
    
    if ( X.isS4() == true)
    {
      
      ncols = XX.cols();
      nrows = XX.rows();
      
      beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);
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
    
    
    Rcpp::NumericMatrix XX = rcpp_parallel_Xy( Rcpp::wrap(eX), y);
    
    
    if ( X.isS4() == true)
    {
      ncols = XX.cols();
      nrows = XX.rows();
      
      beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);  
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



/*
 * CRIDES DIRECTES DES DE R, CODI REFET PER A GENERALITZAR CRIDES
 * INDEPENDENTMENT DEL TIPUS DE DADES UTILITZAT DES DE R, SI DES DE R S'UTILITZA:
 *    . DELAYED ARRAY: Fer conversió + executar codi multiplicació
 *    . HDF5  : Fer conversió + executar codi mult
 *    . Altres tipus de dades :   executar codi mult
 */
/*
// [[Rcpp::export]]
Rcpp::NumericMatrix xwxt(Rcpp::NumericMatrix mat, Rcpp::NumericVector w) {

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


// [[Rcpp::export]]
Rcpp::NumericMatrix xtwx(Rcpp::NumericMatrix mat, Rcpp::NumericVector w) {

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
*/





/*** R

library(RcppParallel)
library(rbenchmark) 
set.seed(5)
X <- matrix(sample(1:10,5000, replace = TRUE), ncol = 50);X
D <- DelayedArray(X);D
w <- sample(1:10,50, replace = TRUE);w
pmX <- parxwxt(X,w);pmX
pmX <- parxwxt(D,w);pmX

X%*%diag(w)%*%t(X)



RcppParallel::setThreadOptions(numThreads = 40)
n <- 1000
m <- matrix(runif(n*900), ncol = n)
w <- runif(n)



res <- benchmark(X%*%diag(w)%*%t(X),
                 parxwxt(D,w),
                 parxwxt(X,w),
                 order="relative", replications = c(10))
res[,1:4] 


res <- benchmark(t(X)%*%diag(w)%*%X,
                 parxtwx(D,w),
                 parxtwx(X,w),
                 order="relative", replications = c(10))
res[,1:4] 
dim(X)



n <- 5000
p <- 1000
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]

stopifnot(identical(t(M)%*%Y,parXy(t(M),Y)))
res <- benchmark(t(M)%*%Y,
                 parXy(t(M),Y),
                 order="relative", replications = c(10))
res[,1:4] 

parXy(t(M),Y)
dim(t(M))


*/


