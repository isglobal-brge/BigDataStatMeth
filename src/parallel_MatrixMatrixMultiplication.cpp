#include "include/parallel_MatrixMatrixMultiplication.h"
#include "include/ReadDelayedData.h"


/** MATRIX-MATRIX MULTIPLICATION **/

struct XYProd : public RcppParallel::Worker  {
  
  // input matrix to read from
  const RcppParallel::RMatrix<double> matX;
  const RcppParallel::RMatrix<double> matY;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // other variables
  std::size_t numcol;

  // Constructor
  XYProd(const Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY, Rcpp::NumericMatrix rmat, std::size_t numcol)
    : matX(matX), matY(matY), rmat(rmat), numcol(numcol) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    // size_t ncmat = mat.ncol();
    for (std::size_t i = begin; i < end; i++) 
    {
      RcppParallel::RMatrix<double>::Row rowmat = matX.row(i); // rows we will operate on
      for (std::size_t j = 0; j < numcol ; j++) 
      {
        RcppParallel::RMatrix<double>::Column coltmat = matY.column(j); // cols we will operate on 
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



struct XYProdEig : public RcppParallel::Worker  {
  
  // input matrix to read from
  //const RcppParallel::RMatrix<double> matX;
  //const RcppParallel::RMatrix<double> matY;
  Eigen::MatrixXd matX;
  Eigen::MatrixXd matY;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // other variables
  std::size_t numcol;
  
  // Constructor
  XYProdEig(const Eigen::MatrixXd matX, Eigen::MatrixXd matY, Rcpp::NumericMatrix rmat, std::size_t numcol)
    : matX(matX), matY(matY), rmat(rmat), numcol(numcol) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    // size_t ncmat = mat.ncol();
    for (std::size_t i = begin; i < end; i++) 
    {
      // RcppParallel::RMatrix<double>::Row rowmat = matX.row(i); // rows we will operate on
      Eigen::VectorXd rowmat = matX.row(i);
      for (std::size_t j = 0; j < numcol ; j++) 
      {
        // RcppParallel::RMatrix<double>::Column coltmat = matY.column(j); // cols we will operate on 
        ColumnVector coltmat = matY.col(j);
        size_t selem = matY.col(j).size();
        
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


struct XYProdBlock : public RcppParallel::Worker  {
  
  // input matrix to read from
  const RcppParallel::RVector<double> vecX;
  const RcppParallel::RVector<double> vecY;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // other variables
  std::size_t ncY;
  std::size_t ncX;
  
  // Constructor
  XYProdBlock(const Rcpp::NumericVector vecX, Rcpp::NumericVector vecY, Rcpp::NumericMatrix rmat, std::size_t ncY, std::size_t ncX)
    : vecX(vecX), vecY(vecY), rmat(rmat), ncY(ncY), ncX(ncX) {}
  
  void operator()(std::size_t begin, std::size_t end) 
  {
    for (std::size_t i = begin; i < end; i++) 
    {
      for(std::size_t j = 0; j<ncY; j++) 
      {
        for( std::size_t k=0; k<ncX; k++)
        {
          rmat(i,j) =  rmat(i,j) + (vecX[i*ncX + k] * vecY[k*ncY + j]);
        }
      }
    }
  }
};



Rcpp::NumericMatrix rcpp_parallel_XYProd(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY) {
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(matX.nrow(), matY.ncol());

  try
  {
    if(matX.ncol()==matY.nrow())
    {
      try{
        
        // create the worker
        XYProd xyprod(matX, matY, rmat, matY.ncol());
        RcppParallel::parallelFor(0, matX.nrow(), xyprod);
        
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


Eigen::MatrixXd rcpp_parallel_XYProd_eigen(Eigen::MatrixXd matX, Eigen::MatrixXd matY) {
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(matX.rows(), matY.cols());
  
  try
  {
    if(matX.cols()==matY.rows())
    {
      try{
        
        // create the worker
        XYProdEig xyprod(matX, matY, rmat, matY.cols());
        RcppParallel::parallelFor(0, matX.rows(), xyprod);
        
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
  
  return (Rcpp::as<Eigen::MatrixXd>(rmat));
  
}




Rcpp::NumericMatrix rcpp_parallel_XYProdBlock(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY) {
  
  // other varialbes
  double_t ncX = matX.ncol();
  double_t nrX = matX.nrow();
  double_t ncY = matY.ncol();
  double_t nrY = matY.nrow();
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(nrX, ncY);

  try
  {
    if( ncX == nrY)
    {
      try{
        
        // create the worker
        XYProdBlock xyprodblock(flatmatrm(matX), flatmatrm(matY), rmat, ncY, ncX );
        RcppParallel::parallelFor(0, nrX, xyprodblock);
        
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



Rcpp::NumericMatrix rcpp_parallel_XtYProd(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY) {
  
  // Transpose matrix X
  Rcpp::NumericMatrix tmat(Rcpp::transpose(matX));
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(tmat.nrow(), matY.ncol());
  
  
  try
  {
    if(tmat.ncol()==matY.nrow())
    {
      try{
        
        // create the worker
        XYProd xyprod(tmat, matY, rmat, matY.ncol());
        RcppParallel::parallelFor(0, tmat.nrow(), xyprod);
        
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

Rcpp::NumericMatrix rcpp_parallel_XYtProd(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY) {
  
  // Transpose matrix X
  Rcpp::NumericMatrix tmat(Rcpp::transpose(matY));
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(matX.nrow(), tmat.ncol());
  
  try
  {
    if(matX.ncol()==tmat.nrow())
    {
      try{
        
        // create the worker
        XYProd xyprod(matX, tmat, rmat, tmat.ncol());
        RcppParallel::parallelFor(0, matX.nrow(), xyprod);
        
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

//' Matrix multiplication with Delayed Array Object (RcppParallel)
//' 
//' This function performs a block matrix-matrix multiplication with numeric matrix or Delayed Arrays
//' 
//' @param a a double matrix.
//' @param b a double matrix.
//' @param op, (optional, default = "xy"), if op="xy" then performs the x\%*\%y matrix multiplication, if op = "xty" preforms t(X)\%*\% Y, if op = "xyt" performs X\%*\%t(Y)
//' @return numerical matrix
//' @examples
//' # with numeric matrix
//' m <- 500
//' k <- 300
//' n <- 400
//' A <- matrix(rnorm(n*p), nrow=n, ncol=k)
//' B <- matrix(rnorm(n*p), nrow=k, ncol=n)
//' 
//' parXYProd(A,B,128, TRUE)
//' 
//' # with Delaeyd Array
//' AD <- DelayedArray(A)
//' BD <- DelayedArray(B)
//' 
//' parXYProd(AD,BD,128, TRUE)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject parXYProd(Rcpp::RObject X, Rcpp::RObject Y, Rcpp::Nullable<std::string> op = R_NilValue)
{

  size_t ncols = 0, nrows=0;
  
  Eigen::MatrixXd eX;
  Eigen::MatrixXd eY;
  Rcpp::List ret;
  Rcpp::NumericMatrix XY;


  if ( X.isS4() == true){
    eX = read_DelayedArray(X);
  }else {
    try{
      eX = Rcpp::as<Eigen::MatrixXd >(X);
    }catch(std::exception &ex) {
      eX = Rcpp::as<Eigen::VectorXd >(X);
    }
  }
  
  
  if ( Y.isS4() == true){
    eY = read_DelayedArray(Y);
  }else {
    try{
      eY = Rcpp::as<Eigen::MatrixXd >(Y);  
    }catch(std::exception &ex) {
      eY = Rcpp::as<Eigen::VectorXd >(Y);
    }
    
  }
  
  if (op.isNotNull())
  {
    std::string oper = Rcpp::as<std::string>(op); 
    if( oper == "xy" )
      XY =  rcpp_parallel_XYProd(Rcpp::wrap(eX), Rcpp::wrap(eY));
    else if(oper == "xty" )
      XY =  rcpp_parallel_XtYProd(Rcpp::wrap(eX), Rcpp::wrap(eY));
    else if(oper == "xyt" )
      XY =  rcpp_parallel_XYtProd(Rcpp::wrap(eX), Rcpp::wrap(eY));
    else
      throw std::runtime_error("operation not recognized");
  }else {
    XY =  rcpp_parallel_XYProd(Rcpp::wrap(eX), Rcpp::wrap(eY));  
  }
  
  
  if ( X.isS4() == true || Y.isS4() == true)
  {
    ncols = XY.cols();
    nrows = XY.rows();
    
    // beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);
    beachmat::output_param oparam(X);
    auto out_dmat = beachmat::create_numeric_output(nrows, ncols, oparam);
    
    Rcpp::IntegerVector vint;
    
    for (size_t ncol=0; ncol < ncols; ++ncol) {
      vint =  XY.column(ncol) ;
      out_dmat -> set_col(ncol, vint.begin());
    }
    
    return(out_dmat->yield());
    
  }else {
    
    return(wrap(XY));
  }
}



  
//' Block matrix multiplication with Delayed Array Object
//' 
//' This function performs a block matrix-matrix multiplication with numeric matrix or Delayed Arrays
//' 
//' @param a a double matrix.
//' @param b a double matrix.
//' @param op, (optional, default = "xy"), if op="xy" then performs the x\%*\%y matrix multiplication, if op = "xty" preforms t(X)\%*\% Y, if op = "xyt" performs X\%*\%t(Y)
//' @return numerical matrix
//' @examples
//' # with numeric matrix
//' m <- 500
//' k <- 1500
//' n <- 400
//' A <- matrix(rnorm(n*p), nrow=n, ncol=k)
//' B <- matrix(rnorm(n*p), nrow=k, ncol=n)
//' 
//' Bblockmult(A,B,128, TRUE)
//' 
//' # with Delaeyd Array
//' AD <- DelayedArray(A)
//' BD <- DelayedArray(B)
//' 
//' Bblockmult(AD,BD,128, TRUE)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject parXYProdBlock(Rcpp::RObject X, Rcpp::RObject Y, Rcpp::Nullable<std::string> op = R_NilValue)
{
    auto dmtypex = beachmat::find_sexp_type(X);
    auto dmtypey = beachmat::find_sexp_type(Y);
    
    size_t ncols = 0, nrows=0;
    
    Eigen::MatrixXd eX;
    Eigen::MatrixXd eY;
    Rcpp::List ret;
    Rcpp::NumericMatrix XY;
    
    
    if ( X.isS4() == true){
      eX = read_DelayedArray(X);
    }else {
      try{
        eX = Rcpp::as<Eigen::MatrixXd >(X);
      }catch(std::exception &ex) {
        eX = Rcpp::as<Eigen::VectorXd >(X);
      }
    }
    
    
    if ( Y.isS4() == true){
      eY = read_DelayedArray(Y);
    }else {
      try{
        eY = Rcpp::as<Eigen::MatrixXd >(Y);  
      }catch(std::exception &ex) {
        eY = Rcpp::as<Eigen::VectorXd >(Y);
      }
      
    }
    
    if (op.isNotNull())
    {
      std::string oper = Rcpp::as<std::string>(op); 
      if( oper == "xy" )
        XY =  rcpp_parallel_XYProdBlock(Rcpp::wrap(eX), Rcpp::wrap(eY));
      else if(oper == "xty" )
        XY =  rcpp_parallel_XYProdBlock(Rcpp::wrap(eX), Rcpp::wrap(eY));
      else if(oper == "xyt" )
        XY =  rcpp_parallel_XYProdBlock(Rcpp::wrap(eX), Rcpp::wrap(eY));
      else
        throw std::runtime_error("operation not recognized");
    }else {
      XY =  rcpp_parallel_XYProdBlock(Rcpp::wrap(eX), Rcpp::wrap(eY));  
    }
    
    
    if ( X.isS4() == true || Y.isS4() == true)
    {
      ncols = XY.cols();
      nrows = XY.rows();
      
      // beachmat::output_param oparam(beachmat::DELAYED, FALSE, TRUE);
      beachmat::output_param oparam(X);
      auto out_dmat = beachmat::create_numeric_output(nrows, ncols, oparam);
      
      Rcpp::IntegerVector vint;
      
      for (size_t ncol=0; ncol < ncols; ++ncol) {
        vint =  XY.column(ncol) ;
        out_dmat -> set_col(ncol, vint.begin());
      }
      
      return(out_dmat->yield());
      
    }else {
      
      return(wrap(XY));
    }
  }
  
  
  




/*** R

library(rbenchmark)
library(DelayedArray)
library(RcppParallel)



set.seed(5)
X <- matrix(sample(1:10,200, replace = TRUE), ncol = 10);X
set.seed(5)
Y <- matrix(sample(1:10,70, replace = TRUE), nrow = 10);Y



parXYProdBlock(X,Y)
X%*%Y

stopifnot(identical(X%*%Y,parXYProd(X,Y)))
stopifnot(identical(parXYProd(X,Y), parXYProdBlock(X,Y)))

n <- 100
set.seed(4923)
X <- matrix(runif(n*9000), ncol = n)
set.seed(4237)
Y <- matrix(runif(n*9000), nrow = n)
DX <- DelayedArray(X)
DY <- DelayedArray(Y)

  
RcppParallel::setThreadOptions(numThreads = 40)
setThreadOptions(numThreads = defaultNumThreads() / 2)
res <- benchmark(
                 parXYProd(X,Y),
                 parXYProd(DX,DY),
                 parXYProdBlock(X,Y),
                 parXYProdBlock(DX,DY),
                 order="relative", replications = c(3))
res[,1:4] 


parXYProd(X,Y,"xty")
t(X)%*%Y
t(Y)%*%X
parXYProd(DX,DY,"xyt")
X%*%t(Y)


*/