#ifndef svdDecomposition
#define svdDecomposition

  #include <RcppEigen.h>
  #include "parallel_CrossProd.h"
  #include "matrix_utilities.h"
  #include "spectra/SymEigsSolver.h"  // To access symmetric matrix
  #include "spectra/SymGEigsSolver.h"  // To access symmetric matrix
  #include "beachmat/numeric_matrix.h"  // To access numeric matrix
  #include "beachmat/integer_matrix.h"  // To access numeric matrix
  #include "optimizedproduct.h"

  struct svdeig {
    Eigen::VectorXd d;
    Eigen::MatrixXd u;
    Eigen::MatrixXd v;
  };
  
  struct svd {
    Rcpp::NumericVector d;
    Rcpp::NumericMatrix u;
    Rcpp::NumericMatrix v;
  };

  svdeig RcppBDsvd_eig(  Eigen::MatrixXd& X, int k = 0, int nev = 0, bool normalize = true );
  svd RcppBDsvd( const Rcpp::NumericMatrix& X, int k = 0, int nev = 0, bool normalize = true );
  Rcpp::RObject BDCInverse_Cholesky (const Rcpp::RObject & x );


#endif