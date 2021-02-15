#ifndef svdDecomposition
#define svdDecomposition

  #include <RcppEigen.h>
  #include "parallel_CrossProd.h"
  #include "matrix_utilities.h"
  #include "hdf5_to_Eigen.h"
  #include "svdutils.h"
  #include "parallelBlockMult_hdf5.h"
  #include "svdBlockDecomposition_hdf5.h"
  #include "pseudoinv.h"
  #include "spectra/SymEigsSolver.h"   // To access symmetric matrix
  #include "spectra/SymGEigsSolver.h"  // To access symmetric matrix
  #include "beachmat/numeric_matrix.h" // To access numeric matrix
  #include "beachmat/integer_matrix.h" // To access numeric matrix
  #include "optimizedproduct.h"
  #include "tgmath.h"


  #define MAXSVDBLOCK 5120

  // dgesvd_ is a symbol in the LAPACK-BLAS Level 3 
  //    DGESVD computes the singular value decomposition (SVD) of a real M-by-N matrix A, 
  //       optionally computing the left and/or right singular vectors
  extern "C" {
    extern void dgesvd_( char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
  }

  svdeig RcppbdSVD( Eigen::MatrixXd& X, int k = 0, int nev = 0, bool bcenter = true, bool bscale = true);
  svdeig RcppbdSVD_lapack( Eigen::MatrixXd& X,  bool bcenter, bool bscale );
  svdeig RcppbdSVD_hdf5( std::string filename, std::string strsubgroup, std::string strdataset,  
                         int k, int q, int nev, bool bcenter, bool bscale, Rcpp::Nullable<int> ithreads  );
  svdeig RcppbdSVD_hdf5_Block( H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
                               int irows, int icols, Rcpp::Nullable<int> threads);
  
  //..// Rcpp::RObject BDCInverse_Cholesky (const Rcpp::RObject & x );
  
  svdeig RcppCholDec(const Eigen::MatrixXd& X);
  
  Rcpp::RObject bdSVD (const Rcpp::RObject & x, int k=0, int nev=0, bool bcenter=true, bool bscale = true );
  //..// Rcpp::RObject bdSVD_hdf5(const Rcpp::RObject & x, CharacterVector group = R_NilValue, CharacterVector dataset = R_NilValue,
  //..//                                 int parts = 2, int k=0, int nev=0, bool bcenter=true, bool bscale=true, Rcpp::Nullable<int> threads = R_NilValue);
  Rcpp::RObject bdSVD_hdf5(const Rcpp::RObject & x, CharacterVector group = R_NilValue, CharacterVector dataset = R_NilValue,
                           int parts = 2, int k=0, bool bcenter=true, bool bscale=true, Rcpp::Nullable<int> threads = R_NilValue);
  
  Rcpp::RObject bdSVD_lapack(const Rcpp::RObject & x, Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true);

#endif