#ifndef svdDecomposition
#define svdDecomposition

    #include <RcppEigen.h>
    #include "spectra/SymEigsSolver.h"   // To access symmetric matrix
    #include "spectra/SymGEigsSolver.h"  // To access symmetric matrix
    #include "tgmath.h"
    #include "matrix_utilities.h"
    #include "hdf5_to_Eigen.h"
    #include "svdutils.h"
    #include "parallelBlockMult_hdf5.h"
    #include "svdBlockDecomposition_hdf5.h"
    #include "pseudoinv.h"
    #include "beachmat/numeric_matrix.h" // To access numeric matrix
    #include "beachmat/integer_matrix.h" // To access numeric matrix
    #include "optimizedproduct.h"
  


  #define MAXSVDBLOCK 1500

  // dgesvd_ is a symbol in the LAPACK-BLAS Level 3 
  //    DGESVD computes the singular value decomposition (SVD) of a real M-by-N matrix A, 
  //       optionally computing the left and/or right singular vectors
  extern "C" {
    extern void dgesvd_( char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
  }
  
  extern "C" {
      extern void dgesdd_( char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*, int*);
  }

  
  // C++ functions
  svdeig RcppbdSVD( Eigen::MatrixXd& X, int k = 0, int nev = 0, bool bcenter = true, bool bscale = true);
  
  svdeig RcppbdSVD_lapack( Eigen::MatrixXd& X,  bool bcenter, bool bscale, bool complete );
  svdeig RcppbdSVD_lapack_optim( Eigen::MatrixXd& X,  bool bcenter, bool bscale, bool complete );
  
  svdeig RcppbdSVD_hdf5( std::string filename, std::string strsubgroup, std::string strdataset,  
                         int k, int q, int nev, bool bcenter, bool bscale, double dthreshold, 
                         Rcpp::Nullable<int> ithreads);
  
  svdeig RcppbdSVD_hdf5_ptr( H5File* file, std::string strsubgroup, std::string strdataset,  
                             int k, int q, int nev, bool bcenter, bool bscale, bool bstorehdf5, 
                             double dthreshold, Rcpp::Nullable<int> ithreads);
  
  svdeig RcppbdSVD_hdf5_Block( H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
                               int irows, int icols, double dthreshold, Rcpp::Nullable<int> threads );
  
  svdeig RcppCholDec(const Eigen::MatrixXd& X);
  
  
  // R functions (Calls)
  Rcpp::RObject bdSVD (const Rcpp::RObject & X, int k=0, int nev=0, bool bcenter=true, bool bscale = true );
  
  Rcpp::RObject bdSVD_hdf5 (const Rcpp::RObject & file, Rcpp::Nullable<CharacterVector> group, 
                            Rcpp::Nullable<CharacterVector> dataset,
                            Rcpp::Nullable<int> k, Rcpp::Nullable<int> q,
                            Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale,
                            Rcpp::Nullable<double> rankthreshold,
                            Rcpp::Nullable<int> threads);
  
  Rcpp::RObject bdSVD_lapack(const Rcpp::RObject & X, Rcpp::Nullable<bool> bcenter=true, 
                             Rcpp::Nullable<bool> bscale=true, Rcpp::Nullable<bool> complete = false);

#endif
