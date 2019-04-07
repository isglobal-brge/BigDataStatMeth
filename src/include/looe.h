#ifndef looe
#define looe

  #include <RcppEigen.h>
  #include <RcppParallel.h>
  #include "parallel_MatrixVectorMultiplication.h"
  #include "parallel_MatrixMatrixMultiplication.h"
  #include "parallel_CrossProd.h"
  #include "vector_utilities.h"
  #include "ReadDelayedData.h"
  #include "svdDecomposition.h"
  #include "utilities.h"
  #include "matrix_utilities.h"
  #include "parallelBlockMult.h"
  // #include "optimizedproduct.h"


  // define a custom template binary functor
  template<typename cte> struct cteGinvDiagonal {
    EIGEN_EMPTY_STRUCT_CTOR(cteGinvDiagonal)
    double operator()(const cte& a, const cte& b) const { return double( std::pow((a/b),2) ); }
  };

  template<typename lambd> struct islambdazero {
    EIGEN_EMPTY_STRUCT_CTOR(islambdazero)
    double operator()(const lambd& a) const{
      if ((a==0) || Rcpp::NumericVector::is_na(a)==true)
        return 0;
      else
        return a;
    }
  };

  template<typename eval> struct poweval {
    EIGEN_EMPTY_STRUCT_CTOR(poweval)
    double operator()(const eval& a) const{
      return std::pow(a,2);
    }
  };
  

  
  double rcpplooei( double lambda, const Eigen::VectorXd& Lambda, 
                              const Eigen::MatrixXd& Q, const Eigen::VectorXd& Y, 
                              bool paral);

  Eigen::MatrixXd rcppinversecpp( double lambda, const Eigen::VectorXd& Lambda, 
                                             const Eigen::MatrixXd& Q, bool paral );
  
  Rcpp::RObject LOOE(Rcpp::RObject& X, Rcpp::RObject& Y, bool paral, 
                       Rcpp::Nullable<double> nl, Rcpp::Nullable<double> ml,
                       Rcpp::Nullable<Rcpp::RObject> l);
  


#endif