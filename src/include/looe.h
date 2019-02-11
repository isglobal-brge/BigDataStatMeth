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

  // using Rcpp::List;

  struct structlooe {
    Eigen::VectorXd vlooe;
    Eigen::MatrixXd Ginv;
    Eigen::VectorXd lambdas;
    double lambdamin;
  };
  
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


  Eigen::MatrixXd rcppinversecpp (Eigen::MatrixXd X, int lambda, bool eigen, Eigen::VectorXd Lambda, Eigen::MatrixXd Q );
  
  double rcpplooei( double lambda, Eigen::VectorXd Lambda, Eigen::MatrixXd Q, Eigen::VectorXd Y);
    
  Rcpp::RObject LOOE(Rcpp::RObject X, Rcpp::RObject Y, Rcpp::Nullable<double> nl,
                     Rcpp::Nullable<double> ml,
                     Rcpp::Nullable<Rcpp::RObject> l);
  


#endif