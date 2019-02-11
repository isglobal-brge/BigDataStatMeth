#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision


double prova( double a, double b)
{
  Rcpp::Rcout<<"a : "<<a<<"b : "<<b<<"\n";
  return(a+b);
}





/*
// define a custom template binary functor
template<typename cte> struct cteGinvDiagonal {
  EIGEN_EMPTY_STRUCT_CTOR(cteGinvDiagonal)
  double operator()(const cte& a, const cte& b) const { return double( std::pow((a/b),2) ); }
};

 
// [[Rcpp::export]]
void getbinaryExpression() {
  Eigen::VectorXd a (4) ;
  a<< 1,2,3,4;
  Eigen::VectorXd b (4) ;
  b<< 5,6,7,8;
  Eigen::VectorXd c (4) ;
  
  c = a.binaryExpr(b, cteGinvDiagonal<double>());
  
  Rcpp::Rcout<<"a / b : "<<c<<" Sumatori : "<<c.sum()<<"\n";
  
}
*/

// [[Rcpp::export]]
void Proves() {
  double lambda = 0.01;
  Eigen::VectorXd a (4) ;
  a<< 1,2,3,4;
  
  
  
  Rcpp::Rcout<<"a +lambda : "<<a.array()+lambda<<"\n";
  Rcpp::Rcout<<"a +lambda : "<<1/(a.array()+lambda)<<"\n";
}

/*** R
#getbinaryExpression()
Proves()


*/
