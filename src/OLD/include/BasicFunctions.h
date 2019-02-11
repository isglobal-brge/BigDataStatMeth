#ifndef Basic_Functions
#define Basic_Functions

#include <RcppEigen.h>

Rcpp::RObject Normalize_Data (const Rcpp::RObject & x );
Eigen::MatrixXd CrossProduct_eig (const Eigen::MatrixXd x, bool tpc );
Eigen::MatrixXd CrossProduct (const Rcpp::RObject & x, bool tpc );

#endif