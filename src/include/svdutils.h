#ifndef svdutils
#define svdutils


  #include <RcppEigen.h>


  struct svdeig {
    Eigen::VectorXd d;
    Eigen::MatrixXd u;
    Eigen::MatrixXd v;
    bool bokuv;
    bool bokd;
  };


  Rcpp::IntegerVector getInitialPosition( bool transp, int desp );
  Rcpp::IntegerVector getSizetoRead( bool transp, int count, int rows, int cols );

#endif

