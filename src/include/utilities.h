#ifndef utilities
#define utilities

  #include <RcppEigen.h>
  #include <RcppParallel.h>

  bool double_equals(double a, double b);
  void setZeros(Eigen::Map<Eigen::VectorXd>* lambda);
  Eigen::VectorXd equalZero(Eigen::VectorXd v);

#endif