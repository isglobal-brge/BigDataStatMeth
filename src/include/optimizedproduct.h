#ifndef optimizedproduct
#define optimizedproduct

  #include <RcppEigen.h>
  #include "ReadDelayedData.h"

  Eigen::MatrixXd bdcrossproduct (Eigen::MatrixXd& mat);
  Eigen::MatrixXd bdtcrossproduct (Eigen::MatrixXd& mat);
  Eigen::MatrixXd xwxt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w);
  Eigen::MatrixXd xtwx(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w);

#endif

