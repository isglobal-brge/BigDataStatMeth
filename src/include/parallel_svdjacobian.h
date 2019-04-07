#ifndef parallel_svdjacobian
#define parallel_svdjacobian

  #include <RcppEigen.h>
  #include <iostream>
  #include <cmath>
  #include <algorithm>
  #include <stdio.h>
  #include <stdlib.h>
  #include <fstream>
  #include <sys/time.h>
  #include <omp.h>

  #include "redsvd.h"
  #include "ReadDelayedData.h"

  #define epsilon 1.e-8
  #define num 16

  template <typename T> double sgn(T val)
  {
    return (val > T(0)) - (val < T(0));
  }

  void svdjacob (Eigen::MatrixXd U_t, int M, int N, Eigen::MatrixXd& U, Eigen::MatrixXd& V, Eigen::VectorXd& S, double &error, int &iter);
  
#endif
