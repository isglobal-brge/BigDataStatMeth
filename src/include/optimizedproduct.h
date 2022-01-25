#ifndef optimizedproduct
#define optimizedproduct

    // #define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
    
    #include <RcppEigen.h>
    #include "BigDataStatMeth.h"
    //#include <omp.h>
    #include "pkg_omp.h"
    #include <thread>
    #include "ReadDelayedData.h"
    #include "parallelBlockMult_hdf5.h"
    
    Eigen::MatrixXd bdcrossproduct (Eigen::MatrixXd& mat);
    Eigen::MatrixXd bdtcrossproduct (Eigen::MatrixXd& mat);
    Eigen::MatrixXd xwxt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w);
    Eigen::MatrixXd xtwx(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w);
    Eigen::MatrixXd Xwd(const Eigen::MatrixXd& X, const Eigen::VectorXd& w);
    Eigen::MatrixXd Xwd_parallel(const Eigen::MatrixXd& X, const Eigen::VectorXd& w, Rcpp::Nullable<int> threads);
    Eigen::MatrixXd wdX_parallel(const Eigen::MatrixXd& X, const Eigen::VectorXd& w, Rcpp::Nullable<int> threads);
    //  Eigen::MatrixXd bdXwd(Rcpp::RObject X, Rcpp::RObject w, std::string op,
    //                        Rcpp::Nullable<bool> bparal, Rcpp::Nullable<int> threads);
    Eigen::MatrixXd bdScalarwproduct(Rcpp::RObject a, double w, std::string op);
    Eigen::MatrixXd bdwproduct(Rcpp::RObject X, Rcpp::RObject w, std::string op);

#endif

