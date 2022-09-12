#ifndef hdf5_computeMatrixVector
#define hdf5_computeMatrixVector

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    
    // C++ functions
    
    // by Cols
    Eigen::MatrixXd Rcpp_matrixVectorMultiplication_byRow(Eigen::MatrixXd X, Eigen::VectorXd v);
    Eigen::MatrixXd Rcpp_matrixVectorSum_byRow(Eigen::MatrixXd X, Eigen::VectorXd v);
    Eigen::MatrixXd Rcpp_matrixVectorSubstract_byRow(Eigen::MatrixXd X, Eigen::VectorXd v);
    Eigen::MatrixXd Rcpp_matrixVectorDivision_byRow(Eigen::MatrixXd X, Eigen::VectorXd v);
    
    // by Columns
    Eigen::MatrixXd Rcpp_matrixVectorMultiplication_byCol(Eigen::MatrixXd X, Eigen::VectorXd v);
    Eigen::MatrixXd Rcpp_matrixVectorSum_byCol(Eigen::MatrixXd X, Eigen::VectorXd v);
    Eigen::MatrixXd Rcpp_matrixVectorSubstract_byCol(Eigen::MatrixXd X, Eigen::VectorXd v);
    Eigen::MatrixXd Rcpp_matrixVectorDivision_byCol(Eigen::MatrixXd X, Eigen::VectorXd v);
    
    
    
    // R functions
    void bdcomputeMatrixVector_hdf5( std::string filename, std::string group, 
                                     std::string dataset,
                                     std::string vectorgroup, std::string vectordataset,
                                     std::string outdataset, 
                                     std::string func,
                                     Rcpp::Nullable<std::string> outgroup,
                                     Rcpp::Nullable<bool> byrows,
                                     Rcpp::Nullable<int> force);
#endif
