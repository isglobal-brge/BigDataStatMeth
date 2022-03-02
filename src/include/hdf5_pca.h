#ifndef hdf5_pca
#define hdf5_pca

    #include <RcppEigen.h>
    #include "matrix_utilities.h"
    #include "svdDecomposition.h"
    #include "optimizedproduct.h"
    //   #include "hdf5_to_Eigen.h"
    #include "rhdf5Utils.h"
    #include "parallelBlockMult_hdf5.h"
    //#include <omp.h>
    #include "pkg_omp.h"
  

    struct pcaeig {
        Eigen::MatrixXd varCoord;
        Eigen::MatrixXd Y;
        Eigen::MatrixXd var;
        Eigen::MatrixXd percvar;
        Eigen::MatrixXd components;
        std::string hdf5file = "";
    };
    
    
    // C++ functions
    void get_HDF5_PCA_variance_ptr( H5File* file, std::string strdataset);
    void get_HDF5_PCA_variables_ptr(  H5File* file, std::string strdataset);
    
    
    // R functions (Calls)
    void bdPCA_hdf5(std::string filename, std::string group, std::string dataset,
                    Rcpp::Nullable<int> ncomponents,
                    Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale, 
                    Rcpp::Nullable<int> k, Rcpp::Nullable<int> q,
                    Rcpp::Nullable<double> rankthreshold,
                    Rcpp::Nullable<bool> force, Rcpp::Nullable<int> threads);

#endif
