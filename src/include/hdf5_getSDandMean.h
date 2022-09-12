#ifndef hdf5_getSDandMean
#define hdf5_getSDandMean

    #include <RcppEigen.h>
    #include "rhdf5Utils.h"
    
    
    // functions called from c++
    int get_HDF5_mean_sd_by_column_ptr(H5File* file, DataSet* dataset, Eigen::MatrixXd& normalize );
    int get_HDF5_mean_sd_by_row_ptr(H5File* file, DataSet* dataset, Eigen::MatrixXd& normalize );
    
    // functions called from R
    void bdgetSDandMean_hdf5( std::string filename, const std::string group, std::string dataset,
                              Rcpp::Nullable<bool> sd, Rcpp::Nullable<bool> mean,
                              Rcpp::Nullable<bool> byrows,
                              Rcpp::Nullable<int> wsize, Rcpp::Nullable<int> force);

#endif