#ifndef hdf5_imputation
#define hdf5_imputation

#include <RcppEigen.h>
#include "H5Cpp.h"
#include "rhdf5Utils.h"
#include<random>



// functions called from c++
int get_value_to_impute_discrete(std::map<double, double> probMap);
std::map<double, double> VectortoOrderedMap_SNP_counts( Eigen::VectorXd  vdata);
void Impute_snp_HDF5(H5File* file, DataSet* dataset, bool bycols, std::string stroutdataset);



// functions called from R
/**
 void bdImputeSNPHDF5(std::string filename, std::string group, std::string dataset, Rcpp::Nullable<std::string> outgroup, 
 Rcpp::Nullable<std::string> outdataset, Rcpp::Nullable<bool> bycols = true );
 **/
#endif