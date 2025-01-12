#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixPCA.hpp"
#include "Utilities/Utilities.hpp"


//' PCA Descomposition with HDF5 files
//' 
//' Calculates the Principal Component Analysis (PCA) for datasets stored in HDF5 format.
//' 
//' @inheritParams bdNormalize_hdf5
//' @param ncomponents integer, An optional integer specifying the number of 
//' principal components to calculate. Defaults to 0, which computes all components.
//' @param bcenter logical (optional). If TRUE, the data is centered by 
//' subtracting the column means (ignoring NAs) of `x` from their corresponding columns. 
//' If FALSE (default), no centering is performed.
//' @param bscale (optional). If TRUE, the data is scaled by dividing 
//' the (centered) columns of `x` by their standard deviations if `bcenter` 
//' is TRUE, or by the root mean square otherwise. If FALSE (default), no 
//' scaling is performed.
//' @param SVDgroup string. Name of the group where the intermediate SVD results 
//' will be stored or located (If it has been previously calculated). 
//' If falready calculated, this group must contain the d, u and v datasets.
//' @param method optional, The method argument can specify different PCA algorithms, 
//' defalut = "auto" possible values are: "auto", "blocks", "full":
//' 
//'   * `"auto"`:
//'     The option method = "auto" chooses the "full" or "blocks" method depending on 
//'     the size of the matrix to be decomposed
//'   * `"blocks"`:
//'     The PCA can be carried out by blocks, recommended option for large matrices 
//'     that do not fit in memory
//'   * `"full"`:
//'     The PCA is performed directly without partitioning the matrix 
//' 
//' @return original file with results in folder PCA/\<datasetname\>
//' @export
// [[Rcpp::export]]
void bdPCA_hdf5(std::string filename, std::string group, std::string dataset,
                Rcpp::Nullable<int> ncomponents = 0,
                Rcpp::Nullable<bool> bcenter = false, Rcpp::Nullable<bool> bscale = false, 
                Rcpp::Nullable<int> k=2, Rcpp::Nullable<int> q=1,
                Rcpp::Nullable<double> rankthreshold = 0.0,
                Rcpp::Nullable<std::string> SVDgroup = R_NilValue,
                Rcpp::Nullable<bool> overwrite = false, 
                Rcpp::Nullable<Rcpp::CharacterVector> method = R_NilValue,
                Rcpp::Nullable<int> threads = R_NilValue)
{
    

    try
    {
        
        bool bcent, bscal, bforce;
        int ks, qs = 1, incomponents = 0;
        double dthreshold;
        std::string strSVDgroup;
        
        if(ncomponents.isNull())  incomponents = 0 ;
        else    incomponents = Rcpp::as<int>(ncomponents);
        
        if(q.isNull())  qs = 1 ;
        else    qs = Rcpp::as<int>(q);
        
        if(k.isNull())  ks = 2 ;
        else    ks = Rcpp::as<int>(k);
        
        if(bcenter.isNull())  bcent = false ;
        else    bcent = Rcpp::as<bool>(bcenter);
        
        if(bscale.isNull())  bscal = false ;
        else    bscal = Rcpp::as<bool>(bscale);
        
        if(overwrite.isNull())  bforce = false ;
        else    bforce = Rcpp::as<bool>(overwrite);
        
        if(rankthreshold.isNull()) {  
            dthreshold = 0 ;
        } else {
            if( Rcpp::as<double>(rankthreshold) > 0.1 ) {
                Rcpp::Rcout<< "Threshold to big, please set threshold with value lower than 0.1";
                return void();
            } else if( Rcpp::as<double>(rankthreshold) < 0 ) {
                Rcpp::Rcout<< "Threshold must be a positive value near zero";
                return void();
            } else {
                dthreshold = Rcpp::as<double>(rankthreshold);
            }
        }
        
        if(SVDgroup.isNull()) { strSVDgroup = "SVD"; } 
        else {   strSVDgroup = Rcpp::as<std::string>(SVDgroup); }
        
        if( strSVDgroup.substr(strSVDgroup.length(), strSVDgroup.length()) != "/" ){
            strSVDgroup = strSVDgroup + "/";
        }
        
        BigDataStatMeth::RcppPCAHdf5(filename, group, dataset, strSVDgroup, ks, qs, incomponents, bcent, bscal, dthreshold, bforce, false, method, threads);
        
    } catch( H5::FileIException& error ) {
        ::Rf_error( "c++ exception bdPCA_hdf5 (File IException)" );
        return void();
    } catch( H5::DataSetIException& error ) { 
        ::Rf_error( "c++ exception bdPCA_hdf5 (DataSet IException)" );
        return void();
    } catch( H5::DataSpaceIException& error ) { 
        ::Rf_error( "c++ exception bdPCA_hdf5 (DataSpace IException)" );
        return void();
    } catch( H5::DataTypeIException& error ) {
        ::Rf_error( "c++ exception bdPCA_hdf5 (DataType IException)" );
        return void();
    } catch(std::exception &ex) {
        Rcpp::Rcout<< "C++ exception bdPCA_hdf5 : "<< ex.what();
        return void();
    } catch (...) {
        ::Rf_error("C++ exception bdPCA_hdf5 (unknown reason)");
        return void();
    }
    
    return void();
    
}

 /***
  //' @param ncomponents integer, An optional integer specifying the number of 
  //' principal components to calculate. Defaults to 0, which computes all components.
  //' @param bcenter logical (optional). If TRUE, the data is centered by 
  //' subtracting the column means (ignoring NAs) of `x` from their corresponding columns. 
  //' If FALSE (default), no centering is performed.
  //' @param bscale (optional). If TRUE, the data is scaled by dividing 
  //' the (centered) columns of `x` by their standard deviations if `bcenter` 
  //' is TRUE, or by the root mean square otherwise. If FALSE (default), no 
  //' scaling is performed.
  //' @param k number of local SVDs to concatenate at each level. Defaults to 2.
  //' This parameter helps optimize the performance and memory usage during PCA 
  //' calculations. 
  //' @param q number of levels to compute SVD for PCA.
  //' This parameter helps optimize the performance and memory usage during PCA 
  //' calculations. 
  //' @param rankthreshold `double`. Threshold used to determine the range of the matrix. 
  //' The matrix rank is defined as the number of singular values that differ from the 
  //' threshold. By default, `threshold = 0` is used to compute the matrix rank, 
  //' but it can be adjusted to a value close to zero for approximations.
  //' @param SVDgroup string. Name of the group where the intermediate SVD results 
  //' will be stored or located (If it has been previously calculated). 
  //' If falready calculated, this group must contain the d, u and v datasets.
  //' @param overwrite logical value, if true, the SVD is forced to be computed 
  //' although the SVD results exists. 
  //' @param method optional, The method argument can specify different PCA algorithms, 
  //' defalut = "auto" possible values are: "auto", "blocks", "full":
  //' 
  //'   * `"auto"`:
  //'     The option method = "auto" chooses the "full" or "blocks" method depending on 
  //'     the size of the matrix to be decomposed
  //'   * `"blocks"`:
  //'     The PCA can be carried out by blocks, recommended option for large matrices 
  //'     that do not fit in memory
  //'   * `"full"`:
  //'     The PCA is performed directly without partitioning the matrix 
  //' 
  //' @param threads integer (optional), an optional parameter specifying the 
  //' number of threads to use.
  //' @return original file with results in folder PCA/\<datasetname\>
  */