#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixPCA.hpp"
#include "Utilities/Utilities.hpp"


//' PCA Descomposition
//' 
//' Compute PCA
//' 
//' @param filename string, file name where dataset is stored 
//' @param group string group name  where dataset is stored in file
//' @param dataset string dataset name with data to perform PCA
//' @param ncomponents integer, number of components to be computed, by default 
//' ncomponents = 0, all components are computed
//' @param bcenter logical value if true data is centered to zero
//' @param bscale logical value, if true data is scaled
//' @param k number of local SVDs to concatenate at each level, performance parameter 
//' @param q number of levels to compute SVD for PCA, performance parameter
//' @param rankthreshold double, threshold used to determine the range of the array. 
//' The matrix rank is equal to the number of
//' singular values different from the threshold. By default, threshold = 0 is used 
//' to get the matrix rank , but it can be changed to an approximation of 0.
//' @param SVDgroup string. Name of the group where the SVD results are located. 
//' If it has been previously calculated. This group must contain the d, u and v datasets.
//' @param force logical value, if true, the SVD is forced to be computed although 
//' the SVD exists. 
//' @param method (optional, defalut = "auto") possible values are: "auto", 
//' "blocks", "full":
//' 
//'   * `"auto"`:
//'     The option method = "auto" chooses the "full" or "blocks" method depending on 
//'     the size of the matrix to be decomposed }
//'   * `"blocks"`:
//'     The PCA can be carried out by blocks, recommended option for large matrices 
//'     that do not fit in memory
//'   * `"full"`:
//'     The PCA is performed directly without partitioning the matrix 
//' 
//' @param threads integer number of threads used to run PCA
//' @return original file with results in folder PCA/<datasetname>
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
        // pcaeig pcaRes;
        
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
