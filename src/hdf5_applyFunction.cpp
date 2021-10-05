#include "include/hdf5_applyFunction.h"


// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;


//' Apply function to different datasets inside a group
//'
//' Apply function to different datasets inside a group
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param datasets, character array indicating the input dataset to be imputed
//' @param func, character array function to be applyed
//' \describe{
//'     \item{QR}{apply bdQR() function to datasets}
//'     \item{CrossProd}{apply bdCrossprod() function to datasets}
//'     \item{tCrossProd}{apply bdtCrossprod() function to datasets}
//'     \item{invChol}{apply bdInvCholesky() function to datasets}
//' }
//' @param outgroup, optional character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param bforce, boolean if true, previous results in same location inside hdf5 will be overwritten.
//' @return Original hdf5 data file with results after apply function to different datasets
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdapply_Function_hdf5( std::string filename, std::string group, Rcpp::StringVector datasets, 
                                     std::string outgroup, std::string func,
                                     Rcpp::Nullable<bool> force = false )
{
    
    H5File* file;
    DataSet* pdataset = nullptr;
    Rcpp::NumericVector oper = {0, 1, 2, 3};
    oper.names() = Rcpp::CharacterVector({"QR", "CrossProd", "tCrossProd", "invChol"});
    
    try
    {

        bool bforce;
        
        if(force.isNull()) { bforce = false; } 
        else {   bforce = Rcpp::as<bool>(force); }

        
        // Test file
        if( ResFileExist_filestream(filename) ) {
            file = new H5File( filename, H5F_ACC_RDWR ); 
        } else {
            Rcpp::Rcout<<"\nFile not exits, create file before split dataset";
            return wrap(false);
        }
        
        
        // Seek all datasets to perform calculus
        for( int i=0; i < datasets.size(); i++ ) 
        {
            std::string strdataset = group +"/" + datasets(i);
            
            if( exists_HDF5_element_ptr(file, strdataset ) == 0 ) {
                file->close();
                Rcpp::Rcout<<"Group not exists, create the input dataset before proceed";
                return wrap(false);
            }
            
            pdataset = new DataSet(file->openDataSet(strdataset));
            
            // Real data set dimension
            IntegerVector dims_out = get_HDF5_dataset_size(*pdataset);

            // Get block from complete matrix
            Eigen::MatrixXd original = GetCurrentBlock_hdf5_Original( file, pdataset, 0, 0, dims_out[0], dims_out[1]);
            
            if(i==0) {
                prepare_outGroup(file, outgroup, bforce);
            }
            
            
            if( oper.findName( func ) == 0)
            {
                strQR decQR;
                decQR = rcpp_bdQR(original, true);
                
                write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i) + ".Q", Rcpp::wrap(decQR.Q), false);
                write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i) + ".R", Rcpp::wrap(decQR.R), false);
                
                pdataset->close();
                
            }else if( oper.findName( func ) == 1) {
                
                Eigen::MatrixXd results = bdcrossproduct(original);    
                write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i), Rcpp::wrap(results), false);
                pdataset->close();
                
            }else if( oper.findName( func ) == 2) {
                
                Eigen::MatrixXd results = bdtcrossproduct(original);    
                write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i), Rcpp::wrap(results), false);
                pdataset->close();
                
            }else if( oper.findName( func ) == 3) {
                
                svdeig results = RcppCholDec(original);    
                if( results.v == Eigen::MatrixXd::Zero(2,2) && results.u == Eigen::MatrixXd::Zero(2,2)) {
                    pdataset->close();
                    file->close();
                    return wrap(false);
                    
                } else {
                    write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i), Rcpp::wrap(results.v), false);
                    pdataset->close();
                }
                
            }else {
                pdataset->close();
                file->close();
                Rcpp::Rcout<<"Group not exists, create the input dataset before proceed";
                return wrap(false);
                
            }
            
        }
        
        
    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdataset->close();
        file->close();
        ::Rf_error( "c++ exception (File IException)" );
        return(wrap(-1));
    }
    
    file->close();
    
    return(wrap(0));
}



/***R

library(BigDataStatMeth)

setwd("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/UAB/DOCTORAT/BitDataStatMeth - BDSM/Analysis/BigDataStatMeth_Analysis/Cholesterol/test")

bdSplit_matrix_hdf5( "cars.hdf5", "data", "X", "dataoutCols", nblocks = 3, bycols = FALSE, force = TRUE)

bdSplit_matrix_hdf5( "cars.hdf5", "data", "X", "dataoutRows", nblocks = 3, bycols = TRUE, force = TRUE)



*/