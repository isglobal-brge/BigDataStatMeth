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
//' @param datasets, character array indicating the input datasets to be used
//' @param func, character array function to be applyed
//' \describe{
//'     \item{QR}{apply bdQR() function to datasets}
//'     \item{CrossProd}{apply bdCrossprod() function to datasets}
//'     \item{tCrossProd}{apply bdtCrossprod() function to datasets}
//'     \item{invChol}{apply bdInvCholesky() function to datasets}
//'     \item{blockmult}{apply matrix multiplication, in that case, we need the datasets to be used defined
//'     in b_datasets variable, datasets and b_datasets must be of the same lenght, in that case, the operation is performed according to index, for example,
//'     if we have datasets = {"A1", "A2", "A3} and b_datasets = {"B1", "B2", "B3}, the functions performs : A1%*%B1, A2%*%B2 and A3%*%B3 }
//' }
//' @param outgroup, character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param b_datasets, optional character array indicating the input datasets to be used when we need a second dataset in functions like matrix multiplication
//' @param bforce, optional boolean if true, previous results in same location inside hdf5 will be overwritten, by default force = false, data was not overwritten..
//' @return Original hdf5 data file with results after apply function to different datasets
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdapply_Function_hdf5( std::string filename, 
                                     std::string group, 
                                     Rcpp::StringVector datasets, 
                                     std::string outgroup, 
                                     std::string func, 
                                     Rcpp::Nullable<std::string> b_group = R_NilValue, 
                                     Rcpp::Nullable<Rcpp::StringVector> b_datasets = R_NilValue,
                                     Rcpp::Nullable<bool> force = false )
{
    
    H5File* file;
    DataSet* pdataset = nullptr;
    DataSet* pbdataset = nullptr;
    Rcpp::StringVector str_bdatasets;
    std::string str_bgroup;
    Rcpp::NumericVector oper = {0, 1, 2, 3, 4, 11, 22};
    oper.names() = Rcpp::CharacterVector({"QR", "CrossProd", "tCrossProd", "invChol", "blockmult", "CrossProd_double", "tCrossProd_double"});
    
    try
    {

        bool bforce;
        
        if(force.isNull()) { bforce = false; } 
        else {   bforce = Rcpp::as<bool>(force); }

        // Test file
        if( ResFileExist_filestream(filename) ) {
            file = new H5File( filename, H5F_ACC_RDWR ); 
        } else {
            Rcpp::Rcout<<"\nFile not exits, create file before apply function to datasets";
            return wrap(false);
        }


        if( b_datasets.isNotNull() &&  ( oper(oper.findName( func )) == 1 ||  oper(oper.findName( func )) == 2 ||  oper(oper.findName( func )) == 4) ) {
            
            if( as<Rcpp::StringVector>(b_datasets).size() != datasets.size() ){
                Rcpp::Rcout<<"To perform matrix multiplication, CrossProd or tCrossProd "<<
                    "with two matrices b_datasets variable must be defined and the length "<<
                        " of datasets and b_datasets must be equal";
                return wrap(false);  
            }
            str_bdatasets = as<Rcpp::StringVector>(b_datasets);
            
            if( oper.findName( func ) == 1){
                func = "CrossProd_double";
            }else if( oper.findName( func ) == 2){
                func = "tCrossProd_double";
            }
            
        }
        
        

        if(b_group.isNull()) { str_bgroup = group; } 
        else {   str_bgroup = Rcpp::as<std::string>(b_group); }

        // Seek all datasets to perform calculus
        for( int i=0; i < datasets.size(); i++ ) 
        {
            std::string strdataset = group +"/" + datasets(i);
            
            if( exists_HDF5_element_ptr(file, strdataset ) == 0 ) {
                
                Rcpp::Rcout<<"\n Estem buscant ... :"<<strdataset<<"\n";
                
                file->close();
                Rcpp::Rcout<<"Group or dataset does not exists, create the input dataset before proceed";
                return wrap(false);
            }
            
            pdataset = new DataSet(file->openDataSet(strdataset));
            
            // Real data set dimension
            IntegerVector dims_out = get_HDF5_dataset_size(*pdataset);

            // Get block from complete matrix
            Eigen::MatrixXd original = GetCurrentBlock_hdf5_Original( file, pdataset, 0, 0, dims_out[0], dims_out[1]);
            
            if(i==0) {
                prepare_outGroup(file, outgroup, bforce);
                if(oper.findName( func ) == 0){
                    prepare_outDataset(file, outgroup + "/" + datasets(i) + ".Q", bforce);
                    prepare_outDataset(file, outgroup + "/" + datasets(i) + ".R", bforce);
                }else {
                    prepare_outDataset(file, outgroup + "/" + datasets(i), bforce);
                }
            }
            
            
            if( oper(oper.findName( func )) == 0)
            {
                strQR decQR;
                decQR = rcpp_bdQR(original, true);
                
                write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i) + ".Q", Rcpp::wrap(decQR.Q), false);
                write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i) + ".R", Rcpp::wrap(decQR.R), false);
                
                pdataset->close();
                
            } else if( oper(oper.findName( func )) == 1) {
                
                Eigen::MatrixXd results = bdcrossproduct(original);    
                write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i), Rcpp::wrap(results), false);
                pdataset->close();
                
            } else if( oper(oper.findName( func )) == 2) {
                
                Eigen::MatrixXd results = bdtcrossproduct(original);    
                write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i), Rcpp::wrap(results), false);
                pdataset->close();
                
            } else if( oper(oper.findName( func )) == 3) {
                
                svdeig results = RcppCholDec(original);    
                if( results.v == Eigen::MatrixXd::Zero(2,2) && results.u == Eigen::MatrixXd::Zero(2,2)) {
                    pdataset->close();
                    file->close();
                    return wrap(false);
                    
                } else {
                    write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i), Rcpp::wrap(results.v), false);
                    pdataset->close();
                }
                
            } else if( oper(oper.findName( func )) == 4 || oper(oper.findName( func )) == 11 || oper(oper.findName( func )) == 22) {
                
                std::string outputdataset;
                Eigen::MatrixXd originalB;
                    
                std::string b_strdataset = str_bgroup + "/" + str_bdatasets(i);
                
                if( exists_HDF5_element_ptr(file, b_strdataset ) == 0 ) {
                    
                    Rcpp::Rcout<<"\n Estem buscant ... :"<<b_strdataset<<"\n";
                    
                    pdataset->close();
                    file->close();
                    Rcpp::Rcout<<"Group or dataset does not exists, create the input dataset before proceed";
                    return wrap(false);
                }

                pbdataset = new DataSet(file->openDataSet(b_strdataset));
                
                // Real data set dimension
                IntegerVector dims_outB = get_HDF5_dataset_size(*pbdataset);

                originalB = GetCurrentBlock_hdf5_Original( file, pbdataset, 0, 0, dims_outB[0], dims_outB[1]);
                
                // 
                // Rcpp::Rcout<<"\n ORIGINAL : \n"<<original<<"\n";
                // Rcpp::Rcout<<"\n ORIGINAL B : \n"<<originalB<<"\n";
                
                
                if( oper(oper.findName( func )) == 4 ) {
                    outputdataset = outgroup + "/" + datasets(i) + "_x_" + str_bdatasets(i);
                } else if  (oper(oper.findName( func )) == 11) {
                    outputdataset = outgroup + "/Cross_" + datasets(i) + str_bdatasets(i);
                    original = GetCurrentBlock_hdf5( file, pdataset, 0, 0, dims_out[0], dims_out[1]);
                    
                } else if ( oper(oper.findName( func )) == 22) {
                    outputdataset = outgroup + "/tCross_" + datasets(i) + "_x_" + str_bdatasets(i);
                    originalB = GetCurrentBlock_hdf5( file, pbdataset, 0, 0, dims_outB[0], dims_outB[1]);
                }
                
                // Rcpp::Rcout<<"\n  \n \n";
                // 
                // Rcpp::Rcout<<"\n ORIGINAL : \n"<<original<<"\n";
                // Rcpp::Rcout<<"\n ORIGINAL B : \n"<<originalB<<"\n";
                
                // // Get data
                // originalB = GetCurrentBlock_hdf5_Original( file, pbdataset, 0, 0, dims_outB[0], dims_outB[1]);

                
                prepare_outGroup(file, outgroup, bforce);
                prepare_outDataset(file, outputdataset, bforce);
                
                
                Eigen::MatrixXd results = Bblock_matrix_mul_parallel(original, originalB, 128, R_NilValue);
                
                write_HDF5_matrix_from_R_ptr(file, outputdataset, Rcpp::wrap(results), false);
                
                pdataset->close();
                pbdataset->close();

                
            } else {
                pdataset->close();
                file->close();
                Rcpp::Rcout<<"Function does not exists, please use one of the following : 'QR', 'CrossProd', 'tCrossProd', 'invChol', 'blockmult' ";
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