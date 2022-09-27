#include "include/hdf5_applyFunction.h"


// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;


//' Apply function to different datasets inside a group
//'
//' Apply function to different datasets inside a group
//' 
//' @param filename, Character array, indicating the name of the file to create
//' @param group, Character array, indicating the input group where the data set
//' to be imputed is. 
//' @param datasets, Character array, indicating the input datasets to be used
//' @param outgroup, Character, array, indicating group where the data set will 
//' be saved after imputation if `outgroup` is NULL, output dataset is stored 
//' in the same input group. 
//' @param func, Character array, function to be applyed : 
//' QR to apply bdQR() function to datasets
//' CrossProd to apply bdCrossprod() function to datasets
//' tCrossProd to apply bdtCrossprod() function to datasets
//' invChol to apply bdInvCholesky() function to datasets
//' blockmult to apply matrix multiplication, in that case, we need the datasets 
//' to be used defined in b_datasets variable, datasets and b_datasets must be 
//' of the same lenght, in that case, the operation is performed according to 
//' index, for example, if we have datasets = {"A1", "A2", "A3} and 
//' b_datasets = {"B1", "B2", "B3}, the functions performs : A1%*%B1, 
//' A2%*%B2 and A3%*%B3 
//' CrossProd_double to  performs crossprod using two matrices, see blockmult 
//' tCrossProd_double to  performs transposed crossprod using two matrices, 
//' see blockmult 
//' solve to solve matrix equation system, see blockmult for parametrization 
//' sdmean to get sd and mean from de datasets by cols or rows
//' @param b_group, optional Character array indicating the input group where 
//' data are stored when we need a second dataset to operate, for example in 
//' functions like matrix multiplication
//' @param b_datasets, optional Character array indicating the input datasets 
//' to be used when we need a second dataset in functions like matrix 
//' multiplication
//' @param force, optional Boolean if true, previous results in same location 
//' inside hdf5 will be overwritten, by default force = false, data was not 
//' overwritten.
//' @param transp_dataset optional parameter. Boolean if true we use the 
//' transposed dataframe to perform calculus. By default transp_dataset = false, 
//' we use the original dataset stored in hdf5 data file. Currently this option 
//' is only valid with "blockmult", "CrossProd_double" and "tCrossProd_double"
//' @param transp_bdataset optional parameter. Boolean if true we use the 
//' transposed dataframe to perform calculus.By default transp_bdataset = false, 
//' we use the original dataset stored in hdf5 data file. Currently this option 
//' is only valid with "blockmult", "CrossProd_double" and "tCrossProd_double"
//' @param fullMatrix boolean, optional parameter used in Inverse Cholesky, by 
//' default false. If fullMatrix = true, in the hdf5 file the complete matrix 
//' is stored. If false, only the lower triangular matrix is saved
//' @param byrows boolean, optional parameter used in sd and mean calculus, by 
//' default false. If byrows = true, the sd and mean is computed by columns. 
//' If false, sd and mean is computed by rows.
//' @param threads optional parameter. Integer with numbers of threads to be used
//' @return Original hdf5 data file with results after apply function to 
//' different datasets
//' @export
// [[Rcpp::export]]
void bdapply_Function_hdf5( std::string filename, 
                                     std::string group, 
                                     Rcpp::StringVector datasets, 
                                     std::string outgroup, 
                                     std::string func, 
                                     Rcpp::Nullable<std::string> b_group = R_NilValue, 
                                     Rcpp::Nullable<Rcpp::StringVector> b_datasets = R_NilValue,
                                     Rcpp::Nullable<bool> force = false,
                                     Rcpp::Nullable<bool> transp_dataset = false,
                                     Rcpp::Nullable<bool> transp_bdataset = false,
                                     Rcpp::Nullable<bool> fullMatrix = false,
                                     Rcpp::Nullable<bool> byrows = false,
                                     Rcpp::Nullable<int> threads = 2 )
{
    
    H5File* file = nullptr;
    DataSet* pdataset = nullptr;
    DataSet* pbdataset = nullptr;
    Rcpp::StringVector str_bdatasets;
    std::string str_bgroup;
    bool btransdataA, btransdataB, bfullMatrix, bbyrows;
    Rcpp::NumericVector oper = {0, 1, 2, 3, 4, 11, 22, 5, 6, 7};
    oper.names() = Rcpp::CharacterVector({"QR", "CrossProd", "tCrossProd",
               "invChol", "blockmult", "CrossProd_double", "tCrossProd_double",
               "solve", "normalize", "sdmean"});

    try
    {

        bool bforce;
        
        if(force.isNull()) { bforce = false; } 
        else {   bforce = Rcpp::as<bool>(force); }

        if(transp_dataset.isNull()) { btransdataA = false; } 
        else {   btransdataA = Rcpp::as<bool>(transp_dataset); }

        if(transp_bdataset.isNull()) { btransdataB = false; } 
        else {   btransdataB = Rcpp::as<bool>(transp_bdataset); }
        
        if(fullMatrix.isNull()) { bfullMatrix = false; } 
        else {   bfullMatrix = Rcpp::as<bool>(fullMatrix); }
        
        if(byrows.isNull()) { bbyrows = false; } 
        else {   bbyrows = Rcpp::as<bool>(byrows); }
        
        // Test file
        if( ResFileExist_filestream(filename) ) {
            file = new H5File( filename, H5F_ACC_RDWR ); 
        } else {
            Rcpp::Rcout<<"\nFile not exits, create file before apply function to datasets";
            return void();
            // return wrap(false);
        }

        if( b_datasets.isNotNull() &&  ( oper(oper.findName(func)) == 1 ||  
           oper(oper.findName(func)) == 2 || oper(oper.findName(func)) == 4 || 
           oper(oper.findName(func)) == 5 || oper(oper.findName(func)) == 6) ) {
            
            //. 01/01/2022 . // if( as<Rcpp::StringVector>(b_datasets).size() != datasets.size() ){
            //. 01/01/2022 . //      Rcpp::Rcout<<"To perform matrix multiplication, CrossProd, tCrossProd or solve "<<
            //. 01/01/2022 . //          "with two matrices b_datasets variable must be defined and the length "<<
            //. 01/01/2022 . //              " of datasets and b_datasets must be equal";
            //. 01/01/2022 . //      return void();
            // return wrap(false);  
            // }
            str_bdatasets = as<Rcpp::StringVector>(b_datasets);
            
            if( oper.findName( func ) == 1) {
                func = "CrossProd_double";
            } else if( oper.findName( func ) == 2) {
                func = "tCrossProd_double";
            }
            
        }
        
        if( b_group.isNull()) { 
            str_bgroup = group; 
        } else {   
            str_bgroup = Rcpp::as<std::string>(b_group); 
        }
        
        // Seek all datasets to perform calculus
        for( int i=0; i < datasets.size(); i++ ) 
        {
            
            std::string strdataset = group +"/" + datasets(i);
            
            if( exists_HDF5_element_ptr(file, strdataset ) == 0 ) {

                file->close();
                Rcpp::Rcout<<"Group or dataset does not exists, create the input dataset before proceed";
                return void();
                // return wrap(false);
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
                } else if (oper.findName( func ) == 5) {
                    std::string tmp_strdataset = str_bgroup + "/" + str_bdatasets(i);
                    std::string tmp_outputdataset = outgroup + "/solved_" + datasets(i) + "x_eq_" + str_bdatasets(i);
                    prepare_outDataset(file, tmp_outputdataset, bforce);
                } else if (oper.findName( func ) == 7) {
                    prepare_outDataset(file, outgroup + "/mean." + datasets(i), bforce);
                    prepare_outDataset(file, outgroup + "/sd." + datasets(i), bforce);
                } else {
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
                int ithreads;
                
                Rcpp::Nullable<long> elementsBlock = R_NilValue;
                    
                Rcpp_bdInvCholesky_hdf5(file, pdataset, 
                                        outgroup, Rcpp::as<std::string>(datasets[i]), 
                                        bforce, bfullMatrix, threads, elementsBlock );
                    pdataset->close();
                
                // svdeig results = RcppCholDec(original);    
                // if( results.v == Eigen::MatrixXd::Zero(2,2) && results.u == Eigen::MatrixXd::Zero(2,2)) {
                //     pdataset->close();
                //     file->close();
                //     return void();
                //     // return wrap(false);
                //     
                // } else {
                //     write_HDF5_matrix_from_R_ptr(file, outgroup + "/" + datasets(i), Rcpp::wrap(results.v), false);
                //     pdataset->close();
                // }
                
                // bdInvCholesky_hdf5(filename, group, dataset, outgroup, "InverseA", elementsBlock = 100000),
                
            } else if( oper(oper.findName( func )) == 4 || 
                       oper(oper.findName( func )) == 11 || 
                       oper(oper.findName( func )) == 22) {
                
                std::string outputdataset;
                Eigen::MatrixXd originalB;
                    
                std::string b_strdataset = str_bgroup + "/" + str_bdatasets(i);
                
                if( exists_HDF5_element_ptr(file, b_strdataset ) == 0 ) {
                    
                    pdataset->close();
                    file->close();
                    Rcpp::Rcout<<"Group or dataset does not exists, create the input dataset before proceed";
                    return void();
                    // return wrap(false);
                }
                
                pbdataset = new DataSet(file->openDataSet(b_strdataset));

                // Real data set dimension
                IntegerVector dims_outB = get_HDF5_dataset_size(*pbdataset);

                originalB = GetCurrentBlock_hdf5_Original( file, pbdataset, 0, 0, dims_outB[0], dims_outB[1]);

                if( oper(oper.findName( func )) == 4 ) {

                    //. Changed Criteria - names too long ---> Not portable .// 
                    //      outputdataset = outgroup + "/" + datasets(i) + "_x_" + str_bdatasets(i);
                    outputdataset = outgroup + "/" + datasets(i);
                    
                } else if  (oper(oper.findName( func )) == 11) {
                    outputdataset = outgroup + "/Cross_" + datasets(i) + str_bdatasets(i);
                    original = GetCurrentBlock_hdf5( file, pdataset, 0, 0, dims_out[0], dims_out[1]);
                    
                } else if ( oper(oper.findName( func )) == 22) {
                    outputdataset = outgroup + "/tCross_" + datasets(i) + str_bdatasets(i);
                    originalB = GetCurrentBlock_hdf5( file, pbdataset, 0, 0, dims_outB[0], dims_outB[1]);
                }

                prepare_outGroup(file, outgroup, bforce);
                prepare_outDataset(file, outputdataset, bforce);

                Eigen::MatrixXd results;
                
                if (btransdataA == false && btransdataB == false) {
                    results = Bblock_matrix_mul_parallel(original, originalB, 2048, R_NilValue);
                } else if (btransdataA == true && btransdataB == false) {
                    results = Bblock_matrix_mul_parallel(original.transpose(), originalB, 2048, R_NilValue);
                }else if (btransdataA == false && btransdataB == true) {
                    results = Bblock_matrix_mul_parallel(original, originalB.transpose(), 2048, R_NilValue);
                } else {
                    results = Bblock_matrix_mul_parallel(original.transpose(), originalB.transpose(), 2048, R_NilValue);
                }

                write_HDF5_matrix_from_R_ptr(file, outputdataset, Rcpp::wrap(results), false);

                pdataset->close();
                pbdataset->close();

            } else if( oper(oper.findName( func )) == 5) {
                
                std::string outputdataset;
                Eigen::MatrixXd originalB;
                
                
                std::string b_strdataset = str_bgroup + "/" + str_bdatasets(i);
                outputdataset = outgroup + "/solved_" + datasets(i) + "x_eq_" + str_bdatasets(i);
                
                pbdataset = new DataSet(file->openDataSet(b_strdataset));
                
                // Real data set dimension
                IntegerVector dims_outB = get_HDF5_dataset_size(*pbdataset);
                
                original = GetCurrentBlock_hdf5( file, pdataset, 0, 0, dims_out[0], dims_out[1]);
                originalB = GetCurrentBlock_hdf5( file, pbdataset, 0, 0, dims_outB[0], dims_outB[1]);
                
                Rcpp::NumericMatrix results = Rcpp::as<Rcpp::NumericMatrix>(bdSolve(wrap(original), wrap(originalB)));
                
                write_HDF5_matrix_from_R_ptr(file, outputdataset, results, false);
                pdataset->close();
                
            } else if( oper(oper.findName( func )) == 7) {
                
                int blocksize;
                
                // std::string strgroupout = "mean_sd/" + group;
                // prepare_outGroup(file, outgroup, bforce);

                IntegerVector dims_out = get_HDF5_dataset_size_ptr(pdataset);
                Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2,dims_out[0]);
                
                if( bbyrows == false) {
                    
                    if(dims_out[1] > maxElemBlock){
                        blocksize = 1;
                    } else {
                        int maxsize = std::max( dims_out[0], dims_out[1]);
                        blocksize = std::ceil( maxElemBlock / maxsize); }
                    
                    datanormal = Eigen::MatrixXd::Zero(2,dims_out[0]);
                    get_HDF5_mean_sd_by_column_ptr( file, pdataset, datanormal);
                } else {
                    
                    if(dims_out[0] > maxElemBlock) {
                        blocksize = 1;
                    } else {
                        int maxsize = std::max( dims_out[0], dims_out[1]);
                        blocksize = std::ceil( maxElemBlock / maxsize); }
                    
                    datanormal = Eigen::MatrixXd::Zero(2,dims_out[1]);
                    get_HDF5_mean_sd_by_row_ptr( file, pdataset, datanormal);
                }
                
                // if(exists_HDF5_element_ptr(file, strgroupout) == 0) {
                //     create_HDF5_groups_ptr( file, strgroupout); }
                
                // Store center and scale for each column
                write_HDF5_matrix_ptr(file, outgroup + "/mean."+ datasets(i), wrap(datanormal.row(0)));
                write_HDF5_matrix_ptr(file, outgroup + "/sd."+ datasets(i), wrap(datanormal.row(1)));
                pdataset->close();

            } else {
                pdataset->close();
                file->close();
                Rcpp::Rcout<<"Function does not exists, please use one of the following : 'QR', 'CrossProd', 'tCrossProd', 'invChol', 'blockmult' ";
                return void();
                // return wrap(false);
                
            }
            
        }
        
        
    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdataset->close();
        pbdataset->close();
        file->close();
        ::Rf_error( "c++ exception (File IException)" );
        return void();
        // return(wrap(-1));
    }
    
    file->close();
  
    Rcpp::Rcout<< func <<" function has been computed in all blocks\n";  
    return void();
    // return(wrap(0));
}



/***R

library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))
setwd("/Users/mailos/DOCTORAT_Local/BigDataStatMeth/")

# Prepare data and functions
set.seed(123)
Y <- matrix(rnorm(250), 10, 10)
X <- matrix(rnorm(250), 20, 20)
filename <- "test.hdf5"

# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("test.hdf5", Y, "data", "Y", force = T)
bdAdd_hdf5_matrix( X, "test.hdf5",  "data", "X", force = TRUE)


data <- bdgetDatasetsList_hdf5("test.hdf5", group = "data")

bdapply_Function_hdf5(filename = filename,
                      group = "data",datasets = data,
                      outgroup = "KXA",func = "sdmean",
                      force = T)




*/
