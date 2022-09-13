#include "include/hdf5_computeMatrixVector.h"

// Internal calls :

// by Rows

Eigen::MatrixXd Rcpp_matrixVectorMultiplication_byRow(Eigen::MatrixXd X, Eigen::VectorXd v) {
    
    X = X.array().colwise() * v.array();
    return(X);
}

Eigen::MatrixXd Rcpp_matrixVectorSum_byRow(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().colwise() + v.array();
    return(X);
}

Eigen::MatrixXd Rcpp_matrixVectorSubstract_byRow(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().colwise() - v.array();
    return(X);
}

Eigen::MatrixXd Rcpp_matrixVectorDivision_byRow(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().colwise() / v.array();
    return(X);
}


// by Columns

Eigen::MatrixXd Rcpp_matrixVectorMultiplication_byCol(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().rowwise() * v.transpose().array();    
    return(X);
}

Eigen::MatrixXd Rcpp_matrixVectorSum_byCol(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().rowwise() + v.transpose().array();    
    return(X);
}

Eigen::MatrixXd Rcpp_matrixVectorSubstract_byCol(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().rowwise() - v.transpose().array();    
    return(X);
}

Eigen::MatrixXd Rcpp_matrixVectorDivision_byCol(Eigen::MatrixXd X, Eigen::VectorXd v) {
    X = X.array().rowwise() / v.transpose().array();    
    return(X);
}


//' Apply vector calculus to a dataset in hdf5 file
//' 
//' This function applies a calculus with a vector to a matrix. Multiplies, 
//' sums, substract or divide each matrix row/column from a hdf5 dataset 
//' using a vector
//' 
//' @param filename string file name where dataset to apply weights is located
//' @param group string with the path inside the hdf5 data file where matrix 
//' is located
//' @param dataset string with the matrix name
//' @param vectorgroup string with the path inside the hdf5 data file where 
//' vector is located
//' @param vectordataset string with the vector name
//' @param outdataset character array with output dataset name where we want to 
//' store results
//' @param outgroup optional, character array with output group name where we 
//' want to store results if not provided then results are stored in the same 
//' group as original dataset
//' @param func, Character array, function to be applyed : 
//'"+" : to sum a vector to a matrix dataset by columns or rows
//'"-" : to substract a vector to a matrix dataset by columns or rows
//'"*" : to multiply a vector to a matrix dataset by columns or rows
//'"/" : to divide a vector to a matrix dataset by columns or rows
//' @param byrows logical (default = FALSE). By default weights are applied by 
//' columns but if byrows=TRUE then weights are applied by rows 
//' @param force, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
//' @return file with weighted dataset
//' @examples
//'library(BigDataStatMeth)
//'    
//'# Prepare data and functions
//'set.seed(123)
//'Y <- matrix(rnorm(250), 10, 10)
//'X <- matrix(rnorm(250), 10, 1)
//'        
//'# Create hdf5 data file with  data (Y)
//'bdCreate_hdf5_matrix_file("test.hdf5", Y, "data", "Y", force = T)
//'bdAdd_hdf5_matrix( X, "test.hdf5",  "data", "X", force = TRUE)
//'            
//'bdcomputeMatrixVector_hdf5("test.hdf5", 
//'                           group = "data", dataset = "Y",
//'                           vectorgroup = "data", vectordataset = "X", 
//'                           outdataset = "ProdComputed", 
//'                           func = "*",
//'                           byrows = T, force = T)
//'    
//'bdcomputeMatrixVector_hdf5("test.hdf5", 
//'                           group = "data", dataset = "Y",
//'                           vectorgroup = "data", vectordataset = "X", 
//'                           outdataset = "SumComputed", 
//'                           func = "-",
//'                           byrows = T, force = T)
//'    
//'bdcomputeMatrixVector_hdf5("test.hdf5", 
//'                           group = "data", dataset = "Y",
//'                           vectorgroup = "data", vectordataset = "X", 
//'                           outdataset = "SubsComputed", 
//'                           func = "-",
//'                           byrows = F, force = T)
//' @export
// [[Rcpp::export]]
void bdcomputeMatrixVector_hdf5( std::string filename, std::string group, 
                            std::string dataset,
                            std::string vectorgroup, std::string vectordataset,
                            std::string outdataset, 
                            std::string func,
                            Rcpp::Nullable<std::string> outgroup = R_NilValue,
                            Rcpp::Nullable<bool> byrows = R_NilValue,
                            Rcpp::Nullable<int> force  = false)
{
    
    bool bforce, bbyrows;
    int blocksize;
    std::string strgroupout;
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1);
    
    H5File* file = nullptr;
    DataSet* pdatasetin = nullptr;
    DataSet* pvectorin = nullptr;
    DataSet* pdatasetout = nullptr;
    
    Rcpp::NumericVector oper = {0, 1, 2, 3};
    oper.names() = Rcpp::CharacterVector({"+", "-", "*", "/"});
    
    try{
        
        if( byrows.isNull()) {
            bbyrows = false;
        } else {
            bbyrows = Rcpp::as<bool> (byrows);
        }
        
        if( force.isNull()) {
            bforce = true;
        } else {
            bforce = Rcpp::as<bool> (force);
        }
        
        if( outgroup.isNull()) {
            strgroupout = group;
        } else {
            strgroupout = Rcpp::as<std::string> (outgroup);
        }
        
        // Function exists?
        if( oper(oper.findName(func)) != 0 && oper(oper.findName(func)) != 1 &&
            oper(oper.findName(func)) != 2 && oper(oper.findName(func)) != 3) 
        {
            Rcpp::Rcout<<"Function does not exists, please use one of the following : '+', '-', '*', '/' ";
            return void();
        } 
        
        if(!ResFileExist(filename)) {
            Rcpp::Rcout<<"\nFile not exits, create file before normalize dataset";  
            return void();
        }
        
        file = new H5File( filename, H5F_ACC_RDWR );
        
        if(exists_HDF5_element_ptr(file, group)==0) {
            Rcpp::Rcout<<"\nGroup not exits, create file and matrix dataset before proceed";
            file->close();
            return void();
        } else {
            if(!exists_HDF5_element_ptr(file, group + "/" + dataset)) {
                Rcpp::Rcout<<"\nDataset not exits, create matrix dataset before proceed";
                file->close();
                return void();
            }
        }
        
        if(exists_HDF5_element_ptr(file, vectorgroup)==0) {
            Rcpp::Rcout<<"\nGroup not exits, create file and vector dataset before proceed";
            file->close();
            return void();
        } else {
            if(!exists_HDF5_element_ptr(file, vectorgroup + "/" + vectordataset)) {
                Rcpp::Rcout<<"\nDataset not exits, create vector dataset before proceed";
                file->close();
                return void();
            }
        }
        
        if(exists_HDF5_element_ptr(file, strgroupout + "/" + outdataset) && bforce == false) {
            Rcpp::Rcout<<"\n Output dataset exists, please set force = TRUE to overwrite\n";
            file->close();
            return void();
        } else if( exists_HDF5_element_ptr(file, strgroupout) && bforce == true) {
            remove_HDF5_element_ptr(file, strgroupout + "/" + outdataset); 
        }
        
        pdatasetin = new DataSet(file->openDataSet(group + "/" + dataset));
        IntegerVector dims_out = get_HDF5_dataset_size_ptr(pdatasetin);
        
        pvectorin = new DataSet(file->openDataSet(vectorgroup + "/" + vectordataset));
        IntegerVector vdims_out = get_HDF5_dataset_size_ptr (pvectorin);
        
        if( (bbyrows == false && vdims_out[1] != dims_out[1]) || vdims_out[0]!=1 ) {
            Rcpp::Rcout<<"\nNon-conformable dimensions or vector dataset is not a vector";
            pvectorin->close();
            pdatasetin->close();
            file->close();
            return void();
        } else if( bbyrows == true && vdims_out[1] != dims_out[0]) {
            Rcpp::Rcout<<"\nNon-conformable dimensions - byrows = true";
            pvectorin->close();
            pdatasetin->close();
            file->close();
            return void();
        }
        
        Eigen::MatrixXd vWeights_tmp =  GetCurrentBlock_hdf5( file, pvectorin, 0, 0, vdims_out[0], vdims_out[1]);
        Eigen::VectorXd vWeights = vWeights_tmp.row(0);
        pvectorin->close();
        
        if( bbyrows == false) {
            // Define blocksize atending number of elements in rows and cols
            if(dims_out[1] > maxElemBlock){
                blocksize = 1;
            } else {
                int maxsize = std::max( dims_out[0], dims_out[1]);
                blocksize = std::ceil( maxElemBlock / maxsize);
            }
            
        } else {
            // Define blocksize atending number of elements in rows and cols
            if(dims_out[0] > maxElemBlock) {
                blocksize = 1;
            } else {
                int maxsize = std::max( dims_out[0], dims_out[1]);
                blocksize = std::ceil( maxElemBlock / maxsize);
            }
        }
        
        // if not exists -> create output group 
        if(exists_HDF5_element_ptr(file, strgroupout) == 0) {
            create_HDF5_groups_ptr( file, strgroupout);
        }
        
        if( bbyrows == false) {
            for(hsize_t i=0; i*blocksize <= dims_out[0] ; i++)
            {
                int sizetoread = 0;
                if((i+1)*blocksize<dims_out[0]) {
                    sizetoread = blocksize;
                } else {
                    sizetoread = dims_out[0]-(i*blocksize);
                }
                
                // Prepare file and dataset
                if(i==0) {
                    create_HDF5_dataset_ptr(file, strgroupout + "/" + outdataset, dims_out[0], dims_out[1], "real");
                    pdatasetout = new DataSet(file->openDataSet(strgroupout + "/" + outdataset));
                }
                
                IntegerVector offset = IntegerVector::create( i*blocksize, 0);
                IntegerVector count = IntegerVector::create(sizetoread, dims_out[1] );
                // Compute and write data
                Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, pdatasetin, offset[0], offset[1], count[0], count[1]);
                
                if( oper(oper.findName( func )) == 0 ) {
                    X = Rcpp_matrixVectorSum_byCol(X, vWeights); 
                } else if( oper(oper.findName( func )) == 1 ) {
                    X = Rcpp_matrixVectorSubstract_byCol(X, vWeights);
                }else if( oper(oper.findName( func )) == 2 ) {
                    X = Rcpp_matrixVectorMultiplication_byCol(X, vWeights);
                }else if( oper(oper.findName( func )) == 3 ) {
                    X = Rcpp_matrixVectorDivision_byCol(X, vWeights);
                }
                // X = X.array().rowwise() * vWeights.transpose().array();
                write_HDF5_matrix_subset_v2(file, pdatasetout, offset, count, block, stride, wrap(X));
            }
        } else {
            
            for(hsize_t i=0; i*blocksize <= dims_out[1] ; i++) {
                int sizetoread = 0;
                if( (i+1)*blocksize < dims_out[1] ){
                    sizetoread = blocksize;
                } else {
                    sizetoread = dims_out[1]-(i*blocksize);
                }
                
                // Prepare file and dataset
                if(i==0) {
                    create_HDF5_dataset_ptr(file, strgroupout + "/" + outdataset, dims_out[0], dims_out[1], "real");
                    pdatasetout = new DataSet(file->openDataSet(strgroupout + "/" + outdataset));
                }
                
                IntegerVector offset = IntegerVector::create( 0, i*blocksize);
                IntegerVector count = IntegerVector::create(dims_out[0], sizetoread );
                
                // Compute and write data
                Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, pdatasetin, offset[0], offset[1], count[0], count[1]);
                
                if( oper(oper.findName( func )) == 0 ) {
                    X = Rcpp_matrixVectorSum_byRow(X, vWeights); 
                } else if( oper(oper.findName( func )) == 1 ) {
                    X = Rcpp_matrixVectorSubstract_byRow(X, vWeights);
                }else if( oper(oper.findName( func )) == 2 ) {
                    X = Rcpp_matrixVectorMultiplication_byRow(X, vWeights);
                }else if( oper(oper.findName( func )) == 3 ) {
                    X = Rcpp_matrixVectorDivision_byRow(X, vWeights);
                }
                //..// X = X.array().colwise() * vWeights.array();
                write_HDF5_matrix_subset_v2(file, pdatasetout, offset, count, block, stride, wrap(X));
            }
        }
        
    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        ::Rf_error( "c++ exception bdWeightedProduct_hdf5 (File IException)" );
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        ::Rf_error( "c++ exception bdWeightedProduct_hdf5 (DataSet IException)" );
        return void();
    } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        ::Rf_error( "c++ exception bdWeightedProduct_hdf5 (DataSpace IException)" );
        return void();
    } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        ::Rf_error( "c++ exception bdWeightedProduct_hdf5 (DataType IException)" );
        return void();
    }catch(std::exception &ex) {
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        Rcpp::Rcout<< ex.what();
        return void();
    }
    
    pdatasetin->close();
    pdatasetout->close();
    file->close();
    
    Rcpp::Rcout<<"Calculus has been computed\n";
    return void();
}


/****R

library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))
# setwd("/Users/mailos/DOCTORAT_Local/BigDataStatMeth/")

# Prepare data and functions
set.seed(123)
Y <- matrix(rnorm(250), 10, 10)
X <- matrix(rnorm(250), 10, 1)

# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("test.hdf5", Y, "data", "Y", force = T)
bdAdd_hdf5_matrix( X, "test.hdf5",  "data", "X", force = TRUE)


bdcomputeMatrixVector_hdf5("test.hdf5", 
                       group = "data", dataset = "Y",
                       vectorgroup = "data", vectordataset = "X", 
                       outdataset = "ProdComputed", 
                       func = "*",
                       byrows = T, force = T)

bdcomputeMatrixVector_hdf5("test.hdf5", 
                           group = "data", dataset = "Y",
                           vectorgroup = "data", vectordataset = "X", 
                           outdataset = "SumComputed", 
                           func = "-",
                           byrows = T, force = T)

bdcomputeMatrixVector_hdf5("test.hdf5", 
                           group = "data", dataset = "Y",
                           vectorgroup = "data", vectordataset = "X", 
                           outdataset = "SubsComputed", 
                           func = "-",
                           byrows = F, force = T)


setwd("/Users/mailos/DOCTORAT_Local/mgcca")
filename <- "tmp/gettables2.hdf5"
bdcomputeMatrixVector_hdf5(filename, group = "FinalRes", dataset = "M",
                           vectorgroup = "FinalRes", vectordataset = "Ksum05",
                           outdataset = "MMKsum05", func = "*",
                           byrows = T,force = T)






*/