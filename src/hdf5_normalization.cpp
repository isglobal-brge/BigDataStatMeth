#include "include/hdf5_normalization.h"

// Internal call 
Eigen::MatrixXd RcppNormalize_Data_hdf5 ( Eigen::MatrixXd  X, bool bc, bool bs, 
                                          bool btransp, Eigen::MatrixXd normdata )
{
    Eigen::MatrixXd rX;
    
    if( btransp == true) {
        if( bc==true && bs==true )  {
            rX = (X.colwise() - normdata.row(0).transpose() ).array().colwise() / normdata.row(1).transpose().array();
        }   else if (bc == true  && bs==false)   {
            rX = (X.colwise() - normdata.row(0).transpose());
        }  else if ( bc == false && bs == true)   {
            rX = X.array().colwise() / normdata.row(1).transpose().array();
        }
    } else {
        if( bc==true && bs==true )  {
            rX = (X.rowwise() - normdata.row(0)).array().rowwise() / normdata.row(1).array();
        } else if (bc == true  && bs==false) {
            rX = (X.rowwise() - normdata.row(0));
        }  else if ( bc == false && bs == true)   {
            rX = X.array().rowwise() / normdata.row(1).array();
        } 
    }
    
    return(rX);
}



//' Normalize dataset in hdf5 file
//' 
//' This function normalize data scaling, centering or scaling and centering 
//' in a dataset stored in hdf5 file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string Matrix
//' @param dataset string Matrix
//' @param bcenter logical (default = TRUE) if TRUE, centering is done by 
//' subtracting the column means
//' @param bscale logical (default = TRUE) if TRUE, centering is done by 
//' subtracting the column means
//' @param byrows logical (default = FALSE) if TRUE, centering is done by 
//' subtracting the rows means, util when working with hdf5 datasets stored 
//' in Row Major format.
//' @param wsize integer (default = 1000), file block size to read to 
//' perform normalization
//' @param force, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
//' @return file with scaled, centered or scaled and centered dataset
//' @examples
//'   a = "See vignette"
//' @export
// [[Rcpp::export]]
void bdNormalize_hdf5( std::string filename, const std::string group, std::string dataset,
                                Rcpp::Nullable<bool> bcenter = R_NilValue, Rcpp::Nullable<bool> bscale  = R_NilValue,
                                Rcpp::Nullable<bool> byrows = R_NilValue,
                                Rcpp::Nullable<int> wsize  = R_NilValue, Rcpp::Nullable<int> force  = false)
{
  
    bool bc, bs, bforce, bbyrows;
    int blocksize;
    std::string strgroupout;
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1);
    
    Eigen::MatrixXd datanormal;
    
    H5File* file = nullptr;
    DataSet* pdatasetin = nullptr;
    DataSet* pdatasetout = nullptr;

    try{
      
        if( bcenter.isNull())
            bc = true;
        else
            bc = Rcpp::as<bool> (bcenter);
        
        if( bscale.isNull())
            bs = true;
        else
            bs = Rcpp::as<bool> (bscale);
        
        if( byrows.isNull())
            bbyrows = false;
        else
            bbyrows = Rcpp::as<bool> (byrows);
        
        if( force.isNull())
            bforce = true;
        else
            bforce = Rcpp::as<bool> (force);
        
        
        if(!ResFileExist(filename)) {
            Rcpp::Rcout<<"\nFile not exits, create file before normalize dataset\n";  
            return void();
        }
        
        file = new H5File( filename, H5F_ACC_RDWR );
        
        if(exists_HDF5_element_ptr(file, group)==0) {
            Rcpp::Rcout<<"\nGroup not exits, create file and dataset before normalize data\n";
            file->close();
            return void();
        } else {
            if(!exists_HDF5_element_ptr(file, group + "/" + dataset)) {
                Rcpp::Rcout<<"\n Dataset not exits, create dataset before normalize data \n";
                file->close();
                return void();
            }
        }
        
        strgroupout = "NORMALIZED/" + group;
        
        if(exists_HDF5_element_ptr(file, strgroupout + "/" + dataset) && bforce == false) {
            Rcpp::Rcout<<"\n Normalized dataset exists, please set force = TRUE to overwrite\n";
            file->close();
            return void();
        } else if( exists_HDF5_element_ptr(file, strgroupout) && bforce == true) {
            remove_HDF5_element_ptr(file, strgroupout + "/" + dataset); 
            remove_HDF5_element_ptr(file, strgroupout + "/mean." + dataset); 
            remove_HDF5_element_ptr(file, strgroupout + "/sd." + dataset); 
        }
        
        pdatasetin = new DataSet(file->openDataSet(group + "/" + dataset));
        IntegerVector dims_out = get_HDF5_dataset_size(*pdatasetin);
        
        if( bbyrows == false) {
            // Define blocksize atending number of elements in rows and cols
            if( wsize.isNull()) {
                if(dims_out[1] > maxElemBlock){
                    blocksize = 1;
                } else {
                    int maxsize = std::max( dims_out[0], dims_out[1]);
                    blocksize = std::ceil( maxElemBlock / maxsize);
                }
            } else {
                if(dims_out[1] > maxElemBlock){
                    blocksize = 1;
                } else {
                    blocksize = Rcpp::as<int> (wsize);
                }
            }
            
            datanormal = Eigen::MatrixXd::Zero(2,dims_out[0]);
            // Get data to normalize matrix (mean and sd by column)
            get_HDF5_mean_sd_by_column_ptr( file, pdatasetin, datanormal);
        } else {
            // Define blocksize atending number of elements in rows and cols
            if( wsize.isNull()) {
                if(dims_out[0] > maxElemBlock) {
                    blocksize = 1;
                } else {
                    int maxsize = std::max( dims_out[0], dims_out[1]);
                    blocksize = std::ceil( maxElemBlock / maxsize);
                }
            } else {
                if(dims_out[0] > maxElemBlock) {
                    blocksize = 1;
                } else {
                    blocksize = Rcpp::as<int> (wsize);
                }
            }
            
            datanormal = Eigen::MatrixXd::Zero(2,dims_out[1]);
            // Get data to normalize matrix (mean and sd by column)
            get_HDF5_mean_sd_by_row_ptr( file, pdatasetin, datanormal);
        }
        
        // if not exists -> create output group 
        if(exists_HDF5_element_ptr(file, strgroupout) == 0) {
            create_HDF5_groups_ptr( file, strgroupout);
        }
        
        // Store center and scale for each column
        write_HDF5_matrix_ptr(file, strgroupout + "/mean." + dataset, wrap(datanormal.row(0)));
        write_HDF5_matrix_ptr(file, strgroupout + "/sd." + dataset, wrap(datanormal.row(1)));
        
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
                    create_HDF5_dataset_ptr(file, strgroupout + "/" + dataset, dims_out[0], dims_out[1], "real");
                    pdatasetout = new DataSet(file->openDataSet(strgroupout + "/" + dataset));
                }
                
                IntegerVector offset = IntegerVector::create( i*blocksize, 0);
                IntegerVector count = IntegerVector::create(sizetoread, dims_out[1] );
                
                // Normalize and write data
                Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, pdatasetin, offset[0], offset[1], count[0], count[1]);
                X = RcppNormalize_Data_hdf5(X, bc, bs, true, datanormal.block(0,offset[0], 2, count[0]));
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
                    // Create dataset
                    create_HDF5_dataset_ptr(file, strgroupout + "/" + dataset, dims_out[0], dims_out[1], "real");
                    pdatasetout = new DataSet(file->openDataSet(strgroupout + "/" + dataset));
                }
                
                IntegerVector offset = IntegerVector::create( 0, i*blocksize);
                IntegerVector count = IntegerVector::create(dims_out[0], sizetoread );
                
                // Normalize and write data
                Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, pdatasetin, offset[0], offset[1], count[0], count[1]);
                X = RcppNormalize_Data_hdf5(X, bc, bs, false, datanormal.block(0,offset[1], 2, count[1]));
                write_HDF5_matrix_subset_v2(file, pdatasetout, offset, count, block, stride, wrap(X));
            }
        }
    
    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (File IException)" );
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataSet IException)" );
        return void();
    } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataSpace IException)" );
        return void();
    } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        pdatasetin->close();
        pdatasetout->close();
        file->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataType IException)" );
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
    
    Rcpp::Rcout<<"Normalization has been computed\n";
    return void();
}


/****R

    library(BigDataStatMeth)
    # devtools::reload(pkgload::inst("BigDataStatMeth"))
    
    # setwd("/Users/mailos/DOCTORAT_Local/BigDataStatMeth/")
    
    # Prepare data and functions
    Y <- matrix(rnorm(250), 10, 10)
    X <- matrix(rnorm(250), 10, 10)
    
    # Create hdf5 data file with  data (Y)
    bdCreate_hdf5_matrix_file("cca_cars.hdf5", Y, "data", "Y", force = T)
    bdAdd_hdf5_matrix( t(Y), "cca_cars.hdf5",  "data", "Yt", force = TRUE)
    
    bdNormalize_hdf5(filename = "cca_cars.hdf5", 
                     group = "data", dataset = "Y", 
                     bcenter = TRUE, bscale = TRUE)
    
    bdNormalize_hdf5(filename = "cca_cars.hdf5", 
                     group = "data", dataset = "Yt", byrows = T,
                     bcenter = TRUE, bscale = TRUE)
    
*/