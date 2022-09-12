#include "include/hdf5_getSDandMean.h"
using namespace Rcpp;



// Get mean and corrected sd from each column in dataset in the case of n<<m, this information is used
// to normalize data, center or scale.
int get_HDF5_mean_sd_by_column_ptr(H5File* file, DataSet* dataset, Eigen::MatrixXd& normalize )
{
    
    IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
    
    try
    {
        
        // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
        Exception::dontPrint();
        
        int block_size = 500;
        
        IntegerVector stride = IntegerVector::create(1, 1);
        IntegerVector block = IntegerVector::create(1, 1);
        IntegerVector offset = IntegerVector::create(0, 0);
        IntegerVector count = IntegerVector::create(0, 0);
        
        count[1] = dims_out[1];
        if( block_size < dims_out[0] )
            count[0] = block_size;
        else
            count[0] = dims_out[0];
        
        // Read data in blocks of 500 columns
        for(int i=0; (i < floor(dims_out[0]/block_size)) || i==0; i++)
        {
            
            if(i>0){
                offset[0] = offset[0] + block_size;
                
                if( offset[0] + block_size <= dims_out[0] ) {
                    count[0] = block_size;
                }else {
                    count[0] = dims_out[0] - offset[1]+block_size; 
                }
            }
            
            Eigen::MatrixXd X = GetCurrentBlock_hdf5(file, dataset, offset[0], offset[1], count[0], count[1]);
            
            Eigen::VectorXd mean = X.rowwise().mean();
            Eigen::VectorXd std = ((X.colwise() - mean).array().square().rowwise().sum() / (X.cols() - 1)).sqrt();
            
            normalize.block( 0, offset[0] , 1, mean.size()) = mean.transpose();
            normalize.block( 1, offset[0], 1, std.size()) = std.transpose();
            
        }
        
    } catch(FileIException& error) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (File IException)" );
        return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (DataSet IException)" );
        return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (Group IException)" );
        return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (DataSpace IException)" );
        return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_column_ptr (Data TypeIException)" );
        return -1;
    }
    
    return(0);  // successfully terminated
    
}


// Get mean and corrected sd from each row in dataset 
int get_HDF5_mean_sd_by_row_ptr(H5File* file, DataSet* dataset, Eigen::MatrixXd& normalize )
{
    
    IntegerVector dims_out = get_HDF5_dataset_size(*dataset);
    
    try
    {
        
        // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
        Exception::dontPrint();
        
        int block_size = 500;
        
        IntegerVector stride = IntegerVector::create(1, 1);
        IntegerVector block = IntegerVector::create(1, 1);
        IntegerVector offset = IntegerVector::create(0, 0);
        IntegerVector count = IntegerVector::create(0, 0);
        
        count[0] = dims_out[0];
        if( block_size < dims_out[1] ) {
            count[1] = block_size;
        } else{
            count[1] = dims_out[1];
        }
        
        // Read data in blocks of 500 columns
        for(int i=0; (i < floor(dims_out[1]/block_size)) || i==0; i++)
        {
            
            if(i>0){
                offset[1] = offset[1] + block_size;
                
                if( offset[1] + block_size <= dims_out[1] ) {
                    count[1] = block_size;
                }else {
                    count[1] = dims_out[1] - offset[0] + block_size;
                }
            }
            
            Eigen::MatrixXd X = GetCurrentBlock_hdf5(file, dataset, offset[0], offset[1], count[0], count[1]);
            
            Eigen::RowVectorXd mean = X.colwise().mean();
            Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
            
            normalize.block( 0, offset[1], 1, mean.size()) = mean;
            normalize.block( 1, offset[1], 1, std.size()) = std;
            
        }
        
    } catch(FileIException& error) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_row_ptr (File IException)" );
        return -1;
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_row_ptr (DataSet IException)" );
        return -1;
    } catch(GroupIException& error) { // catch failure caused by the Group operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_row_ptr (Group IException)" );
        return -1;
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_row_ptr (DataSpace IException)" );
        return -1;
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception get_HDF5_mean_sd_by_row_ptr (Data TypeIException)" );
        return -1;
    }
    
    return(0);  // successfully terminated
    
}

//' Get sd and Mean by Rows or Columns
//' 
//' This functions gets Standard Deviation (sd) or Mean by Rows or Columns and
//' store results in hdf5 dataset inside the file
//' 
//' @param filename string file name where dataset to normalize is stored
//' @param group string Matrix
//' @param dataset string Matrix
//' @param sd logical (default = TRUE) if TRUE, standard deviation is computed
//' @param mean logical (default = TRUE) if TRUE, mean is computed 
//' @param byrows logical (default = FALSE) if TRUE, sd and mean are computed
//' by columns, if byrows=TRUE then sd and mean are computed by Rows.
//' @param wsize integer (default = 1000), file block size to read to 
//' perform calculus exitexit
//' @param force, boolean if true, previous results in same location inside 
//' hdf5 will be overwritten.
//' @return hdf5 data file containing a dataset with sd, and mean 
//' @examples
//' 
//' library(BigDataStatMeth)
//' # devtools::reload(pkgload::inst("BigDataStatMeth"))
//' setwd("/Users/mailos/DOCTORAT_Local/BigDataStatMeth/")
//'     
//' # Prepare data and functions
//' et.seed(123)
//' Y <- matrix(rnorm(250), 10, 10)
//' X <- matrix(rnorm(250), 10, 1)
//'     
//' # Create hdf5 data file with  data (Y)
//' bdCreate_hdf5_matrix_file("test.hdf5", Y, "data", "Y", force = T)
//' bdAdd_hdf5_matrix( X, "test.hdf5",  "data", "X", force = TRUE)
//' 
//' # Get mean and sd        
//' bdgetSDandMean_hdf5(filename = "test.hdf5", group = "data", dataset = "Y",
//'                     sd = T, mean = T,byrows = T)
//'         
//' @export
// [[Rcpp::export]]
void bdgetSDandMean_hdf5( std::string filename, const std::string group, 
                    std::string dataset,
                    Rcpp::Nullable<bool> sd = R_NilValue, 
                    Rcpp::Nullable<bool> mean  = R_NilValue,
                    Rcpp::Nullable<bool> byrows = R_NilValue,
                    Rcpp::Nullable<int> wsize  = R_NilValue, 
                    Rcpp::Nullable<int> force  = false)
{
    
    
    bool bsd, bmean, bforce, bbyrows;
    int blocksize;
    std::string strgroupout;
    IntegerVector stride = IntegerVector::create(1, 1);
    IntegerVector block = IntegerVector::create(1, 1);
    
    Eigen::MatrixXd datanormal;
    
    H5File* file = nullptr;
    DataSet* pdatasetin = nullptr;
    
    try {
        
        
        if( mean.isNull()) {
            bmean = true;
        } else {
            bmean = Rcpp::as<bool> (mean); }
        
        if( sd.isNull()) {
            bsd = true;
        } else {
            bsd = Rcpp::as<bool> (sd); }
        
        if( byrows.isNull()) {
            bbyrows = false;
        } else {
            bbyrows = Rcpp::as<bool> (byrows); }
        
        if( force.isNull()) {
            bforce = true;
        } else {
            bforce = Rcpp::as<bool> (force); }
        
        
        if(!ResFileExist(filename)) {
            Rcpp::Rcout<<"\nFile not exits, create file before get mean and sd\n";  
            return void();
        }
        
        file = new H5File( filename, H5F_ACC_RDWR );
        
        if(exists_HDF5_element_ptr(file, group)==0) {
            Rcpp::Rcout<<"\nGroup not exits, create file and dataset before get mean and sd\n";
            file->close();
            return void();
        } else {
            if(!exists_HDF5_element_ptr(file, group + "/" + dataset)) {
                Rcpp::Rcout<<"\n Dataset not exits, create dataset before get mean and sd \n";
                file->close();
                return void();
            }
        }
        
        strgroupout = "mean_sd/" + group;
        
        if(exists_HDF5_element_ptr(file, strgroupout + "/" + dataset) && bforce == false) {
            Rcpp::Rcout<<"\n Mean and sd exists please set force = TRUE to overwrite\n";
            file->close();
            return void();
        } else if( exists_HDF5_element_ptr(file, strgroupout) && bforce == true) {
            remove_HDF5_element_ptr(file, strgroupout + "/" + dataset); 
            remove_HDF5_element_ptr(file, strgroupout + "/" + dataset + ".mean"); 
            remove_HDF5_element_ptr(file, strgroupout + "/" + dataset + ".sd"); 
        }
        
        pdatasetin = new DataSet(file->openDataSet(group + "/" + dataset));
        IntegerVector dims_out = get_HDF5_dataset_size_ptr(pdatasetin);
        
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
        write_HDF5_matrix_ptr(file, strgroupout + "/"+dataset + ".mean", wrap(datanormal.row(0)));
        write_HDF5_matrix_ptr(file, strgroupout + "/"+dataset + ".sd", wrap(datanormal.row(1)));
        
    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdatasetin->close();
        file->close();
        ::Rf_error( "c++ exception bdgetSDandMean_hdf5 (File IException)" );
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdatasetin->close();
        file->close();
        ::Rf_error( "c++ exception bdgetSDandMean_hdf5 (DataSet IException)" );
        return void();
    } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        pdatasetin->close();
        file->close();
        ::Rf_error( "c++ exception bdgetSDandMean_hdf5 (DataSpace IException)" );
        return void();
    } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        pdatasetin->close();
        file->close();
        ::Rf_error( "c++ exception bdgetSDandMean_hdf5 (DataType IException)" );
        return void();
    }catch(std::exception &ex) {
        pdatasetin->close();
        file->close();
        Rcpp::Rcout<< ex.what();
        return void();
    }
    
    pdatasetin->close();
    file->close();
    
}



/****R

library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))
setwd("/Users/mailos/DOCTORAT_Local/BigDataStatMeth/")

# Prepare data and functions
set.seed(123)
Y <- matrix(rnorm(250), 10, 10)
X <- matrix(rnorm(250), 10, 1)

# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("test.hdf5", Y, "data", "Y", force = T)
bdAdd_hdf5_matrix( X, "test.hdf5",  "data", "X", force = TRUE)

bdgetSDandMean_hdf5(filename = "test.hdf5", group = "data", dataset = "Y",
                    sd = T, mean = T,byrows = T)


*/