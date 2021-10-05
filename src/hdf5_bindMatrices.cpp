#include "include/hdf5_bindMatrices.h"


// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;


//' Bind matrices by rows or columns
//'
//' Merge existing matrices inside hdf5 data file by rows or by columns
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param datasets, character array indicating the input dataset to be imputed
//' @param func, character array function to be applyed
//' \describe{
//'     \item{bindRows}{merge datasets by rows}
//'     \item{bindCols}{apply datasets by columns}
//' }
//' @param outgroup, character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, character array indicating the name for the new merged dataset
//' @param bforce, boolean if true, previous results in same location inside hdf5 will be overwritten.
//' @return Original hdf5 data file with results after input datasets
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdBind_hdf5( std::string filename, std::string group, Rcpp::StringVector datasets, 
                           std::string outgroup, std::string outdataset, std::string func,
                           Rcpp::Nullable<bool> force = false )
{
    
    H5File* file;
    DataSet* pdataset = nullptr;
    DataSet* unlimDataset = nullptr;
    Rcpp::NumericVector oper = {0, 1};
    oper.names() = Rcpp::CharacterVector({ "bindCols", "bindRows"});
    
    try
    {
        
        IntegerVector count = IntegerVector::create(0, 0);
        IntegerVector offset = IntegerVector::create(0, 0);
        IntegerVector stride = IntegerVector::create(1, 1);
        IntegerVector block = IntegerVector::create(1, 1);
        
        bool bforce;
        
        if(force.isNull()) { bforce = false; } 
        else {   bforce = Rcpp::as<bool>(force); }
        
        
        // Test file
        if( ResFileExist_filestream(filename) ) {
            file = new H5File( filename, H5F_ACC_RDWR ); 
        } else {
            Rcpp::Rcout<<"\nFile not exits, create file before bind matrices";
            return wrap(false);
        }
        
        
        // Seek all datasets to perform calculus
        for( int i=0; i < datasets.size(); i++ ) 
        {
            std::string strdataset = group +"/" + datasets(i);
            std::string stroutDatasetName = outgroup + "/" + outdataset;
            
            if( exists_HDF5_element_ptr(file, strdataset ) == 0 ) {
                file->close();
                Rcpp::Rcout<<"Group or dataset does not exists, please create the input dataset before proceed";
                return wrap(false);
            }
            
            pdataset = new DataSet(file->openDataSet(strdataset));
            
            // Real data set dimension
            IntegerVector dims_out = get_HDF5_dataset_size(*pdataset);
            
            // Get block from complete matrix
            Eigen::MatrixXd original = GetCurrentBlock_hdf5( file, pdataset, 0, 0, dims_out[0], dims_out[1]);
            
            // Remove dataset if exists ( only if force = TRUE )
            if(i==0) {
                prepare_outGroup(file, outgroup, bforce);
                prepare_outDataset(file, outgroup + "/" + outdataset, bforce);
            }
            
            if( oper.findName( func ) == 0 || oper.findName( func ) == 1) {
                
                // Rcpp::Rcout<< "Files actuals : "<<original.rows()<<"\n";
                // Rcpp::Rcout<< "Files llegides : "<<count[0]<<"\n";
                // Rcpp::Rcout<< "Columnes actuals : "<<original.cols()<<"\n";
                // Rcpp::Rcout<< "Columnes llegides : "<<count[1]<<"\n";
                
                
                if(oper.findName( func ) == 0 ){
                    
                    // Test if dimmensions are correct
                    if( original.cols() != count[1] && i!=0) {
                        pdataset->close();
                        file->close();
                        ::Rf_error( "c++ exception can't bind datasets by columns, number of columns differ between datasets" );
                        return (wrap(false));
                    }
                    offset[0] = offset[0] + count[0];
                } else {
                    
                    // Test if dimmensions are correct
                    if( original.rows() != count[0]  && i!=0) {
                        pdataset->close();
                        file->close();
                        ::Rf_error( "c++ exception can't bind datasets by rows, number of rows differ between datasets" );
                        return (wrap(false));
                    }
                    offset[1] = offset[1] + count[1];
                }
                
                count[0] = original.rows();
                count[1] = original.cols();
                
                if(i == 0) {
                    // If dataset exists --> remove dataset
                    if( exists_HDF5_element_ptr(file,stroutDatasetName))
                        remove_HDF5_element_ptr(file,stroutDatasetName);
                    // Create unlimited dataset in hdf5 file
                    create_HDF5_unlimited_matrix_dataset_ptr(file, stroutDatasetName, count[0], count[1], "numeric");
                }
                    
                // Rcpp::Rcout<<"\n Where peta - 1 \n";
                unlimDataset = new DataSet(file->openDataSet(stroutDatasetName));
                
                Rcpp::Rcout<<"\n Where peta - 2 \n";
                
                if(oper.findName( func ) == 0 && i!=0) {
                    
                    // Rcpp::Rcout<<"\n Where peta - 3.1 \n";
                    extend_HDF5_matrix_subset_ptr(file, unlimDataset, count[0], 0);
                    // Rcpp::Rcout<<"\n Where peta - 3 \n";
                    
                } else if (oper.findName( func ) == 1 && i!=0) {
                    
                    // Rcpp::Rcout<<"\n Where peta - 4.1 \n";
                    extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
                    // Rcpp::Rcout<<"\n Where peta - 4 \n";
                    
                }
                
                // Rcpp::Rcout<<"\n Where peta - 5 \n";
                // Rcpp::Rcout<< "Offset [0] : "<<offset[0]<<"\n";
                // Rcpp::Rcout<< "Offset [1] : "<<offset[1]<<"\n";
                // Rcpp::Rcout<< "Count [0] : "<<count[0]<<"\n";
                // Rcpp::Rcout<< "Count [1] : "<<count[1]<<"\n";
                
                
                write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(original)  );  
                // Rcpp::Rcout<<"\n Where peta - 6 \n";
                unlimDataset->close();
                
                
                pdataset->close();

            }else {
                pdataset->close();
                file->close();
                Rcpp::Rcout<<"Group not exists, create the input dataset before proceed";
                return wrap(false);
                
            }
            
        }
        
        
    }
    catch( FileIException& error ) { // catch failure caused by the H5File operations
        unlimDataset->close();
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




data("mtcars")
library(BigDataStatMeth)
library(rhdf5)

#..# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/UAB/DOCTORAT/BitDataStatMeth - BDSM/Analysis/BigDataStatMeth_Analysis/Cholesterol/test")

bdSplit_matrix_hdf5( "cars.hdf5", "data", "X", "dataoutCols", nblocks = 3, bycols = FALSE, force = TRUE)


devtools::reload(pkgload::inst("BigDataStatMeth"))

x.blocks <- BigDataStatMeth::bdinGroupDatasetin_hdf5("cars.hdf5", "Xrows")
y.blocks <- BigDataStatMeth::bdinGroupDatasetin_hdf5("cars.hdf5", "Yrows")

bdBind_hdf5("cars.hdf5", "Xrows", x.blocks, "merges", "outMergedRows", "bindRows", force = TRUE )
bdBind_hdf5("cars.hdf5", "Yrows", y.blocks, "merges", "outMergedCols", "bindCols", force = TRUE )





bdapply_Function_hdf5( "cars.hdf5", "Xrows", x.blocks, "Xrows_CrossProd", "CrossProd", force = TRUE )
bdapply_Function_hdf5( "cars.hdf5", "Xrows", x.blocks, "Xrows_tCrossProd", "tCrossProd", force = TRUE )



# Test results hdf5 : 
# Examine hierarchy before open file
h5ls("cars.hdf5")

# Open file and get data, all data is stored under SVD group
h5f = H5Fopen("cars.hdf5")

mergedRows <- h5f$Xcols_merge$outMergedCols

h5closeAll()

mergedRows
*/