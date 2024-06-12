#include <BigDataStatMeth.hpp>
#include "hdf5Utilities/hdf5BindDatasets.hpp"



//' Bind matrices by rows or columns
//'
//' Merge existing matrices inside hdf5 data file by rows or by columns
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param datasets, character array indicating the input dataset to be imputed
//' @param outgroup, character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, character array indicating the name for the new merged dataset
//' @param func, character array function to be applyed
//' \describe{
//'     \item{bindRows}{merge datasets by rows}
//'     \item{bindCols}{merge datasets by columns}
// //'     \item{bindRowsbyIndex}{merge datasets by rows taking in to accoutn an index}
//' }
//' @param overwrite, boolean if true, previous results in same location inside hdf5 will be overwritten.
//' @return Original hdf5 data file with results after input datasets
//' @export
// [[Rcpp::export]]
void bdBind_hdf5_datasets( std::string filename, std::string group, Rcpp::StringVector datasets, 
                  std::string outgroup, std::string outdataset, std::string func,
                  Rcpp::Nullable<bool> overwrite = false )
{
    
    try
    {
        
        Rcpp::NumericVector oper = {0, 1, 2};
        oper.names() = Rcpp::CharacterVector({ "bindCols", "bindRows", "bindRowsbyIndex"});
        
        bool boverwrite;
        
        if( overwrite.isNull()) { boverwrite = false; } 
        else {   boverwrite = Rcpp::as<bool>(overwrite); }

        if (func.compare("bindCols") != 0 && func.compare("bindRows") != 0  && func.compare("bindRowsbyIndex") != 0 ) {
            throw std::range_error( "Function to apply must be \"bindRows\", \"bindCols\" or \"bindRowsbyIndex\" other values are not allowed" );
            return void();
        }
        
        std::string stroutDatasetName = outgroup + "/" + outdataset;
        
        
        int bindFunction = oper.findName( func );
        
        
        BigDataStatMeth::hdf5Dataset* dsOut = new BigDataStatMeth::hdf5Dataset(filename, outgroup, outdataset, boverwrite);
        // dsOut->createUnlimitedDataset(count[0], count[1], "real");
        
        RcppBind_datasets_hdf5( filename, group, datasets, dsOut, bindFunction, false);
        
        delete dsOut;
        
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        Rcpp::Rcout<<"c++ exception bdBind_hdf5_datasets (File IException)";
        return void();
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdBind_hdf5_datasets (Group IException)";
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        Rcpp::Rcout << "c++ exception bdBind_hdf5_datasets (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        Rcpp::Rcout << "c++ exception bdBind_hdf5_datasets" << ex.what();
        return void();
    }
    
    // file->close();
    Rcpp::Rcout<< outdataset <<" dataset has been recomposed from blocks\n";
    return void();
    
    
    
    
    // H5File* file = nullptr;
    // DataSet* pdataset = nullptr;
    // DataSet* unlimDataset = nullptr;
    // Rcpp::NumericVector oper = {0, 1, 2};
    // oper.names() = Rcpp::CharacterVector({ "bindCols", "bindRows", "bindRowsbyIndex"});
    // 
    // try
    // {
    //     
    //     IntegerVector count = IntegerVector::create(0, 0);
    //     IntegerVector offset = IntegerVector::create(0, 0);
    //     IntegerVector stride = IntegerVector::create(1, 1);
    //     IntegerVector block = IntegerVector::create(1, 1);
    //     
    //     bool bforce;
    //     
    //     if(force.isNull()) { bforce = false; } 
    //     else {   bforce = Rcpp::as<bool>(force); }
    //     
    //     
    //     // Test file
    //     if( ResFileExist_filestream(filename) ) {
    //         file = new H5File( filename, H5F_ACC_RDWR ); 
    //     } else {
    //         Rcpp::Rcout<<"\nFile not exits, create file before bind matrices";
    //         return void();
    //     }
    //     
    //     // Seek all datasets to perform calculus
    //     for( int i=0; i < datasets.size(); i++ ) 
    //     {
    //         std::string strdataset = group +"/" + datasets(i);
    //         std::string stroutDatasetName = outgroup + "/" + outdataset;
    //         
    //         if( exists_HDF5_element_ptr(file, strdataset ) == 0 ) {
    //             file->close();
    //             Rcpp::Rcout<<"Group or dataset does not exists, please create the input dataset before proceed";
    //             return void();
    //         }
    //         
    //         pdataset = new DataSet(file->openDataSet(strdataset));
    //         
    //         // Real data set dimension
    //         IntegerVector dims_out = get_HDF5_dataset_size(*pdataset);
    //         
    //         // Get block from complete matrix
    //         Eigen::MatrixXd original = GetCurrentBlock_hdf5( file, pdataset, 0, 0, dims_out[0], dims_out[1]);
    //         
    //         // Remove dataset if exists ( only if force = TRUE )
    //         if(i==0) {
    //             prepare_outGroup(file, outgroup, bforce);
    //             prepare_outDataset(file, outgroup + "/" + outdataset, bforce);
    //         }
    //         
    //         if( oper.findName( func ) == 0 || oper.findName( func ) == 1) {
    //             
    //             if(oper.findName( func ) == 0 ){
    //                 
    //                 // Test if dimmensions are correct
    //                 if( original.cols() != count[1] && i!=0) {
    //                     // Append needed cols to merge by cols
    //                     int iappend = count[1] - original.cols();
    //                     original.conservativeResize(original.rows(), original.cols() + iappend);
    //                 }
    //                 offset[0] = offset[0] + count[0];
    //             } else {
    //                 
    //                 // Test if dimmensions are correct
    //                 if( original.rows() != count[0]  && i!=0) {
    //                     // Append needed rows to merge by rows
    //                     int iappend = count[0] - original.rows();
    //                     original.conservativeResize(original.rows() + iappend, original.cols());
    //                 }
    //                 offset[1] = offset[1] + count[1];
    //             }
    //             
    //             count[0] = original.rows();
    //             count[1] = original.cols();
    //             
    //             if(i == 0) {
    //                 // If dataset exists --> remove dataset
    //                 if( exists_HDF5_element_ptr(file,stroutDatasetName))
    //                     remove_HDF5_element_ptr(file,stroutDatasetName);
    //                 // Create unlimited dataset in hdf5 file
    //                 create_HDF5_unlimited_matrix_dataset_ptr(file, stroutDatasetName, count[0], count[1], "numeric");
    //             }
    //             
    //             unlimDataset = new DataSet(file->openDataSet(stroutDatasetName));
    //             
    //             
    //             if(oper.findName( func ) == 0 && i!=0) {
    //                 extend_HDF5_matrix_subset_ptr(file, unlimDataset, count[0], 0);
    //                 
    //             } else if (oper.findName( func ) == 1 && i!=0) {
    //                 extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
    //                 
    //             }
    //             
    //             write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(original)  );  
    //             unlimDataset->close();
    //             pdataset->close();
    //             
    //         } else {
    //             pdataset->close();
    //             file->close();
    //             Rcpp::Rcout<<"Group not exists, create the input dataset before proceed";
    //             return void();
    //         }
    //     }
    // }
    // catch( FileIException& error ) { // catch failure caused by the H5File operations
    //     unlimDataset->close();
    //     pdataset->close();
    //     file->close();
    //     Rcpp::Rcout<<"c++ exception (File IException)";
    //     return void();
    // }
    // 
    // file->close();
    // Rcpp::Rcout<<outdataset<<" dataset has been recomposed from blocks\n";
    // return void();
    
}



/***R

library(BigDataStatMeth)
library(rhdf5)
library(data.table)

setwd("/Users/mailos/DOCTORAT_Local/BigDataStatMeth_Analysis/cca/")

# devtools::reload(pkgload::inst("BigDataStatMeth"))

# Prepare data and functions
X <- matrix(rnorm(150), 50, 3)
Y <- matrix(rnorm(250), 50, 5)


# Create hdf5 data file with  data (Y)
bdCreate_hdf5_matrix_file("cca_cars.hdf5", Y, "data", "Y", force = TRUE)

# Create hdf5 data file with data (X)
bdAdd_hdf5_matrix( X, "cca_cars.hdf5",  "data", "X", force = TRUE)




# Prepare data - Normalize data (only Center)
# 
bdNormalize_hdf5(filename = "cca_cars.hdf5", 
                 group = "data", dataset = "X", 
                 bcenter = TRUE, bscale = FALSE)

bdNormalize_hdf5(filename = "cca_cars.hdf5", 
                 group = "data", dataset = "Y", 
                 bcenter = TRUE, bscale = FALSE)


# Set number of partitions
m <- 10



# # Step 1 :
# # x.block <- splitMatByRow(X, m)
# # y.block <- splitMatByRow(matrix(Y, ncol = 1), m)



# Split datasets X abd Y by rows and store data to data file
bdSplit_matrix_hdf5( filename = "cca_cars.hdf5", 
                     group = "NORMALIZED/data", dataset = "X", 
                     outgroup = "Step1/Xrows", 
                     nblocks = m, bycols = FALSE, force = TRUE)

bdSplit_matrix_hdf5( filename = "cca_cars.hdf5", 
                     group = "NORMALIZED/data", dataset = "Y", 
                     outgroup = "Step1/Yrows", 
                     nblocks = m, bycols = FALSE, force = TRUE)


# Step 2 :

# Get splitted dataset names
x.blocks <- BigDataStatMeth::bdgetDatasetsList_hdf5("cca_cars.hdf5", "Step1/Xrows")
bdapply_Function_hdf5( filename = "cca_cars.hdf5", group = "Step1/Xrows", 
                       datasets = x.blocks, 
                       outgroup = "Step2/Xrows", 
                       func = "QR", 
                       force = TRUE )

y.blocks <- BigDataStatMeth::bdgetDatasetsList_hdf5("cca_cars.hdf5", "Step1/Yrows")
bdapply_Function_hdf5( filename = "cca_cars.hdf5", group = "Step1/Yrows", 
                       datasets = y.blocks, 
                       outgroup = "Step2/Yrows", 
                       func = "QR", 
                       force = TRUE )




#  Step 3 :

# Merge R in Rt from X and Y

x.blocks.qr <- bdgetDatasetsList_hdf5("cca_cars.hdf5", "Step2/Xrows")
bdBind_hdf5(filename = "cca_cars.hdf5", group = "Step2/Xrows", datasets = x.blocks.qr[which(x.blocks.qr %like% ".R")],
            outgroup = "Step3/merged", outdataset = "XRt", 
            func = "bindRows", force = TRUE )
bdapply_Function_hdf5( "cca_cars.hdf5", "Step3/merged", "XRt", "Step3/Final_QR", "QR", force = TRUE )


y.blocks.qr <- bdgetDatasetsList_hdf5("cca_cars.hdf5", "Step2/Yrows")
bdBind_hdf5(filename = "cca_cars.hdf5", group = "Step2/Yrows", datasets = y.blocks.qr[which(y.blocks.qr %like% ".R")],
            outgroup = "Step3/merged", outdataset = "YRt", 
            func = "bindRows", force = TRUE )
bdapply_Function_hdf5( "cca_cars.hdf5", "Step3/merged", "YRt", "Step3/Final_QR", "QR", force = TRUE )





# Step 4 :

bdSplit_matrix_hdf5("cca_cars.hdf5", "Step3/Final_QR", "XRt.Q", 
                    outgroup = "Step4/splitted", 
                    nblocks = m, 
                    bycols = FALSE, force = TRUE )

bdSplit_matrix_hdf5("cca_cars.hdf5", "Step3/Final_QR", "YRt.Q", 
                    outgroup = "Step4/splitted", 
                    nblocks = m, 
                    bycols = FALSE, force = TRUE )

# Step 5 :

# Get splitted matrices names
tmp <- bdgetDatasetsList_hdf5("cca_cars.hdf5", "Step4/splitted")
X.Rt.Q.divide <- tmp[which(tmp %like% "XRt.Q")]
# multiply previous splitted matrices with Q descomposed matrices from model (X)
bdapply_Function_hdf5(  filename = "cca_cars.hdf5", group = "Step2/Xrows", 
                        datasets = x.blocks.qr[which(x.blocks.qr %like% ".Q")], 
                        outgroup = "Step5", func = "blockmult",
                        b_group = "Step4/splitted", b_datasets = X.Rt.Q.divide,
                        force = TRUE )

tmp <- bdgetDatasetsList_hdf5("cca_cars.hdf5", "Step4/splitted")
Y.Rt.Q.divide <- tmp[which(tmp %like% "YRt.Q")]
# multiply previous splitted matrices with Q descomposed matrices from Y
bdapply_Function_hdf5(  filename = "cca_cars.hdf5", group = "Step2/Yrows", 
                        datasets = y.blocks.qr[which(y.blocks.qr %like% ".Q")], 
                        outgroup = "Step5", func = "blockmult",
                        b_group = "Step4/splitted", b_datasets = Y.Rt.Q.divide,
                        force = TRUE )


# Step 6 : 
#   Merge all blocks to create complete QX and QY

blocks.Q <- bdgetDatasetsList_hdf5("cca_cars.hdf5", "Step5")
bdBind_hdf5(filename = "cca_cars.hdf5", group = "Step5", datasets = blocks.Q[which(blocks.Q %like% "X.")],
            outgroup = "Step6", outdataset = "XQ", 
            func = "bindRows", force = TRUE )

bdBind_hdf5(filename = "cca_cars.hdf5", group = "Step5", datasets = blocks.Q[which(blocks.Q %like% "Y.")],
            outgroup = "Step6", outdataset = "YQ", 
            func = "bindRows", force = TRUE )


# Step 7
#   tQXQY <- crossprod(t(QX), QY)[1:ncol(x), ]

res <- bdCrossprod_hdf5(filename = "cca_cars.hdf5", 
                        group = "Step6", A = "XQ",
                        groupB = "Step6", B = "YQ", 
                        outgroup = "Step7")



# Step 8 : 
# z <- svd( tQXQY )
# 
res <- bdSVD_hdf5(file = "cca_cars.hdf5", 
                  group = "Step7", dataset = "CrossProd_XQxYQ",
                  bcenter = FALSE, bscale = FALSE)


# Step 9 :
#   We can solve data on memory (data is small)
#   Solve( QX[1L:dx, 1L:dx, drop = FALSE], z$u)
#   Solve( QY[1L:dy, 1L:dy, drop = FALSE], z$v)

h5ls(res$file)

h5f = H5Fopen(res$file)
XQ <- h5f$Step6$XQ[1:ncol(X), 1:ncol(X)]
YQ <- h5f$Step6$YQ[1:ncol(Y), 1:ncol(Y)]
XR <- h5f$Step3$Final_QR$XRt.R
YR <- h5f$Step3$Final_QR$YRt.R
u <- h5f$SVD$CrossProd_XQxYQ$u
d <- h5f$SVD$CrossProd_XQxYQ$d
v <- h5f$SVD$CrossProd_XQxYQ$v
h5closeAll()

# Get qr compact (more or less)
XR[lower.tri(XR, diag = F)] <- 0
XQ[upper.tri(XQ, diag = TRUE)] <- 0
XQR <- XR + XQ


YR[lower.tri(YR, diag = F)] <- 0
YQ[upper.tri(YQ, diag = TRUE)] <- 0
YQR <- YR + YQ


xcoef.hdf5 <- bdSolve(XQR, u)
ycoef.hdf5 <- bdSolve(YQR, v)

*/
