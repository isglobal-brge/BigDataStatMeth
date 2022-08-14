#include "include/hdf5_sortMatrices.h"


//' Sort existing dataset 
//'
//' Sort an existing dataset taking in to account a list with sorted positions
//' 
//' @param filename, character array indicating the name of the file to be sorted
//' @param group, character array indicating the input group where the data set 
//' to be sorted is stored.
//' @param datasets, character array indicating the input dataset to be sorted
//' @param outdataset, character array indicating the name for the new sorted 
//' dataset
//' @param blockedSortlist, a list with blocks with sorted positions, see example
//' $`1`
//'                       chr order newOrder Diagonal
//' TCGA-OR-A5J1 TCGA-OR-A5J1     1        1        1
//' TCGA-OR-A5J2 TCGA-OR-A5J2     2        2        1
//' TCGA-OR-A5J3 TCGA-OR-A5J3     3        3        1
//' TCGA-OR-A5J4 TCGA-OR-A5J4     4        4        1
//' 
//' $`2`
//'                       chr order newOrder
//' TCGA-OR-A5J5 TCGA-OR-A5JA    10        5        1
//' TCGA-OR-A5J6 TCGA-OR-A5JB    11        6        1
//' TCGA-OR-A5J7 TCGA-OR-A5JC    12        7        0
//' TCGA-OR-A5J8 TCGA-OR-A5JD    13        8        1
//' 
//' $`3`
//'                       chr order newOrder
//' TCGA-OR-A5J9 TCGA-OR-A5J5     5        9        1
//' TCGA-OR-A5JA TCGA-OR-A5J6     6       10        1
//' TCGA-OR-A5JB TCGA-OR-A5J7     7       11        1
//' TCGA-OR-A5JC TCGA-OR-A5J8     8       12        1
//' TCGA-OR-A5JD TCGA-OR-A5J9     9       13        0
//' 
//' where rownames is the current rowname char is the new rowname, order is the
//' current position and newOrder is the new position
//' @param func, character array function to be applyed
//' \describe{
//'     \item{sortRows}{sort datasets rows}
//'     \item{sortCols}{sort datasets columns}
//' }
//' @param outgroup, optional, character array indicating group where the data 
//' set will be saved after imputation if `outgroup` is NULL, output dataset is 
//' stored in the same input group. 
//' @param force, boolean if true, previous results in same location inside hdf5
//' will be overwritten.
//' @return Original hdf5 data file with sorted dataset
//' @examples
//' 
//' @export
// [[Rcpp::export]]
void bdSort_hdf5_dataset( std::string filename, std::string group, 
                          std::string dataset, std::string outdataset, 
                          Rcpp::List blockedSortlist, std::string func, 
                          Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                          Rcpp::Nullable<bool> force = false )
{
    
    
    H5File* file = nullptr;
    DataSet* pdataset = nullptr;
    DataSet* poutdataset = nullptr;
    Rcpp::NumericVector oper = {0, 1};
    oper.names() = Rcpp::CharacterVector({ "sortRows", "sortCols"});
    std::string strOutgroup;
    
    try
    {
        
        IntegerVector count = IntegerVector::create(0, 0);
        IntegerVector offset = IntegerVector::create(0, 0);
        IntegerVector stride = IntegerVector::create(1, 1);
        IntegerVector block = IntegerVector::create(1, 1);
        
        bool bforce;
        
        if( blockedSortlist.length()<=0 ) {
            Rcpp::Rcout<<"\nList is empty, please create a list with the new sort";
            return void();
        }
        
        
        if( force.isNull() ) { bforce = false; } 
        else { bforce = Rcpp::as<bool>(force); }
        
        if( outgroup.isNull() ) { strOutgroup = group; } 
        else {   strOutgroup = Rcpp::as<std::string>(outgroup); }
        
        // Test file
        if( ResFileExist_filestream(filename) ) {
            file = new H5File( filename, H5F_ACC_RDWR ); 
        } else {
            Rcpp::Rcout<<"\nFile not exits, create file before bind matrices";
            return void();
        }
        
        // Test source dataset
        std::string strsoureDataset = group + "/" + dataset;
        std::string stroutDatasetName = strOutgroup + "/" + outdataset;
        
        // Test if input dataset exists
        if( exists_HDF5_element_ptr(file, group) == 0 ) {
            Rcpp::Rcout<<"\nGroup "<< group<<" not exits, create file and dataset before get inverse of Cholesky\n";
            file->close();
            return void();
        } else {
            if(!exists_HDF5_element_ptr(file, strsoureDataset)) {
                Rcpp::Rcout<<"\n Dataset "<< strsoureDataset <<" not exits, create dataset before get inverse of Cholesky\n";
                file->close();
                return void();
            }
        }
        pdataset = new DataSet(file->openDataSet(strsoureDataset));
        
        
        
        // Test if output dataset exists and remove it if force = TRUE
        // if( exists_HDF5_element_ptr(file, stroutDatasetName) && bforce == false) {
        //     Rcpp::Rcout<<"\n Dataset "<< stroutDatasetName <<"  also exists, please set force = TRUE to overwrite\n";
        //     file->close();
        //     return void();
        // } else if(exists_HDF5_element_ptr(file, outdataset) && bforce == true) {
        //     remove_HDF5_element_ptr(file, stroutDatasetName); 
        // }
        if( !exists_HDF5_element_ptr(file, stroutDatasetName)) {
            Rcpp::Rcout<<"\n Dataset "<< stroutDatasetName <<" not exits, create empty output dataset \n";
            pdataset->close();
            file->close();
            return void();
        }
        poutdataset = new DataSet(file->openDataSet(stroutDatasetName));
        // Real data set dimension
        IntegerVector dims_out = get_HDF5_dataset_size_ptr(pdataset);
        int nrows = dims_out[0];
        int ncols = dims_out[1];
        
        
        for( int i = 0; i < blockedSortlist.length(); i++) {
            
            Rcpp::DataFrame df(blockedSortlist[i]);
            std::vector<double> order = df[1];
            std::vector<double> neworder = df[2];
            std::vector<double> diagonal = df[3];
            
            auto indices_0 = find_all(diagonal.begin(), diagonal.end(), 0);
            
            if( indices_0.size() > 0) {
            
                for(int t=0; t<indices_0.size(); t++){
                    Rcpp::Rcout<<"Indices val : " <<&indices_0[t]<<"\n";    
                }
                
                
            } else {
                
                if( oper.findName( func ) == 0 ) {
                    offset[1] = order[0] - 1;
                    count[0] = order[order.size() - order[0]];
                    count[1] = dims_out[1]; 
                } else if( oper.findName( func ) == 1 ) {
                    offset[0] = order[0] - 1;
                    count[1] = dims_out[1]; 
                    count[0] = order[order.size() - order[0]];
                } 
                
                Rcpp::Rcout<<"\n-> Llegirem des de : "<<offset[0]<<" - "<<offset[1]<<"\n";
                Rcpp::Rcout<<"La quantitat de files - columnes  : "<<count[0]<<" - "<<count[1]<<"\n";
                
                Eigen::MatrixXd A = GetCurrentBlock_hdf5(file, pdataset, offset[0], offset[1], count[0], count[1]);
                
                if( oper.findName( func ) == 0 ) {
                    offset[1] = neworder[0]-1;
                } else if( oper.findName( func ) == 1 ) {
                    offset[0] = neworder[0]-1;
                }
                
                Rcpp::Rcout<<"\n-> Escriurem des de : "<<offset[0]<<" - "<<offset[1]<<"\n";
                Rcpp::Rcout<<"La quantitat de files - columnes  : "<<count[0]<<" - "<<count[1]<<"\n";
                write_HDF5_matrix_subset_v2( file, poutdataset, offset, count, stride, block, Rcpp::wrap( A ) );
                
            }

            // 
            // if( oper.findName( func ) == 0 ) {
            //     offset[1] = order[0];
            //     count[0] = dims_out[1]; 
            //     count[1] = order[order.size() - order[0]];
            //     
            //     Eigen::MatrixXd A = GetCurrentBlock_hdf5(file, pdataset, offset[0], offset[1], count[0], count[1]);
            //     offset[1] = neworder[0];
            //     write_HDF5_matrix_subset_v2( file, poutdataset, offset, count, stride, block, Rcpp::wrap( A ) );
            //     
            // } else if( oper.findName( func ) == 1 ) {
            //     offset[0] = order[0];
            //     count[1] = dims_out[1]; 
            //     count[0] = order[order.size() - order[0]];
            //     
            //     
            //     
            //     Eigen::MatrixXd A = GetCurrentBlock_hdf5(file, pdataset, offset[0], offset[1], count[0], count[1]);
            //     offset[0] = neworder[0];
            //     write_HDF5_matrix_subset_v2( file, poutdataset, offset, count, stride, block, Rcpp::wrap( A ) );
            // }    
            // 
            
        }
        
        
        

    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        poutdataset->close();
        pdataset->close();
        file->close();
        Rcpp::Rcout<<"c++ exception bdSort_hdf5_dataset (File IException)";
        return void();
    } catch( GroupIException & error ) { // catch failure caused by the DataSet operations
        poutdataset->close();
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdSort_hdf5_dataset (Group IException)";
        return void();
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        poutdataset->close();
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdSort_hdf5_dataset (DataSet IException)";
        return void();
    } catch(std::exception& ex) {
        poutdataset->close();
        pdataset->close();
        file->close();
        Rcpp::Rcout << "c++ exception bdSort_hdf5_dataset" << ex.what();
        return void();
    }
    
    pdataset->close();
    poutdataset->close();
    file->close();
    Rcpp::Rcout<<outdataset<<" dataset has been sorted \n";
    return void();
    
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
bdCreate_hdf5_matrix_file("cca_cars.hdf5", Y, "data", "Y", force = T)

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
XQ[upper.tri(XQ, diag = T)] <- 0
XQR <- XR + XQ


YR[lower.tri(YR, diag = F)] <- 0
YQ[upper.tri(YQ, diag = T)] <- 0
YQR <- YR + YQ


xcoef.hdf5 <- bdSolve(XQR, u)
ycoef.hdf5 <- bdSolve(YQR, v)

*/
