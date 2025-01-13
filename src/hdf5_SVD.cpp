#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixSvd.hpp"



//' Block SVD decomposition with HDF5 files
//'
//' This function computes the singular values and left singular vectors of a 
//' real nxp matrix using Block SVD decomposition with an incremental algorithm. 
//' The input matrix is read from an HDF5 file, and the results are saved in the 
//' same file.
//' 
//' @inheritParams bdNormalize_hdf5
//' @param k numerical, number of local SVDs to concatenate at each level. 
//' Defaults is set to 2.
//' This parameter helps optimize the performance and memory usage during PCA 
//' calculations. 
//' @param q numerical, number of levels to compute SVD for PCA.
//' This parameter helps optimize the performance and memory usage during PCA 
//' calculations. 
//' @param bcenter logical (optional). If TRUE (default), the data is centered 
//' by subtracting the column means (ignoring NAs) of the `dataset` from their 
//' corresponding columns. If FALSE, no centering is performed.
//' @param bscale (optional). If TRUE (default), the data is scaled by dividing 
//' the (centered) columns of `x` by their standard deviations if `bcenter` 
//' is TRUE, or by the root mean square otherwise. If FALSE, no scaling is 
//' performed.
//' @param rankthreshold `double`. Threshold used to determine the range of 
//' the matrix. The matrix rank is defined as the number of singular values that 
//' differ from the threshold. By default, `threshold = 0` is used to compute 
//' the matrix rank, but it can be adjusted to a value close to zero for 
//' approximations.
//' @param overwrite logical value, If TRUE, forces the recalculation of results 
//' even if they already exist.
//' @param method optional, defalut is "auto" possible values are: "auto", 
//' "blocks", "full":
//'     * `"auto"`:
//'       The option method = "auto" chooses the "full" or 
//'       "blocks" method depending on the size of the matrix to be decomposed 
//'     * `"blocks"`:
//'       The SVD decomposition can be carried out by blocks, recommended option 
//'       for large matrices that do not fit in memory
//'     * `"full"`:
//'       The SVD decomposition is performed directly without partitioning the matrix
//' @param threads integer (optional), an optional parameter specifying the 
//' number of threads to use.
//' @return three dataset inside HDF5 data files with the singular values and 
//' left and right singular vectors of the dataset:
//' 
//'   * `"u"`:
//'     eigenvectors of AA^t, mxn and column orthogonal matrix 
//'   * `"v`:
//'     eigenvectors of A^tA, nxn orthogonal matrix
//'   * `"d"`:
//'     singular values, nxn diagonal matrix (non-negative real values) 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSVD_hdf5 ( Rcpp::RObject filename, Rcpp::Nullable<Rcpp::CharacterVector> group = R_NilValue, 
                       Rcpp::Nullable<Rcpp::CharacterVector> dataset = R_NilValue,
                       Rcpp::Nullable<int> k=2, Rcpp::Nullable<int> q=1,
                       Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true,
                       Rcpp::Nullable<double> rankthreshold = 0.0,
                       Rcpp::Nullable<bool> overwrite = R_NilValue,
                       Rcpp::Nullable<Rcpp::CharacterVector> method = R_NilValue,
                       Rcpp::Nullable<int> threads = R_NilValue)
{
 
     std::string str_filename;
     double dthreshold;
     
     try {
         
         int ks, qs, nvs = 0;
         bool bcent, bscal, bforce, bRowMajor = false; //, bbyblocks = true;
         Rcpp::CharacterVector strgroup, strdataset;
         
         if(k.isNull())  ks = 2 ;
         else    ks = Rcpp::as<int>(k);
         
         if(q.isNull())  qs = 1 ;
         else    qs = Rcpp::as<int>(q);
         
         if(bcenter.isNull())  bcent = true ;
         else    bcent = Rcpp::as<bool>(bcenter);
         
         if(bscale.isNull())  bscal = true ;
         else    bscal = Rcpp::as<bool>(bscale);
         
         if(overwrite.isNull())  bforce = false ;
         else    bforce = Rcpp::as<bool>(overwrite);
         
         if(group.isNull())  strgroup = "" ;
         else    strgroup = Rcpp::as<std::string>(group);
         
         if(dataset.isNull())  strdataset = "";
         else    strdataset = Rcpp::as<std::string>(dataset);
         
         if(Rcpp::is<Rcpp::CharacterVector>(filename)) {
             str_filename = Rcpp::as<std::string>(filename);
         } else {
             Rcpp::Rcout<< "File name must be character string";
             return Rcpp::List::create(Rcpp::Named("file") = "");
         }
         
         if(rankthreshold.isNull()) {  
             dthreshold = 0 ;
         } else {
             if( Rcpp::as<double>(rankthreshold) > 0.1 ) {
                 Rcpp::Rcout<< "Threshold to big, please set threshold with value lower than 0.1";
                 return Rcpp::List::create(Rcpp::Named("file") = str_filename);
             } else if( Rcpp::as<double>(rankthreshold) < 0 ) {
                 Rcpp::Rcout<< "Threshold must be a positive value near zero";
                 return Rcpp::List::create(Rcpp::Named("file") = str_filename);
             } else {
                 dthreshold = Rcpp::as<double>(rankthreshold);
             }
         }
         
         // retsvd = BigDataStatMeth::RcppbdSVD_hdf5( filename, Rcpp::as<std::string>(strgroup), Rcpp::as<std::string>(strdataset), ks, qs, nvs, bcent, bscal, dthreshold, threads );
         BigDataStatMeth::RcppbdSVD_hdf5( str_filename, Rcpp::as<std::string>(strgroup), Rcpp::as<std::string>(strdataset), ks, qs, nvs, bcent, bscal, dthreshold, bforce, bRowMajor, method, threads );
         
     } catch(std::exception &ex) {
         Rcpp::Rcout<<"c++ exception bdSVD_hdf5 \n"<< ex.what();
         return Rcpp::List::create(Rcpp::Named("file") = R_NilValue);
     }
     
     return Rcpp::List::create(Rcpp::Named("file") = str_filename);
 
}

/**
 //' @param file a real nxp matrix in hdf5 file
 //' @param group group in hdf5 data file where dataset is located
 //' @param dataset matrix dataset with data to perform SVD
 //' @param k number of local SVDs to concatenate at each level 
 //' @param q number of levels
 //' @param rankthreshold double, threshold used to determine the range of the array. 
 //' The matrix rank is equal to the number of singular values different from the 
 //' threshold. By default, threshold = 0 is used to get the matrix rank , but it 
 //' can be changed to an approximation of 0.
 */