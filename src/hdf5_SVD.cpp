#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixSvd.hpp"



//' Block SVD decomposition for hdf5 files using an incremental algorithm.
//'
//' Singular values and left singular vectors of a real nxp matrix 
//' Block SVD decomposition using an incremental algorithm.
//' @param file a real nxp matrix in hdf5 file
//' @param group group in hdf5 data file where dataset is located
//' @param dataset matrix dataset with data to perform SVD
//' @param k number of local SVDs to concatenate at each level 
//' @param q number of levels
//' @param bcenter (optional, defalut = TRUE) . If center is TRUE then centering 
//' is done by subtracting the column means (omitting NAs) of x from their 
//' corresponding columns, and if center is FALSE, no centering is done.
//' @param bscale (optional, defalut = TRUE) .  If scale is TRUE then scaling is 
//' done by dividing the (centered) columns of x by their standard deviations if 
//' center is TRUE, and the root mean square otherwise. If scale is FALSE, no 
//' scaling is done.
//' @param rankthreshold double, threshold used to determine the range of the array. 
//' The matrix rank is equal to the number of singular values different from the 
//' threshold. By default, threshold = 0 is used to get the matrix rank , but it 
//' can be changed to an approximation of 0.
//' @param method (optional, defalut = "auto") possible values are: "auto", 
//' "blocks", "full":
//' \itemize{
//'   \item{"auto"}{ The option method = "auto" chooses the "full" or 
//'   "blocks" method depending on the size of the matrix to be decomposed }
//'   \item{"blocks"}{ The SVD decomposition can be carried out by blocks, 
//'   recommended option for large matrices that do not fit in memory }
//'   \item{"full"}{ The SVD decomposition is performed directly without partitioning the matrix }
//' } 
//' @param threads (optional) only used in some operations inside function. If 
//' threads is null then threads =  maximum number of threads available - 1.
//' @return a list of three components with the singular values and left and 
//' right singular vectors of the matrix
//' @return A List with : 
//' \itemize{
//'   \item{"u"}{ eigenvectors of AA^t, mxn and column orthogonal matrix }
//'   \item{"v"}{ eigenvectors of A^tA, nxn orthogonal matrix }
//'   \item{"d"}{ singular values, nxn diagonal matrix (non-negative real values) }
//' }
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSVD_hdf5 ( Rcpp::RObject file, Rcpp::Nullable<Rcpp::CharacterVector> group = R_NilValue, 
                       Rcpp::Nullable<Rcpp::CharacterVector> dataset = R_NilValue,
                       Rcpp::Nullable<int> k=2, Rcpp::Nullable<int> q=1,
                       Rcpp::Nullable<bool> bcenter=true, Rcpp::Nullable<bool> bscale=true,
                       Rcpp::Nullable<double> rankthreshold = 0.0,
                       Rcpp::Nullable<bool> force = R_NilValue,
                       Rcpp::Nullable<Rcpp::CharacterVector> method = R_NilValue,
                       Rcpp::Nullable<int> threads = R_NilValue)
{
 
     std::string filename;
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
         
         if(force.isNull())  bforce = false ;
         else    bforce = Rcpp::as<bool>(force);
         
         if(group.isNull())  strgroup = "" ;
         else    strgroup = Rcpp::as<std::string>(group);
         
         if(dataset.isNull())  strdataset = "";
         else    strdataset = Rcpp::as<std::string>(dataset);
         
         if(Rcpp::is<Rcpp::CharacterVector>(file)) {
             filename = Rcpp::as<std::string>(file);
         } else {
             Rcpp::Rcout<< "File name must be character string";
             return Rcpp::List::create(Rcpp::Named("file") = "");
         }
         
         if(rankthreshold.isNull()) {  
             dthreshold = 0 ;
         } else {
             if( Rcpp::as<double>(rankthreshold) > 0.1 ) {
                 Rcpp::Rcout<< "Threshold to big, please set threshold with value lower than 0.1";
                 return Rcpp::List::create(Rcpp::Named("file") = filename);
             } else if( Rcpp::as<double>(rankthreshold) < 0 ) {
                 Rcpp::Rcout<< "Threshold must be a positive value near zero";
                 return Rcpp::List::create(Rcpp::Named("file") = filename);
             } else {
                 dthreshold = Rcpp::as<double>(rankthreshold);
             }
         }
         
         // retsvd = BigDataStatMeth::RcppbdSVD_hdf5( filename, Rcpp::as<std::string>(strgroup), Rcpp::as<std::string>(strdataset), ks, qs, nvs, bcent, bscal, dthreshold, threads );
         BigDataStatMeth::RcppbdSVD_hdf5( filename, Rcpp::as<std::string>(strgroup), Rcpp::as<std::string>(strdataset), ks, qs, nvs, bcent, bscal, dthreshold, bforce, bRowMajor, method, threads );
         
     } catch(std::exception &ex) {
         Rcpp::Rcout<<"c++ exception bdSVD_hdf5 \n"<< ex.what();
         return Rcpp::List::create(Rcpp::Named("file") = R_NilValue);
     }
     
     return Rcpp::List::create(Rcpp::Named("file") = filename);
 
}

