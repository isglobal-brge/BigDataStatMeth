#ifndef BIGDATASTATMETH_HDF5_MATRIXSVD_HPP
#define BIGDATASTATMETH_HDF5_MATRIXSVD_HPP

#include <RcppEigen.h>
#include "Utilities/Utilities.hpp"
#include "hdf5Algebra/matrixNormalization.hpp"
#include "H5Cpp.h"
#include "matrixSvdBlock.hpp"

namespace BigDataStatMeth {

extern "C" {
    extern void dgesdd_( char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*, int*);
}


// Lapack SVD decomposition - Optimized algorithm with dgesdd
template <class T>
extern svdeig RcppbdSVD_lapack( T X, bool bcenter, bool bscale, bool complete ) {

    svdeig retsvd;
    
    char Schar='S';
    char Achar='A';
    int info = 0;
    
    try {
        
        if(bcenter == true || bscale == true) {
            X = RcppNormalize_Data(X, bcenter, bscale);
        }
        
        int m = X.rows();
        int n = X.cols();
        int lda = std::max(1,m);
        int ldu = std::max(1,m);
        int ldvt = std::min(m, n);
        int k = std::min(m,n);
        int lwork;
        
        if( complete == false ) {
            lwork = 4*std::min(m,n)*std::min(m,n) + 7*std::min(m,n);
        } else {
            lwork = 4*std::min(m,n)*std::min(m,n) + 6*std::min(m,n) + std::max(m,n);
        }
        
        Eigen::VectorXd s = Eigen::VectorXd::Zero(k);
        Eigen::VectorXd work = Eigen::VectorXd::Zero(lwork);
        Eigen::VectorXi iwork(8*std::min(m,n));
        Eigen::MatrixXd u;
        Eigen::MatrixXd vt = Eigen::MatrixXd::Zero(ldvt,n);
        
        if( complete == false ) {
            u = Eigen::MatrixXd::Zero(ldu,k);
            dgesdd_( &Schar, &m, &n, X.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);
        } else {
            u = Eigen::MatrixXd::Zero(ldu,m);
            dgesdd_( &Achar, &m, &n, X.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);
        }
        
        retsvd.d = s;
        retsvd.u = u;
        retsvd.v = vt.transpose();
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< "C++ exception RcppbdSVD_lapack : "<< ex.what();
        return retsvd;
    } catch (...) {
        ::Rf_error("C++ exception RcppbdSVD_lapack (unknown reason)");
        return retsvd;
    }
    
    return retsvd;
    
}



// ##' @param k number of local SVDs to concatenate at each level 
// ##' @param q number of levels
extern svdeig RcppbdSVD_hdf5_Block( BigDataStatMeth::hdf5Dataset* dsA, int k, int q, int nev, bool bcenter, bool bscale, 
                             int irows, int icols, double dthreshold, Rcpp::Nullable<int> threads = R_NilValue )
{
    
    Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1, 1);
    Rcpp::IntegerVector block = Rcpp::IntegerVector::create(1, 1);
    svdeig retsvd;
    Eigen::MatrixXd nX;
    bool transp = false;
    std::string strGroupName  = "tmpgroup";
    std::string strPrefix;
    
    Rcpp::CharacterVector strvmatnames = {"A","B","C","D","E","F","G","H","I","J","K",
                                    "L","M","N","O","P","Q","R","S","T","U","V",
                                    "W","X","Y","Z"};
    strPrefix = strvmatnames[q-1];
    
    try{
        
        if(irows >= icols)
            transp = true;
        
        First_level_SvdBlock_decomposition_hdf5( dsA, k, q, nev, bcenter, bscale, irows, icols, dthreshold, threads);
 
 /***        
        for(int j = 1; j < q; j++) { // For each decomposition level :
            Next_level_SvdBlock_decomposition_hdf5(file, strGroupName, k, j, bcenter, bscale, dthreshold, threads);
        }
        
        // Get dataset names
        StringVector joindata =  get_dataset_names_from_group(file, strGroupName, strPrefix);
        
        // 1.- Join matrix and remove parts from file
        std::string strnewdataset = std::string((joindata[0])).substr(0,1);
        join_datasets(file, strGroupName, joindata, strnewdataset);
        remove_HDF5_multiple_elements_ptr(file, strGroupName, joindata);
        
        // 2.- Get SVD from Blocks full mattrix
        DataSet datasetlast = file->openDataSet(strGroupName + "/" + strnewdataset);
        IntegerVector dims_out = get_HDF5_dataset_size(datasetlast);
        
        Eigen::MatrixXd matlast;
        matlast = GetCurrentBlock_hdf5(file, &datasetlast, 0, 0, dims_out[0],dims_out[1]);
        
        
        if( nev < std::min( int(matlast.rows()), int(matlast.cols()) )) {
            retsvd = RcppbdSVD(matlast, nev, 0, false, false);
        } else {
            retsvd = RcppbdSVD_lapack(matlast, false, false, false);
        }
        
        
        // Write results to hdf5 file : in folder "SVD" and dataset "SVD".<name input dataset>
        // Create structure and write d 
        StringVector name = get_dataset_names_from_dataset_ptr(dataset);
        
        create_HDF5_groups_ptr(file,"SVD/"+ name[0]);
        write_HDF5_matrix_ptr(file, "SVD/"+ name[0]+"/d", wrap(retsvd.d));
        
        Eigen::MatrixXd v;    
        
        
        // 3.- crossprod initial matrix and svdA$u
        
        Eigen::MatrixXd A;
        
        
        if( bcenter == true || bscale == true) {
            DataSet* normalizedData = nullptr;
            normalizedData = new DataSet(file->openDataSet(strGroupName + "/normalmatrix"));
            IntegerVector dims_out_normal = get_HDF5_dataset_size(*normalizedData);
            A = GetCurrentBlock_hdf5_Original(file, normalizedData, 0, 0,dims_out_normal[0], dims_out_normal[1] );
            
            normalizedData->close();
            delete(normalizedData);
        } else {
            
            IntegerVector dims_out_normal = get_HDF5_dataset_size(*dataset);
            A = GetCurrentBlock_hdf5(file, dataset, 0, 0,dims_out_normal[0], dims_out_normal[1] );
        }
        
        
        v = Bblock_matrix_mul_parallel(A, retsvd.u, 1024, threads); //  PARALLEL ==> NOT PARALLEL
        
        // 4.- resuls / svdA$d
        v = v.array().rowwise()/(retsvd.d).transpose().array();
        
        
        
        if (transp == true)  {
            write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/u", wrap(retsvd.v));
            write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/v", wrap(retsvd.u));
        } else {
            write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/u", wrap(retsvd.u));
            write_HDF5_matrix_transposed_ptr(file, "SVD/"+ name[0]+"/v", wrap(retsvd.v));
        }
        
        
        // Clean data
        if( bcenter == true || bscale == true) {
            remove_HDF5_multiple_elements_ptr(file, strGroupName, "normalmatrix");
            remove_HDF5_element_ptr(file, strGroupName);
        }
 ***/       
        
    } ccatch(std::exception &ex) {
        Rcpp::Rcout<< "C++ exception RcppbdSVD_lapack : "<< ex.what();
        return retsvd;
    } catch (...) {
        ::Rf_error("C++ exception RcppbdSVD_lapack (unknown reason)");
        return retsvd;
    }
    
    return retsvd;
}











// SVD decomposition with hdf5 file
//    input data : hdf5 file (object from crossproduct matrix) datagroup = 'strsubgroupIN'
//    output data : hdf5 file svd data in datagroup svd 
//                        svd/d 
//                        svd/u 
//                        svd/v 
//                        
//  https://github.com/isglobal-brge/svdParallel/blob/8b072f79c4b7c44a3f1ca5bb5cba4d0fceb93d5b/R/generalBlockSVD.R
//  @param k number of local SVDs to concatenate at each level 
//  @param q number of levels
//  
// extern svdeig RcppbdSVD_hdf5( std::string filename, std::string strsubgroup, std::string strdataset,  
//                        int k, int q, int nev, bool bcenter, bool bscale, double dthreshold, 
//                        Rcpp::Nullable<int> ithreads = R_NilValue )
extern void RcppbdSVD_hdf5( std::string filename, std::string strsubgroup, std::string strdataset,  
                              int k, int q, int nev, bool bcenter, bool bscale, double dthreshold, 
                              bool bforce, bool asRowMajor, Rcpp::Nullable<int> ithreads = R_NilValue)
{
    
    try {
        
        svdeig retsvd;
        Eigen::MatrixXd X;
        hdf5Dataset* dsu;
        hdf5Dataset* dsv;
        hdf5Dataset* dsd;
        
        std::vector<hsize_t> stride = {1, 1},
                             block = {1, 1},
                             offset = {0, 0},
                             count = {0, 0};
        
        hdf5Dataset* dsA = new hdf5Dataset(filename, strsubgroup, strdataset, false);
        dsA->openDataset();
        
        std::vector<hsize_t> dims_out = {dsA->nrows(), dsA->ncols()};;
        count = { dims_out[0], dims_out[1]};
        
        // Small matrices ==> Direct SVD (lapack)
        if( dims_out[0] < MAXSVDBLOCK &&  dims_out[1] < MAXSVDBLOCK ) {
            
            std::vector<double> vdA( count[0] * count[1] ); 
            dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
            
            if(asRowMajor == true) {
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdA.data(), count[0], count[1] );    
                retsvd = RcppbdSVD_lapack(X, bcenter, bscale, false);
            } else {
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> X (vdA.data(), count[1], count[0] );
                retsvd = RcppbdSVD_lapack(X, bcenter, bscale, false);
            }
            
            std::string stroutgroup = "SVD/"+ strdataset;
            
            dsu = new hdf5Dataset(filename, stroutgroup, "u", bforce);
            dsu->createDataset( retsvd.u.rows(), retsvd.u.cols(), "real");
            dsu->writeDataset( Rcpp::wrap(retsvd.u) );
            
            dsv = new hdf5Dataset(filename, stroutgroup, "v", bforce);
            dsv->createDataset( retsvd.v.rows(), retsvd.v.cols(), "real");
            dsv->writeDataset( Rcpp::wrap(retsvd.v) );
            
            dsd = new hdf5Dataset(filename, stroutgroup, "d", bforce);
            dsd->createDataset( retsvd.d.size(), 1, "real");
            dsd->writeDataset( Rcpp::wrap(retsvd.d) );
            
        } else {

/***            
            // data stored transposed in hdf5
            int xdim = (unsigned long long)dims_out[1];
            int ydim = (unsigned long long)dims_out[0];
            
            retsvd = RcppbdSVD_hdf5_Block( file, dataset, k, q, nev, bcenter, bscale, xdim, ydim, dthreshold, wrap(ithreads));
***/
            
        }
        
        delete dsu;
        delete dsv;
        delete dsd;
        delete dsA;
        
    }  catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception RcppbdSVD_hdf5 (File IException)" );
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception RcppbdSVD_hdf5 (DataSet IException)" );
        return void();
    } catch(std::exception &ex) {
        Rcpp::Rcout<<"c++ exception RcppbdSVD_hdf5 \n"<< ex.what();
        return void();
    }
    return void();
}


}

#endif // BIGDATASTATMETH_HDF5_MATRIXSVD_HPP

