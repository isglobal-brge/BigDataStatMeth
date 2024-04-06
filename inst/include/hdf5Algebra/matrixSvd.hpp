#ifndef BIGDATASTATMETH_HDF5_MATRIXSVD_HPP
#define BIGDATASTATMETH_HDF5_MATRIXSVD_HPP

// #include <RcppEigen.h>
// #include "H5Cpp.h"

#include "hdf5Algebra/matrixNormalization.hpp"
#include "hdf5Algebra/matrixSvdBlock.hpp"
#include "hdf5Algebra/multiplication.hpp"
#include "memAlgebra/memOptimizedProducts.hpp"
#include "memAlgebra/memMultiplication.hpp"

#include "hdf5Utilities/hdf5Utilities.hpp"
#include "hdf5Utilities/hdf5Methods.hpp"

#include "Utilities/Utilities.hpp"
#include <spectra/SymEigsSolver.h>


namespace BigDataStatMeth {

    
    // SVD decomposition 
    extern inline svdeig RcppbdSVD( Eigen::MatrixXd& X, int k, int ncv, bool bcenter, bool bscale )
    {
        
        svdeig retsvd;
        Eigen::MatrixXd nX;
        int nconv;
        
        if( k==0 )    k = (std::min(X.rows(), X.cols()))-1;
        else if (k > (std::min(X.rows(), X.cols()))-1 ) k = (std::min(X.rows(), X.cols()))-1;
        
        if(ncv == 0)  ncv = k + 1 ;
        if(ncv<k) ncv = k + 1;
        
        {
            Eigen::MatrixXd Xtcp;
            if(bcenter ==true || bscale == true)  {
                nX = RcppNormalize_Data(X, bcenter, bscale, false);
                Xtcp =  bdtcrossproduct(nX);
            } else {
                Xtcp =  bdtcrossproduct(X);
            }
            
            Spectra::DenseSymMatProd<double> op(Xtcp);
            Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, ncv);
            // Spectra::SymEigsSolver< Spectra::DenseSymMatProd<double> > eigs(op, k, ncv);
            
            // Initialize and compute
            eigs.init();
            nconv = eigs.compute();
            
            // if(eigs.info() == Spectra::CompInfo::Successful) {
            if(eigs.info() == Spectra::SUCCESSFUL) {
                retsvd.d = eigs.eigenvalues().cwiseSqrt();
                retsvd.u = eigs.eigenvectors();
                retsvd.bokuv = true;
            } else {
                retsvd.bokuv = false;
            }
        }
        if(retsvd.bokuv == true)
        {
            Eigen::MatrixXd Xcp;
            if(bcenter ==true || bscale==true ) {
                Xcp =  bdcrossproduct(nX);  
            } else {
                Xcp =  bdcrossproduct(X);
            }  
            
            Spectra::DenseSymMatProd<double> opv(Xcp);
            Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigsv(&opv, k, ncv);
            // Spectra::SymEigsSolver< Spectra::DenseSymMatProd<double> > eigsv(opv, k, ncv);
            
            // Initialize and compute
            eigsv.init();
            nconv = eigsv.compute();
            
            // Retrieve results
            // if(eigsv.info() == Spectra::CompInfo::Successful) {
            if(eigsv.info() == Spectra::SUCCESSFUL) {
                retsvd.v = eigsv.eigenvectors();
            } else {
                retsvd.bokd = false;
            }
        }
        
        return retsvd;
    }
    
    
    
    // ##' @param k number of local SVDs to concatenate at each level 
    // ##' @param q number of levels
    extern inline void RcppbdSVD_hdf5_Block( BigDataStatMeth::hdf5Dataset* dsA, 
                                        BigDataStatMeth::hdf5Dataset* dsu, BigDataStatMeth::hdf5Dataset* dsv, 
                                        BigDataStatMeth::hdf5Dataset* dsd, int k, int q, int nev, bool bcenter, bool bscale, 
                                        int irows, int icols, double dthreshold, Rcpp::Nullable<int> threads = R_NilValue )
    {
        
        try{
            
            std::vector<hsize_t> stride = {1, 1},
                block = {1, 1};
            
            svdeig retsvd;
            Eigen::MatrixXd nX;
            Eigen::MatrixXd matlast;
            Eigen::MatrixXd v;    
            bool transp = false;
            std::string strGroupName  = "tmpgroup",
                        strPrefix; // Outpu group name for temporal data
            
            Rcpp::CharacterVector strvmatnames = {"A","B","C","D","E","F","G","H","I","J","K",
                                                  "L","M","N","O","P","Q","R","S","T","U","V",
                                                  "W","X","Y","Z"};
            strPrefix = strvmatnames[q-1];
            
            typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
            
            if(irows >= icols) {
                transp = true;
            }
            
            First_level_SvdBlock_decomposition_hdf5( dsA, strGroupName, k, q, nev, bcenter, bscale, irows, icols, dthreshold, threads);
            
            for(int j = 1; j < q; j++) { // For each decomposition level :
                Next_level_SvdBlock_decomposition_hdf5( dsA, strGroupName, k, j, dthreshold, threads);
            }
            
            // Get dataset names
            Rcpp::StringVector joindata =  dsA->getDatasetNames(strGroupName, strPrefix);
            
            // 1.- Join matrix and remove parts from file
            std::string strnewdataset = std::string((joindata[0])).substr(0,1);
            hdf5Dataset* dsJoined = new hdf5DatasetInternal(dsA->getFileName(), strGroupName, strnewdataset, true);
            join_datasets( dsJoined, strGroupName, joindata, false, true );
            
            // 2.- Get SVD from Blocks full mattrix
            hsize_t* dims_out = dsJoined->dim();
            
            std::vector<double> vdreaded( dims_out[0] * dims_out[1] ); 
            dsJoined->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, {1, 1}, {1, 1}, vdreaded.data() );
            
            
            matlast = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdreaded.data(), dims_out[0], dims_out[1] );
            delete dsJoined;
            
            retsvd = RcppbdSVD_lapack(matlast, false, false, false);
            
            // Write results to hdf5 file : in folder "SVD" and dataset "SVD".<name input dataset>
            dsd->createDataset( 1, retsvd.d.size(), "real");
            dsd->writeDataset( Rcpp::wrap(retsvd.d) );
            
            // 3.- crossprod initial matrix and svdA$u
            
            if( bcenter == true || bscale == true || (dsA->getGroupName().find("NORMALIZED_T") != std::string::npos) ) {
                
                if(bcenter == true || bscale == true) {
                    
                    BigDataStatMeth::hdf5Dataset* normalizedData = new BigDataStatMeth::hdf5Dataset( dsA->getFileName(), strGroupName, "normalmatrix", false);
                    normalizedData->openDataset();
                    
                    dims_out = normalizedData->dim();
                    
                    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dims_out[1], dims_out[0]);
                    normalizedData->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, stride, block, A.data() );
                    
                    delete normalizedData ;
                    
                    v = Rcpp_block_matrix_mul_parallel(A, retsvd.u, R_NilValue, threads); //  PARALLEL ==> NOT PARALLEL
                    
                } else {
                    
                    BigDataStatMeth::hdf5DatasetInternal* normalizedData = new BigDataStatMeth::hdf5DatasetInternal( dsA->getFileName(), dsA->getGroupName(), dsA->getDatasetName(), false);
                    normalizedData->openDataset();
                    
                    dims_out = normalizedData->dim();
                    
                    std::vector<double> vdA( dims_out[0] * dims_out[1] );
                    normalizedData->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, stride, block, vdA.data() );
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A (vdA.data(), dims_out[0], dims_out[1]);
                    delete normalizedData ;
                    
                    v = Rcpp_block_matrix_mul_parallel(A, retsvd.u, R_NilValue, threads); //  PARALLEL ==> NOT PARALLEL
                }
                
            } else {
                
                dims_out = dsA->dim();    
                
                Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dims_out[0], dims_out[1]);
                dsA->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, stride, block, A.data() );
                
                A.transposeInPlace();
                v = Rcpp_block_matrix_mul_parallel(A, retsvd.u, R_NilValue, threads); //  PARALLEL ==> NOT PARALLEL
                
            }
                
            // 4.- resuls / svdA$d
            v = v.array().rowwise()/(retsvd.d).transpose().array();
            
            if (transp == false)  {
                dsu->createDataset( v.rows(), v.cols(), "real");
                dsv->createDataset( retsvd.u.rows(), retsvd.u.cols(), "real");
                
                dsu->writeDataset(Rcpp::wrap(v));
                dsv->writeDataset(Rcpp::wrap(retsvd.u));
            } else {
                dsu->createDataset( retsvd.u.rows(), retsvd.u.cols(), "real");
                dsv->createDataset( v.rows(), v.cols(), "real");
                
                dsu->writeDataset(Rcpp::wrap(retsvd.u));
                dsv->writeDataset(Rcpp::wrap(v));
            }
            
            remove_elements(dsA->getFileptr(), strGroupName);
            
        } catch(std::exception &ex) {
            Rcpp::Rcout<< "C++ exception RcppbdSVD_hdf5_Block : "<< ex.what();
            return void();
        } catch (...) {
            ::Rf_error("C++ exception RcppbdSVD_hdf5_Block (unknown reason)");
            return void();
        }
        
        // Rcpp::Rcout<< "\n ==>Retornem";
        
        return void();
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
    extern inline void RcppbdSVD_hdf5( std::string filename, std::string strsubgroup, std::string strdataset,  
                                int k, int q, int nev, bool bcenter, bool bscale, double dthreshold, 
                                bool bforce, bool asRowMajor, 
                                Rcpp::Nullable<Rcpp::CharacterVector> method = R_NilValue,
                                Rcpp::Nullable<int> ithreads = R_NilValue)
    {
        
        try {
            
            hdf5Dataset* dsu;
            hdf5Dataset* dsv;
            hdf5Dataset* dsd;
            std::string strMethod;
            
            std::vector<hsize_t> stride = {1, 1},
                                 block = {1, 1},
                                 offset = {0, 0},
                                 count = {0, 0};
            
            std::vector<std::string> strMethods = {"auto", "blocks", "full"};
            
            if(method.isNull())  strMethod = "auto" ;
            else    strMethod = Rcpp::as<std::string>(method);
            
            if (std::find(strMethods.begin(), strMethods.end(), strMethod) != strMethods.end())
            {
                std::string strmessage = "SVD decomposition by: " + strMethod + " method";
                Rcpp::message(Rcpp::wrap(strmessage));
            } else {
                std::string strWarning = "Method " + strMethod + " not allowed, computing SVD decomposition by 'auto' method";
                strMethod = "auto";
                Rcpp::warning(strWarning);
            }
            
            hdf5Dataset* dsA = new hdf5Dataset(filename, strsubgroup, strdataset, false);
            dsA->openDataset();
            
            // Create results folder
            std::string stroutgroup = "SVD/"+ strdataset;
            
            std::vector<hsize_t> dims_out = {dsA->nrows(), dsA->ncols()};;
            count = { dims_out[0], dims_out[1]};
            
            // Small matrices ==> Direct SVD (lapack)
            if( (dims_out[0] * dims_out[1] < (MAXELEMSINBLOCK / 20) && strMethod == "auto") || strMethod == "full" ) {
                
                Eigen::MatrixXd X;
                svdeig retsvd;
                
                std::vector<double> vdA( count[0] * count[1] ); 
                dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
                
                if(asRowMajor == true) {
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdA.data(), count[0], count[1] );    
                    retsvd = RcppbdSVD_lapack(X, bcenter, bscale, false);
                } else {
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> X (vdA.data(), count[1], count[0] );
                    retsvd = RcppbdSVD_lapack(X, bcenter, bscale, false);
                }
                
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
                
                dsu = new hdf5Dataset(filename, stroutgroup, "u", true);
                dsv = new hdf5Dataset(filename, stroutgroup, "v", true);
                dsd = new hdf5Dataset(filename, stroutgroup, "d", true);
                
                RcppbdSVD_hdf5_Block( dsA, dsu, dsv, dsd, k, q, nev, bcenter, bscale, count[1], count[0], dthreshold, ithreads );
                
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

