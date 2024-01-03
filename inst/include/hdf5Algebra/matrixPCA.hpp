#ifndef BIGDATASTATMETH_HDF5_MATRIXPCA_HPP
#define BIGDATASTATMETH_HDF5_MATRIXPCA_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

#include "Utilities/openme-utils.hpp"
#include "Utilities/Utilities.hpp"
#include "memAlgebra/memOtherFunctions.hpp"
#include "hdf5Algebra/matrixSvd.hpp"
#include "hdf5Algebra/vectormatrix.hpp"
#include "hdf5Utilities/hdf5Utilities.hpp"


namespace BigDataStatMeth {


// Get's variance and cumulative variance from svd decomposition
// var.contr, C, var.coord and var.cos^2 and write results to hdf5 file

// void get_HDF5_PCA_variables_ptr(  H5File* file, std::string strdataset)
void RcppGetPCAVariablesHdf5( std::string strPCAgroup, 
                              BigDataStatMeth::hdf5Dataset* dsd, 
                              BigDataStatMeth::hdf5Dataset* dsv, 
                              bool overwrite )
{
    
    try
    {
        H5::Exception::dontPrint();
        
        int ielements = 0;
        std::vector<hsize_t> stride = {1, 1},
                             block = {1, 1},
                             offset_d = {0, 0},
                             count_d = {dsd->nrows(), dsd->ncols()},
                             offset_v = {0, 0},
                             count_v = {dsv->nrows(), dsv->ncols()};
        
        
            std::vector<double> vdd( count_d[0] * count_d[1] );
            dsd->readDatasetBlock( {offset_d[0], offset_d[1]}, {count_d[0], count_d[1]}, stride, block, vdd.data() );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> d (vdd.data(), count_d[0], count_d[1]);
            
            // Rcpp::Rcout<<"\n d: "<<d.rows()<<" x "<<d.cols()<<"\n";
            
        {    
            // lambda
            Eigen::VectorXd vvar = d.array().pow(2);
            BigDataStatMeth::hdf5Dataset* dslambda = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "lambda" ,overwrite );
            dslambda->createDataset( d.rows(), d.cols(), "real");
            dslambda->writeDataset( vvar.data() );
            delete dslambda;
            
            // Variance
            vvar = (vvar/vvar.array().sum()) * 100;
            BigDataStatMeth::hdf5Dataset* dsvar = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "variance" ,overwrite );
            dsvar->createDataset( d.rows(), d.cols(), "real");
            dsvar->writeDataset( vvar.data() );
            delete dsvar;
            
            // cumulative variance (max 1000 elements (firsts))
            if(vvar.size()>1000){  ielements = 1000;    }
            else{ ielements = vvar.size();    }
            
            BigDataStatMeth::hdf5Dataset* dscumvar = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "cumvar" ,overwrite );
            dscumvar->createDataset( d.rows(), d.cols(), "real");
            dscumvar->writeDataset( (cumsum(vvar)).data());
            delete dscumvar;
        }
        
        {
            Eigen::MatrixXd var_coord;
            {
                // Coord vars
                std::vector<double> vdv( count_v[0] * count_v[1] );
                dsv->readDatasetBlock( {offset_v[0], offset_v[1]}, {count_v[0], count_v[1]}, stride, block, vdv.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> v (vdv.data(), count_v[0], count_v[1]);
                
                var_coord = Rcpp_matrixVectorMultiplication_byRow(v, d);
            }
            
            BigDataStatMeth::hdf5Dataset* dscoord = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "var.coord" ,overwrite );
            dscoord->createDataset( var_coord.rows(), var_coord.cols(), "real");
            dscoord->writeDataset( var_coord.transpose().data() );
            delete dscoord;

            // Cos2
            Eigen::MatrixXd var_cos2 = var_coord.unaryExpr([](double d) {return std::pow(d, 2);});
            BigDataStatMeth::hdf5Dataset* dscos2 = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "var.cos2" ,overwrite );
            dscos2->createDataset( var_cos2.rows(), var_cos2.cols(), "real");
            dscos2->writeDataset( var_cos2.transpose().data() );
            delete dscos2;
             
        }
        
    }catch( H5::FileIException& error ) {
        ::Rf_error( "c++ RcppGetPCAVariablesHdf5 exception (File IException )" );
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ RcppGetPCAVariablesHdf5 exception (DataSet IException )" );
        return void();
    } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ RcppGetPCAVariablesHdf5 exception (DataSpace IException )" );
        return void();
    } 
    
    return void();
}




void RcppGetPCAIndividualsHdf5( std::string strPCAgroup, 
                                BigDataStatMeth::hdf5Dataset* dsX,
                                BigDataStatMeth::hdf5Dataset* dsd, 
                                BigDataStatMeth::hdf5Dataset* dsu, 
                                bool overwrite )
{
    
    try
    {
        H5::Exception::dontPrint();
        
        int ielements = 0;
        std::vector<hsize_t> stride = {1, 1},
            block = {1, 1},
            offset = {0, 0},
            count_x = {dsX->nrows(), dsX->ncols()},
            count_d = {dsd->nrows(), dsd->ncols()},
            count_u = {dsu->nrows(), dsu->ncols()};
        
        Eigen::VectorXd dist2;
        Eigen::MatrixXd adjust;
        double weights = 1/(double)dsX->ncols();
        {
            
            std::vector<double> vdX( count_x[0] * count_x[1] );
            dsX->readDatasetBlock( {0,0}, count_x, stride, block, vdX.data() );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdX.data(), count_x[0], count_x[1]);

            X = RcppNormalizeRowwise(X, true, false);
            
            adjust = Eigen::MatrixXd::Constant(1, X.cols(), weights);
            
            Eigen::RowVectorXd ecart =  adjust / adjust.sum() ;
            ecart = ecart * (X.unaryExpr([](double d) {return std::pow(d, 2);})).transpose();
            ecart = ecart.unaryExpr([](double d) {return std::sqrt(d);});
            
            X = X.array().colwise() / ecart.transpose().array();

            dist2 =  (X.unaryExpr([](double d) {return std::pow(d, 2);})).colwise().sum();
            Eigen::VectorXd dist = dist2.array().sqrt();
            
            BigDataStatMeth::hdf5Dataset* dsdist2 = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "ind.dist" ,overwrite );
            dsdist2->createDataset( dist2.rows(), dist2.cols(), "real");
            dsdist2->writeDataset( dist.data() );
            delete dsdist2;
        }

        // Load d
        std::vector<double> vdd( count_d[0] * count_d[1] );
        dsd->readDatasetBlock( offset, count_d, stride, block, vdd.data() );
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> d (vdd.data(), count_d[0], count_d[1]); 
        
        Eigen::MatrixXd var_coord;
        {
            // Load u
            std::vector<double> vdu( count_x[0] * count_x[1] );
            dsu->readDatasetBlock( {0, 0}, count_x, stride, block, vdu.data() );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> u (vdu.data(), count_x[0], count_x[1]);
            
            u = u * sqrt(1/weights); 
        
            // Components
            BigDataStatMeth::hdf5Dataset* dsComp = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "components" ,overwrite );
            dsComp->createDataset( u.cols(), u.rows(), "real");
            dsComp->writeDataset( u.data() );
            delete dsComp;    
            
            var_coord = Rcpp_matrixVectorMultiplication_byRow(u, d);
        }
        
        // Coord inds
        BigDataStatMeth::hdf5Dataset* dscoord = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "ind.coord" ,overwrite );
        dscoord->createDataset( var_coord.cols(), var_coord.rows(), "real");
        dscoord->writeDataset( Rcpp::wrap(var_coord.transpose()));
        delete dscoord;
        
        // Cos2 inds
        Eigen::MatrixXd coord2 = var_coord.unaryExpr([](double d) {return std::pow(d, 2);});
        
        {
            Eigen::MatrixXd ind_cos2 = coord2.array().rowwise() / dist2.transpose().array();
            
            BigDataStatMeth::hdf5Dataset* dscos2 = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "ind.cos2" ,overwrite );
            dscos2->createDataset( ind_cos2.cols(), ind_cos2.rows(), "real");
            dscos2->writeDataset( Rcpp::wrap(ind_cos2.transpose()) );
            delete dscos2;
        }
        
        Eigen::MatrixXd ind_contrib = coord2.array().rowwise() * adjust.row(0).array();
        
        d = d.unaryExpr([](double xd) {return std::pow(xd, 2);});
        
        ind_contrib = ind_contrib.array().colwise() / d.col(0).array();
        ind_contrib = ind_contrib.unaryExpr([](double d) {return d*100;});
        
        BigDataStatMeth::hdf5Dataset* dscontrib = new BigDataStatMeth::hdf5Dataset(dsd->getFileName(), strPCAgroup, "ind.contrib" ,overwrite );
        dscontrib->createDataset( ind_contrib.cols(), ind_contrib.rows(), "real");
        dscontrib->writeDataset( Rcpp::wrap(ind_contrib.transpose()) );
        delete dscontrib;
            
        // Components
        // createHardLink(dsu->getFileptr(), dsu->getGroupName() + "/" + dsu->getDatasetName(), strPCAgroup + "/components");
        
    }catch( H5::FileIException& error ) {
        ::Rf_error( "c++ RcppGetPCAIndividualsHdf5 exception (File IException )" );
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ RcppGetPCAIndividualsHdf5 exception (DataSet IException )" );
        return void();
    } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ RcppGetPCAIndividualsHdf5 exception (DataSpace IException )" );
        return void();
    } 
    
    return void();
}





extern void RcppPCAHdf5( std::string filename, std::string strgroup, std::string strdataset,  
                         std::string strSVDgroup, int k, int q, int nev, 
                         bool bcenter, bool bscale, double dthreshold, 
                         bool bforce, bool asRowMajor, Rcpp::Nullable<int> ithreads = R_NilValue)
{
    
    try{
        
        // 
        //  Que he de fer?
        //      1. Mirar si el SVD està calculat (Falta testejar existencia datasets)
        //      2. Depenent SVD
        //          2.1 No està calculat-> Calcular-lo
        //          2.2 Està calculat -> (3)
        //      3. Calcular resta PCA amb resultats SVD
        //
        
        
        std::string strPCAgroup = "PCA/" + strdataset;
        
        // Check for svd decomposition (u, v and d matrices) in hdf5 file or if we 
        // need to compute again the SVD ( foce = true )
        BigDataStatMeth::hdf5File* file = new BigDataStatMeth::hdf5File(filename, false);
        file->openFile("r");
        
            bool bexistsSVD = exists_HDF5_element(file->getFileptr(), strSVDgroup);
            bool bexistsPCA = exists_HDF5_element(file->getFileptr(), strPCAgroup);
        
        delete file;
        
        if( !bexistsSVD ||  bforce == true ) {
            
            BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, strgroup, strdataset, false);
            dsA->openDataset();
            RcppTypifyNormalizeHdf5( dsA, bcenter, bscale, false); // Normalize and tipify data ( ((x-mu)/(sd)) * 1/sqrt(n-1) )
            delete dsA;
            
            BigDataStatMeth::RcppbdSVD_hdf5( filename, "NORMALIZED_T/" + strgroup, strdataset, k, q, nev, false, false, dthreshold, bforce, asRowMajor, ithreads );
            
            strSVDgroup = "SVD/" +  strdataset;
            
        }
        
        // Check if PCA decomposition exists
        if( bexistsPCA && bforce == false) {
            Rcpp::Rcout<<"PCA decomposition exits, please set overwrite = true to overwrite the existing results";
            return void();
        }
        
        
        
        BigDataStatMeth::hdf5Dataset* dsd = new BigDataStatMeth::hdf5Dataset(filename, strSVDgroup, "d", false );
        dsd->openDataset();
        
        // ------------ Variables ----------------
        
        BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, "NORMALIZED_T/" + strgroup, strdataset, false );
        dsA->openDataset();
        
        BigDataStatMeth::hdf5Dataset* dsv = new BigDataStatMeth::hdf5Dataset(filename, strSVDgroup, "v", false );
        dsv->openDataset();
        
        RcppGetPCAVariablesHdf5( strPCAgroup, dsd, dsv, bforce );
        
        delete dsv;
        delete dsA;
        
        // ------------ Individuals ----------------
        
        BigDataStatMeth::hdf5Dataset* dsX = new BigDataStatMeth::hdf5Dataset(filename, strgroup, strdataset, false);
        dsX->openDataset();
        
        BigDataStatMeth::hdf5Dataset* dsu = new BigDataStatMeth::hdf5Dataset(filename, strSVDgroup, "u", false );
        dsu->openDataset();
        
        Rcpp::Rcout<<"\nnComputing PCA";
        RcppGetPCAIndividualsHdf5( strPCAgroup, dsX, dsd, dsu, bforce );
        Rcpp::Rcout<<"\nPCA Computed";
        
        delete dsd;
        delete dsu;
        delete dsX;
        
    }catch( H5::FileIException& error ) {
        ::Rf_error( "c++ exception RcppPCAHdf5 (File IException)" );
        return void();
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception RcppPCAHdf5 (DataSet IException)" );
        return void();
    } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception RcppPCAHdf5 (DataSpace IException)" );
        return void();
    } catch( H5::DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception RcppPCAHdf5 (DataType IException)" );
        return void();
    } catch(std::exception &ex) {
        Rcpp::Rcout<< "C++ exception RcppPCAHdf5 : "<< ex.what();
    } catch (...) {
        ::Rf_error("C++ exception RcppPCAHdf5 (unknown reason)");
        return void();
    }
    
    return void();
    
}




}

#endif // BIGDATASTATMETH_HDF5_MATRIXPCA_HPP