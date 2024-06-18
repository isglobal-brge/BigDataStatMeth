#ifndef BIGDATASTATMETH_HDF5_MATRIXSVDBLOCK_HPP
#define BIGDATASTATMETH_HDF5_MATRIXSVDBLOCK_HPP

// #include <RcppEigen.h>
#include "hdf5Utilities/hdf5Utilities.hpp"
#include "hdf5Algebra/matrixSdMean.hpp"
#include "hdf5Algebra/matrixNormalization.hpp"
// #include "hdf5Algebra/matrixSvd.hpp"
#include "memAlgebra/memMultiplication.hpp"
// #include "H5Cpp.h"


namespace BigDataStatMeth {

// dgesvd_ is a symbol in the LAPACK-BLAS Level 3 
//    DGESVD computes the singular value decomposition (SVD) of a real M-by-N matrix A, 
//       optionally computing the left and/or right singular vectors
extern "C" {
    extern void dgesvd_( char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
}

extern "C" {
    extern void dgesdd_( char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*, int*);
}



// Lapack SVD decomposition - Optimized algorithm with dgesdd
template <class T>
extern inline svdeig RcppbdSVD_lapack( T X, bool bcenter, bool bscale, bool complete ) {
    
    svdeig retsvd;
    
    char Schar='S';
    char Achar='A';
    int info = 0;
    
    try {
        
        if(bcenter == true || bscale == true) {
            X = RcppNormalize_Data(X, bcenter, bscale, false);
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




// Reads big matrix from hdf5 file in blocks and perform a svd descomposition from
// each block,results are saved in hdf5 datasets under temporal group to be processed
// if necessary
// int First_level_SvdBlock_decomposition_hdf5(H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
//                                             int irows, int icols, Rcpp::Nullable<int> threads = R_NilValue)
template <class T>
extern inline void First_level_SvdBlock_decomposition_hdf5( T* dsA, std::string strGroupName, int k, int q, int nev, bool bcenter, bool bscale, 
                                             int irows, int icols, double dthreshold, Rcpp::Nullable<int> threads = R_NilValue)
{
    
    static_assert(std::is_same<T*, BigDataStatMeth::hdf5Dataset* >::value || 
                  std::is_same<T*, BigDataStatMeth::hdf5DatasetInternal* >::value,
                  "Error - type not allowed");
    
    BigDataStatMeth::hdf5Dataset* normalizedData;
    BigDataStatMeth::hdf5Dataset* unlimDataset;
    std::vector<hsize_t> stride = {1, 1},
                         block = {1, 1};
    int  M, p, n;
    int maxsizetoread, 
        ithreads, 
        cummoffset;
    bool transp = false;
    
    Rcpp::Nullable<int> wsize = R_NilValue;
    
    try{
        
        irows = dsA->ncols();
        icols = dsA->nrows();
        
        ithreads = get_number_threads(threads, R_NilValue);
        
        // Work with transposed matrix
        if( irows >= icols ) {
            transp = true;
            n = icols;
            p = irows;
        } else {
            n = irows;
            p = icols;
        }
        
        Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2, n);
        
        //..2023/10/02..// datanormal = Eigen::MatrixXd::Zero( 2, dsA->nrows());
        if(!transp) {
            get_HDF5_mean_sd_by_row( dsA, datanormal, wsize);
        } else {
            get_HDF5_mean_sd_by_column(dsA, datanormal, wsize);    
        }
        
        
        // get_HDF5_mean_sd_by_column_partial_ptr( file, dataset, datanormal);
        
        M = pow(k, q);
        if(M>p)
            throw std::runtime_error("k^q must not be greater than the number of columns in the matrix");
        
        //.. 2023/10/04 ..//double block_size = std::floor((double)p/(double)M); 
        double block_size = std::floor((double)icols/(double)M); 
        int realsizeread;
        
        if (bcenter==true || bscale==true) {
            // Create dataset to store normalized data
            normalizedData = new BigDataStatMeth::hdf5DatasetInternal (dsA->getFileName(), strGroupName, "normalmatrix", true);
            //..2023/10/02 ..// normalizedData->createDataset( n, p, "real");
            normalizedData->createDataset( irows, icols, "real");

        }
        
        #pragma omp parallel num_threads(ithreads)
        {
            
             // int tid = omp_get_thread_num();
             // if( tid == 0 )   {
             //     Rcpp::Rcout << "Number of threads C++ : " << omp_get_num_threads() << "\n";
             //     Rcpp::Rcout << "Thread number C++ : " << omp_get_thread_num() << "\n";
             // }
            
            int chunks = M/ithreads;
            
            // Get data from M blocks in initial matrix
            //.. 2024/04/01 ..// #pragma omp for ordered schedule (dynamic, chunks) 
            #pragma omp for schedule (dynamic) ordered
            for( int i = 0; i< M ; i++)  
            {
                
                Eigen::MatrixXd restmp;
                std::string strDatasetName = strGroupName + "/A" + std::to_string(i/(M/k));
                
                // 1.- Get SVD from all blocks
                //    a) Get all blocks from initial matrix 
                std::vector<hsize_t>  offset = getInitialPosition( transp, (unsigned long long)(i*block_size) ); // Initial read position
                std::vector<hsize_t> count = {0, 0};
                
                maxsizetoread = block_size;
                
                // Get max block size to read - for blocks smaller than default block size 
                
                if( ((i+1)*block_size) > icols)
                    maxsizetoread = icols - (i*block_size);
                
                if( i+1 == M && icols - maxsizetoread!=0)
                    realsizeread = icols - (i*block_size);
                else
                    realsizeread = maxsizetoread;
                
                Eigen::MatrixXd X;
                
                count = getSizetoRead(transp, (unsigned long long)(realsizeread), icols, irows );
                
                if(transp==false){
                    std::vector<double> vdX( count[0] * count[1] ); 
                    #pragma omp critical(accessFile)
                    {
                        dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdX.data() );
                    }
                    X = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> (vdX.data(), count[1], count[0] );
                    
                } else {
                    count[0] = irows; // Mirar de canviar i posar d'integrar dins de la funció getSizetoRead()
                    std::vector<double> vdX( count[0] * count[1] );
                    #pragma omp critical(accessFile)
                    {
                        dsA->readDatasetBlock( {offset[1], offset[0]}, {count[1], count[0]}, stride, block, vdX.data() );
                    }
                    X = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> (vdX.data(), count[0], count[1] );
                    
                }
        
        // Rcpp::Rcout<<"\nDebug - 1";
                // Normalize data
                if (bcenter==true || bscale==true) 
                {
                    std::vector<hsize_t> offset_tmp = {offset[1],offset[0]};
                    std::vector<hsize_t> count_tmp = {count[1], count[0]};
                    
                    if(transp == true) {
                        offset_tmp = {offset[0], offset[1]};
                        count_tmp = {count[0], count[1]};
                    } 
                        
                    if( transp == false) {
                        X = RcppNormalize_Data(X, bcenter, bscale, transp, datanormal.block(0, offset_tmp[1], 2, count_tmp[1]));
                        #pragma omp critical(accessFile)
                        {   
                            normalizedData->writeDatasetBlock( Rcpp::wrap(X), offset_tmp, count_tmp, stride, block, false);
                        }
                    } else {
                        X = RcppNormalize_Data(X, bcenter, bscale, transp, datanormal.block(0, offset_tmp[1], 2, count_tmp[1]));
                        count_tmp = {(unsigned long long)X.rows(), (unsigned long long)X.cols()};
                        #pragma omp critical(accessFile)
                        {   
                            normalizedData->writeDatasetBlock( Rcpp::wrap(X), offset_tmp, count_tmp, stride, block, false);
                        }
                    }
                
                }
                // Rcpp::Rcout<<"\nDebug - 2";
        
                {
                    //    b) SVD for each block
                    svdeig retsvd;
                    retsvd = RcppbdSVD_lapack(X, false, false, false);
                    
                    int size_d = (retsvd.d).size();
                    int nzeros = 0;
                    
                    if( (retsvd.d)[size_d - 1] <= dthreshold ){
                        nzeros = 1;
                        for( int j = (size_d - 2); ( j>1 && (retsvd.d)[i] <= dthreshold ); j-- ) {
                            nzeros++;
                        }
                    }
                    // Rcpp::Rcout<<"\nDebug - 2.1";
                    
                    //    c)  U*d
                    // Create diagonal matrix from svd decomposition d
                    int isize = (retsvd.d).size() - nzeros;
                    
                    if( isize < 2 ) {
                        isize = 2;
                    }
                    
                    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
                    d.diagonal() = (retsvd.d).head(isize);
                    
                    // Rcpp::Rcout<<"\nDebug - 2.2";
                    Eigen::MatrixXd u = (retsvd.u).block(0, 0, (retsvd.u).rows(), isize);
                    //..2024/03/27 ..// restmp = block_matrix_mul( u, d, 1024);
                    
                    // 
                    //          !!!!!!!!!!   REVISAR-HO !!!!!!!!!!!!
                    // 
                    //  Quina operació fa en realitat el u*d?? és una multiplicació
                    //  de matrius o és una multiplicació matriu * vector?
                    //  
                    //      MOLT IMPORTANT!!! Confirmar resultats !!!!
                    // 
                    if( u.size() < MAXELEMSINBLOCK ) {
                        restmp = u*d;
                    } else{
                        restmp = Rcpp_block_matrix_mul( u, d, R_NilValue);
                    }
                    
                }
                
                // Rcpp::Rcout<<"\nDebug - 3";
                
                //    d) Write results to hdf5 file
                offset[0] = 0; offset[1] = 0;
                count[0] = restmp.rows();
                count[1] = restmp.cols();
        
                #pragma omp ordered
                {
                    // #pragma omp critical(accessFile)
                    // {
                    if( i%(M/k) == 0 || ( (i%(M/k) > 0 &&  !exists_HDF5_element(dsA->getFileptr(),  strDatasetName)) ) ) {
                        // Create unlimited dataset in hdf5 file
                        unlimDataset = new BigDataStatMeth::hdf5DatasetInternal(dsA->getFileName(), strDatasetName, true );
                        #pragma omp critical(accessFile)
                        {
                            unlimDataset->createUnlimitedDataset(count[0], count[1], "real");
                        }
                        delete unlimDataset;
                        cummoffset = 0;
                    }
                    // Rcpp::Rcout<<"\nDebug - 3.1";
                    // Get write position
                    offset[1] = cummoffset;
                    cummoffset = cummoffset + restmp.cols();
                    
                    #pragma omp critical(accessFile)
                    {
                        unlimDataset = new BigDataStatMeth::hdf5DatasetInternal(dsA->getFileName(), strDatasetName, true );
                        unlimDataset->openDataset();
                        
                        // Extend dataset before put data
                        if((i%(M/k)) != 0 && cummoffset > 0) {
                                unlimDataset->extendUnlimitedDataset(0, count[1] );
                        }
                        
                        unlimDataset->writeDatasetBlock( Rcpp::wrap(restmp), offset, count, stride, block, false);
                    }
                    delete unlimDataset;
                    // Rcpp::Rcout<<"\nDebug - 3.2";
                    
                }
            }
        }
        
        // Rcpp::Rcout<< "\n ===  FINALITZEM EXECUCIÓ !!  ===\n";
        // Rcpp::Rcout<<"\nDebug - 4";
        if (bcenter==true || bscale==true) { 
            Rcpp::Rcout<<"\nSuposo que no se m'ha acudit entrar al delete del normalizeData... \n";
            delete normalizedData; 
        }
        // Rcpp::Rcout<<"\nDebug - 5";
        
    } catch(std::exception& ex) {
        Rcpp::Rcout<< "c++ exception First_level_SvdBlock_decomposition_hdf5: "<<ex.what()<< " \n";
        return void();
    }
    
    // Rcpp::Rcout<< "\n ======>  POSEM PUNT I FINAL!!!!!! !!  <======\n";
    
    return void();
    
}



// Reads small datasets from hdf5 and perform a svd descomposition from each block,
// results are saved in hdf5 datasets under temporal group to be processed if necessary
// int Next_level_SvdBlock_decomposition_hdf5(H5File* file, std::string strGroupName, int k, int q,
//                                            bool bcenter, bool bscale, Rcpp::Nullable<int> threads = R_NilValue)
template <class T>
extern inline void Next_level_SvdBlock_decomposition_hdf5( T* dsA, std::string strGroupName, int k, int q, 
                                                           double dthreshold, Rcpp::Nullable<int> threads = R_NilValue)
// extern inline void Next_level_SvdBlock_decomposition_hdf5( T* dsA, std::string strGroupName, int k, int q, bool bcenter, bool bscale, 
//                                                     int irows, int icols, double dthreshold, Rcpp::Nullable<int> threads = R_NilValue)
    {
    
    
    static_assert(std::is_same<T*, BigDataStatMeth::hdf5Dataset* >::value || 
                  std::is_same<T*, BigDataStatMeth::hdf5DatasetInternal* >::value,
                  "Error - type not allowed");

    BigDataStatMeth::hdf5Dataset* unlimDataset;
    int cummoffset = 0; //.. MODIFICAT 07/09/2020 ..//
    int ithreads;
    std::vector<hsize_t> stride = {1, 1},
                         block = {1, 1},
                         offset = {0, 0},
                         count = {0, 0};

    int M;
    
    /*** TODO:
     *
     * Revisar el cas que no existeixi el dataset per a fer join perquè aquell dataset no s'ha creat perquè tots els
     * valors singulars << threshold !!!
     *
     * ***/
    
    
    Rcpp::CharacterVector strvmatnames = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};

    try {

        // Get dataset names
        Rcpp::StringVector joindata =  dsA->getDatasetNames(strGroupName, (std::string)strvmatnames[q-1]);

        M = joindata.size();
        
        ithreads = get_number_threads(threads, R_NilValue);

        
        #pragma omp parallel num_threads(ithreads)

        // Get data from M blocks in initial matrix
        #pragma omp for ordered schedule (dynamic)
        for( int i = 0; i< M ; i++)
        {
            std::string strDatasetName = strGroupName + "/" + strvmatnames[q] + std::to_string(i/(M/k));
            
            Eigen::MatrixXd restmp;
            Eigen::MatrixXd X;

            #pragma omp critical(accessFile)
            {
                //    a) Get dataset
                BigDataStatMeth::hdf5Dataset* dsCur = new BigDataStatMeth::hdf5Dataset(dsA->getFileName(), strGroupName + "/" + joindata[i], false);
                dsCur->openDataset();
                hsize_t* dims_out = dsCur->dim();
                
                std::vector<double> vdCurDataset( dims_out[0] * dims_out[1] ); 
                dsCur->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, stride, block, vdCurDataset.data() );
                // Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> X (vdCurDataset.data(), dims_out[0], dims_out[1] );
                X = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdCurDataset.data(), dims_out[0], dims_out[1] );
                
                delete dsCur;
            }

            {
                //    b) SVD for each block
                svdeig retsvd;
                retsvd = RcppbdSVD_lapack(X, false, false, false);
            
                int size_d = (retsvd.d).size();
                int nzeros = 0;
            
                if( (retsvd.d)[size_d - 1] <= dthreshold ) {
                    nzeros = 1;
                    for( int j = (size_d - 2); ( j>1 && (retsvd.d)[i] <= dthreshold ); j-- ) {
                        nzeros++;
                    }
                }
            
                //    c)  U*d
                // Create diagonal matrix from svd decomposition d
                int isize = (retsvd.d).size() - nzeros;
                if( isize < 2 ) {
                    isize = 2;
                }
            
                Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
                d.diagonal() = (retsvd.d).head(isize);
            
                //..ORIGINAL..//restmp = Bblock_matrix_mul((retsvd.u).block(0, 0, (retsvd.u).rows(), isize), d, 1024);
                Eigen::MatrixXd u = (retsvd.u).block(0, 0, (retsvd.u).rows(), isize);
                //..2024/03/27 ..// restmp = block_matrix_mul( u, d, 1024);
            
                if( u.size() < MAXELEMSINBLOCK ) {
                    restmp = u*d;
                } else{
                    restmp = Rcpp_block_matrix_mul( u, d, R_NilValue);    
                } 
                
            }
            
            //    d) Write results to dataset
            count[0] = restmp.rows();
            count[1] = restmp.cols();

            #pragma omp ordered
            {
                #pragma omp critical(accessFile)
                {
                    
                    if( i%(M/k) == 0 || ( (i%(M/k) > 0 &&  !BigDataStatMeth::exists_HDF5_element(dsA->getFileptr(),  strDatasetName)) ) ) {
                        // Create unlimited dataset in hdf5 file
                        unlimDataset = new BigDataStatMeth::hdf5DatasetInternal(dsA->getFileName(), strDatasetName, true );
                        unlimDataset->createUnlimitedDataset(count[0], count[1], "real");
                        delete unlimDataset;
                        
                        cummoffset = 0;
                    }
                    
                    unlimDataset = new BigDataStatMeth::hdf5DatasetInternal(dsA->getFileName(), strDatasetName, true );
                    unlimDataset->openDataset();
                    
                    // Get write position
                    offset[1] = cummoffset;
                    cummoffset = cummoffset + restmp.cols();
                    
                    // Extend dataset before put data
                    if((i%(M/k)) != 0 && cummoffset > 0) {
                        unlimDataset->extendUnlimitedDataset(0, count[1] );
                    }
                    unlimDataset->writeDatasetBlock( Rcpp::wrap(restmp), offset, count, stride, block, false);
                    delete unlimDataset;
                }
            }
        }

        //..ONLY DEBUG !!!...//remove_HDF5_multiple_elements_ptr(file, strGroupName, joindata);

    } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (File IException)" );
        return void();
    } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
        ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (DataSet IException)" );
        return void();
    } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
        ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (Group IException)" );
        return void();
    } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (DataSpace IException)" );
        return void();
    } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
        ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (Data TypeIException)" );
        return void();
    }catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return void();
    }

    return void();

}

}

#endif // BIGDATASTATMETH_HDF5_MATRIXSVDBLOCK_HPP