#ifndef BIGDATASTATMETH_HDF5_MATRIXSVDBLOCK_HPP
#define BIGDATASTATMETH_HDF5_MATRIXSVDBLOCK_HPP

#include <RcppEigen.h>
#include "hdf5Utilities/hdf5Utilities.hpp"
// #include "hdf5Algebra/matrixNormalization.hpp"
// #include "H5Cpp.h"


namespace BigDataStatMeth {


// Reads big matrix from hdf5 file in blocks and perform a svd descomposition from
// each block,results are saved in hdf5 datasets under temporal group to be processed
// if necessary
// int First_level_SvdBlock_decomposition_hdf5(H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
//                                             int irows, int icols, Rcpp::Nullable<int> threads = R_NilValue)
extern void First_level_SvdBlock_decomposition_hdf5( BigDataStatMeth::hdf5Dataset* dsA, int k, int q, int nev, bool bcenter, bool bscale, 
                                             int irows, int icols, double dthreshold, Rcpp::Nullable<int> threads = R_NilValue)
{
    
    BigDataStatMeth::hdf5Dataset* normalizedData;
    
    Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1, 1);
    Rcpp::IntegerVector block = Rcpp::IntegerVector::create(1, 1);
    
    
    int  M, p, n;
    int maxsizetoread;
    bool transp = false;
    std::string strGroupName  = "tmpgroup";
    Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2,icols);
    int cummoffset, ret;
    // DataSet* normalizedData;
    int ithreads;
    Rcpp::Nullable<int> wsize  = R_NilValue;
    
    
    // if( exists_HDF5_element_ptr(file,strGroupName))
    //     remove_HDF5_element_ptr(file,strGroupName);
    // 
    // ret = create_HDF5_group_ptr(file, strGroupName);
    
    // // Create group
    // BigDataStatMeth::hdf5Group* tmpGroup(dsA->getFileptr(), strGroupName, true);
    
    
    
    try{
        
        irows = dsA->nrows();
        icols = dsA->ncols();
        
        if(threads.isNotNull()) 
        {
            if (Rcpp::as<int>(threads) < std::thread::hardware_concurrency())
                ithreads = Rcpp::as<int> (threads);
            else 
                ithreads = std::thread::hardware_concurrency()/2;
        }
        else    ithreads = std::thread::hardware_concurrency()/2; 
        
        
        // Work with transposed matrix
        if( irows >= icols ) {
            transp = true;
            n = icols;
            p = irows;
            
        } else {
            n = irows;
            p = icols;
        }
        

        datanormal = Eigen::MatrixXd::Zero( 2, dsA->ncols());
        get_HDF5_mean_sd_by_row( dsA, datanormal, wsize);
        
        // get_HDF5_mean_sd_by_column_partial_ptr( file, dataset, datanormal);
        
        M = pow(k, q);
        if(M>p)
            throw std::runtime_error("k^q must not be greater than the number of columns in the matrix");
        
        double block_size = std::floor((double)p/(double)M); 
        int realsizeread;
        
        
        if (bcenter==true || bscale==true) {
            // Create dataset to store normalized data
            normalizedData = new BigDataStatMeth::hdf5Dataset(dsA->getFileName(), dsA->getGrouName(), "normalmatrix", bforce);
            normalizedData->createDataset( n, p, "real");

        }
        
        
        //
        // !!!! ESTIC AQUÍ !!! 
        //      EL QUE EM TOCA FER ARA ÉS MIRAR EL CODI AL COMPLET,
        //      QUINES MATRIUS HE DE CREAR, PER A QUÈ LES HE DE CREAR
        //      I QUE HE DE FER AMB CADA UNA D'ELLES PER TAL
        //      D'ESTALVIAR PASSOS MEMÒRIA I MALS DE CAP.
        //
        
#pragma omp parallel num_threads(getDTthreads(ithreads, false))
        
        /*
         int tid = omp_get_thread_num();
         if( tid == 0 )   {
         Rcpp::Rcout << "Number of threads C++ : " << omp_get_num_threads() << "\n";
         Rcpp::Rcout << "Thread number C++ : " << omp_get_thread_num() << "\n";
         }
         */
        
        
        // Get data from M blocks in initial matrix
#pragma omp for ordered schedule (static) 
        for( int i = 0; i< M ; i++)  
        {
            
            Eigen::MatrixXd restmp;
            std::string strDatasetName = strGroupName + "/A" + std::to_string(i/(M/k));
            
            // 1.- Get SVD from all blocks
            //    a) Get all blocks from initial matrix 
            
            Rcpp::IntegerVector offset = getInitialPosition( transp, (unsigned long long)(i*block_size) ); // Initial read position
            Rcpp::IntegerVector count; // Blocks size
            maxsizetoread = block_size;
            
            // Get max block size to read - for blocks smaller than default block size 
            if(transp == true)
            {
                if( ((i+1)*block_size) > irows)
                    maxsizetoread = irows - (i*block_size);
                
                if( i+1 == M && irows - maxsizetoread!=0)
                    realsizeread = irows - (i*block_size);
                else
                    realsizeread = maxsizetoread;
                
            } else {
                
                if( ((i+1)*block_size) > icols)
                    maxsizetoread = icols - (i*block_size);
                
                if( i+1 == M && icols - maxsizetoread!=0)
                    realsizeread = icols - (i*block_size);
                else
                    realsizeread = maxsizetoread;
            }
            
            
            Eigen::MatrixXd X;
         
            
            //
            //  **************************************************************************************************
            //  **************************************************************************************************
            //
            //
            // !!!! ARA ESTIC AQUÍ !!! 
            //      VEIENT QUE REDIMONIS HE DE FER AMB AQUESTES LESCTURES....
            //      TRANSPOSAT, NO TRANSPOSAT, MIDA DEL BLOCK,..... HO HE DE REANALITZAR TOT
            //      DES DEL PRINCIPI PERÒ HO HE DE FER QUAN NO EM FACI MAL EL CAP TAL I COM HO ESTÀ FENT EN 
            //      AQUESTS MOMENTS.... TOT EL QUE HI HA A SOBRE S'HA DE REVISAR PERQUÈ NO ESTÀ TESTEJAT
            //      DE CAP DE LES MANERES... LA FORMA DE FER-HO SERIA FENT EL TEST JUNTAMENT AMB LA VERSIÓ
            //      ANTERIOR DE BDSM I COMPROVAR QUE TOTS ELS RESULTATS SON EXACTAMENT IGUALS !!!!
            //
            //
            //  **************************************************************************************************
            //  **************************************************************************************************
            //
            
#pragma omp critical(accessFile)
{
    
    count = getSizetoRead(transp, (unsigned long long)(realsizeread), icols, irows );
    
    /********** Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, dataset, offset[0], offset[1], count[0], count[1]); **********/
    
    if(transp==false){
        X = GetCurrentBlock_hdf5_Original( file, dataset, offset[0], offset[1], count[0], count[1]);
    } else {
        X = GetCurrentBlock_hdf5( file, dataset, offset[0], offset[1], count[0], count[1]);
    }
    
}

// Normalize data
if (bcenter==true || bscale==true) {
    Rcpp::IntegerVector offset_tmp = Rcpp::IntegerVector::create(offset[0], offset[1]);
    Rcpp::IntegerVector count_tmp = Rcpp::IntegerVector::create(count[0], count[1]);
    
    if(transp == false) {
        offset_tmp[0] = offset[1]; offset_tmp[1] = offset[0];
        count_tmp[0] = count[1]; count_tmp[1] = count[0];
    }
    
#pragma omp critical(accessFile)
{   
    X = RcppNormalize_Data_hdf5(X, bcenter, bscale, transp, datanormal);
    write_HDF5_matrix_subset_v2(file, normalizedData, offset_tmp, count_tmp, stride, block, Rcpp::wrap(X));
}

}

/***
 
 WITHOUT MATRIX REDUCTION --> Works OK
 
 {
 
 //    b) SVD for each block
 svdeig retsvd;
 retsvd = RcppbdSVD_lapack(X, false, false, false);
 
 //    c)  U*d
 // Create diagonal matrix from svd decomposition d
 int isize = (retsvd.d).size();
 Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
 d.diagonal() = retsvd.d;
 
 //..// Eigen::MatrixXd restmp = Bblock_matrix_mul_parallel(retsvd.u, d, 1024,1);
 restmp = Bblock_matrix_mul(retsvd.u, d, 1024);
 
 }
 
 ***/


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
    
    //    c)  U*d
    // Create diagonal matrix from svd decomposition d
    int isize = (retsvd.d).size() - nzeros;
    
    if( isize < 2 ) {
        isize = 2;
    }
    
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
    d.diagonal() = (retsvd.d).head(isize);
    
    //..// Eigen::MatrixXd restmp = Bblock_matrix_mul_parallel(retsvd.u, d, 1024,1);
    restmp = Bblock_matrix_mul((retsvd.u).block(0, 0, (retsvd.u).rows(), isize), d, 1024);
    
}




//    d) Write results to hdf5 file
offset[0] = 0; offset[1] = 0;
count[0] = restmp.rows();
count[1] = restmp.cols();

#pragma omp ordered
{
    
#pragma omp critical(accessFile)
{
    
    if(i%(M/k) == 0 || ( (i%(M/k) > 0 &&  !exists_HDF5_element_ptr(file,strDatasetName)) ) ) {
        
        // If dataset exists --> remove dataset
        if( exists_HDF5_element_ptr(file,strDatasetName)){
            remove_HDF5_element_ptr(file,strDatasetName);
        }
        
        // Create unlimited dataset in hdf5 file
        create_HDF5_unlimited_matrix_dataset_ptr(file, strDatasetName, count[0], count[1], "numeric");
        cummoffset = 0;
        
    } 
    
    // Get write position
    offset[1] = cummoffset;
    cummoffset = cummoffset + restmp.cols();
    
    DataSet* unlimDataset = new DataSet(file->openDataSet(strDatasetName));
    
    // Extend dataset before put data
    if((i%(M/k)) != 0 && cummoffset > 0) {
        extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
    }
    
    write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp)  );  
    unlimDataset->close();
}
}

        }
        
    } catch(FileIException& error) { // catch failure caused by the H5File operations
        if (bcenter==true || bscale==true){ normalizedData->close(); }
        file->close();
        ::Rf_error( "c++ exception First_level_SvdBlock_decomposition_hdf5 (File IException)" );
        return void();
    } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
        if (bcenter==true || bscale==true){ normalizedData->close(); }
        file->close();
        ::Rf_error( "c++ exception First_level_SvdBlock_decomposition_hdf5 (DataSet IException)" );
        return void();
    } catch(GroupIException& error) { // catch failure caused by the Group operations
        if (bcenter==true || bscale==true){ normalizedData->close(); }
        file->close();
        ::Rf_error( "c++ exception First_level_SvdBlock_decomposition_hdf5 (Group IException)" );
        return void();
    } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
        if (bcenter==true || bscale==true){ normalizedData->close(); }
        file->close();
        ::Rf_error( "c++ exception First_level_SvdBlock_decomposition_hdf5 (DataSpace IException)" );
        return void();
    } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
        if (bcenter==true || bscale==true){ normalizedData->close(); }
        file->close();
        ::Rf_error( "c++ exception First_level_SvdBlock_decomposition_hdf5 (Data TypeIException)" );
        return void();
    }catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        if (bcenter==true || bscale==true){ normalizedData->close(); }
        file->close();
        return void();
    }
    
    if (bcenter==true || bscale==true) { 
        normalizedData->close(); 
    }
    return void();
    
    
}

// 
// 
// // Reads small datasets from hdf5 and perform a svd descomposition from each block,
// // results are saved in hdf5 datasets under temporal group to be processed if necessary
// // int Next_level_SvdBlock_decomposition_hdf5(H5File* file, std::string strGroupName, int k, int q, 
// //                                            bool bcenter, bool bscale, Rcpp::Nullable<int> threads = R_NilValue)
// void Next_level_SvdBlock_decomposition_hdf5(H5File* file, std::string strGroupName, int k, int q,
//                                             bool bcenter, bool bscale, double dthreshold, 
//                                             Rcpp::Nullable<int> threads = R_NilValue)
// {
//     
//     Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1, 1);
//     Rcpp::IntegerVector block = Rcpp::IntegerVector::create(1, 1);
//     Rcpp::IntegerVector count = Rcpp::IntegerVector::create(0, 0);
//     Rcpp::IntegerVector offset = Rcpp::IntegerVector::create(0, 0);
//     int cummoffset = 0; //.. MODIFICAT 07/09/2020 ..//
//     int ithreads;
//     
//     int M;
//     // int nconv, M, p, n;
//     // int maxsizetoread;
//     
//     /*** TODO:
//      * 
//      * 
//      * Revisar el cas que no existeixi el dataset per a fer join perquè aquell dataset no s'ha creat perquè tots els 
//      * valors singulars << threshold !!! 
//      * 
//      * 
//      * ***/
//     
//     CharacterVector strvmatnames = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
//     
//     try {
//         
//         // Get dataset names
//         StringVector joindata =  get_dataset_names_from_group(file, strGroupName, (std::string)strvmatnames[q-1]);
//         
//         M = joindata.size();
//         
//         if(threads.isNotNull()) 
//         {
//             if (Rcpp::as<int>(threads) < std::thread::hardware_concurrency()) {
//                 ithreads = Rcpp::as<int> (threads);
//             } else {
//                 ithreads = std::thread::hardware_concurrency()/2;
//             }
//         } else{    
//             ithreads = std::thread::hardware_concurrency()/2; //omp_get_max_threads()
//         }
//         
//         
//         //.OpenMP.// omp_set_num_threads(getDTthreads(ithreads, false));
//         
// #pragma omp parallel num_threads(getDTthreads(ithreads, false))
//         
//         
// #pragma omp for ordered schedule (static)
//         
//         // Get data from M blocks in initial matrix
//         for( int i = 0; i< M ; i++) 
//         {
//             
//             std::string strDatasetName = strGroupName + "/" + strvmatnames[q] + std::to_string(i/(M/k));
//             Eigen::MatrixXd restmp;
//             Eigen::MatrixXd X;
//             
// #pragma omp critical(accessFile)         
// {      
//     //    a) Get dataset
//     DataSet currentdataset = file->openDataSet(strGroupName + "/" + joindata[i]);
//     Rcpp::IntegerVector dims_out = get_HDF5_dataset_size(currentdataset);
//     
//     X = GetCurrentBlock_hdf5( file, &currentdataset, 0, 0, dims_out[0], dims_out[1]);
// }      
// 
// /***
//  {
//  
//  svdeig retsvd;
//  
//  //    b) Get dataset svd
//  retsvd = RcppbdSVD_lapack(X, false, false, false);
//  
//  
//  //    c) U*d
//  int isize = (retsvd.d).size();
//  Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
//  d.diagonal() = retsvd.d;
//  
//  //..// restmp = Bblock_matrix_mul_parallel(retsvd.u, d, 1024, threads);
//  restmp = Bblock_matrix_mul(retsvd.u, d, 1024);
//  }
//  ***/
// 
// {
//     
//     //    b) SVD for each block
//     svdeig retsvd;
//     retsvd = RcppbdSVD_lapack(X, false, false, false);
//     
//     int size_d = (retsvd.d).size();
//     int nzeros = 0;
//     
//     if( (retsvd.d)[size_d - 1] <= dthreshold ) {
//         nzeros = 1;
//         for( int j = (size_d - 2); ( j>1 && (retsvd.d)[i] <= dthreshold ); j-- ) {
//             nzeros++;
//         }
//     }
//     
//     //    c)  U*d
//     // Create diagonal matrix from svd decomposition d
//     
//     int isize = (retsvd.d).size() - nzeros;
//     
//     if( isize < 2 ) {
//         isize = 2;
//     }
//     
//     Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
//     d.diagonal() = (retsvd.d).head(isize);
//     
//     //..// Eigen::MatrixXd restmp = Bblock_matrix_mul_parallel(retsvd.u, d, 1024,1);
//     restmp = Bblock_matrix_mul((retsvd.u).block(0, 0, (retsvd.u).rows(), isize), d, 1024);
//     
// }
// 
// //    d) Write results to dataset
// count[0] = restmp.rows();
// count[1] = restmp.cols();
// /***
// #pragma omp ordered
//  {  
//  
// #pragma omp critical(accessFile)
//  {
//  if(i%(M/k) == 0 || (i%(M/k) > 0 &&  !exists_HDF5_element_ptr(file,strDatasetName)) ) {
//  // If dataset exists --> remove dataset
//  if( exists_HDF5_element_ptr(file,strDatasetName)){
//  remove_HDF5_element_ptr(file,strDatasetName);
//  }
//  // Create unlimited dataset in hdf5 file
//  create_HDF5_unlimited_matrix_dataset_ptr(file, strDatasetName, restmp.rows(), restmp.cols(), "numeric");
//  cummoffset = 0;
//  }
//  
//  DataSet* unlimDataset = new DataSet(file->openDataSet(strDatasetName));
//  
//  // Extend dataset before put data
//  if((i%(M/k)) != 0 && cummoffset > 0) {
//  extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
//  }
//  
//  // Get write position
//  offset[1] = cummoffset;
//  cummoffset = cummoffset + restmp.cols();
//  
//  write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp)  );  
//  unlimDataset->close();
//  }
//  }      
//  ***/
// 
// #pragma omp ordered
// {
//     
// #pragma omp critical(accessFile)
// {   
//     if(i%(M/k) == 0 || ( (i%(M/k) > 0 &&  !exists_HDF5_element_ptr(file,strDatasetName)) ) ) {
//         // If dataset exists --> remove dataset
//         if( exists_HDF5_element_ptr(file,strDatasetName)) {
//             remove_HDF5_element_ptr(file,strDatasetName);
//         }
//         
//         // Create unlimited dataset in hdf5 file
//         create_HDF5_unlimited_matrix_dataset_ptr(file, strDatasetName, count[0], count[1], "numeric");
//         cummoffset = 0;
//     } 
//     
//     // Get write position
//     offset[1] = cummoffset;
//     cummoffset = cummoffset + restmp.cols();
//     
//     DataSet* unlimDataset = new DataSet(file->openDataSet(strDatasetName));
//     
//     // Extend dataset before put data
//     if((i%(M/k)) != 0 && cummoffset > 0) {
//         extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
//     }
//     
//     write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp)  );  
//     unlimDataset->close();
// }
// }
// 
// 
// 
//         }
//         
//         //..ONLY DEBUG !!!...//remove_HDF5_multiple_elements_ptr(file, strGroupName, joindata);
//         
//     } catch(FileIException& error) { // catch failure caused by the H5File operations
//         file->close();
//         ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (File IException)" );
//         return void();
//     } catch(DataSetIException& error) { // catch failure caused by the DataSet operations
//         file->close();
//         ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (DataSet IException)" );
//         return void();
//     } catch(GroupIException& error) { // catch failure caused by the Group operations
//         file->close();
//         ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (Group IException)" );
//         return void();
//     } catch(DataSpaceIException& error) { // catch failure caused by the DataSpace operations
//         file->close();
//         ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (DataSpace IException)" );
//         return void();
//     } catch(DataTypeIException& error) { // catch failure caused by the DataSpace operations
//         file->close();
//         ::Rf_error( "c++ exception Next_level_SvdBlock_decomposition_hdf5 (Data TypeIException)" );
//         return void();
//     }catch(std::exception &ex) {
//         Rcpp::Rcout<< ex.what();
//         file->close();
//         return void();
//     }
//     
//     return void();
//     
// }
// 




}

#endif // BIGDATASTATMETH_HDF5_MATRIXSVDBLOCK_HPP