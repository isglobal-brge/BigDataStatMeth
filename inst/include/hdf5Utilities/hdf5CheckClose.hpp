/**
 * @file hdf5CheckClose.hpp
 * @brief Safe file closing utilities for HDF5 datasets
 * 
 * This file provides a set of overloaded functions for safely closing HDF5 datasets.
 * It implements comprehensive error handling and ensures proper cleanup of resources
 * even in error conditions. The functions support closing multiple datasets
 * simultaneously with proper null pointer checking.
 * 
 * Key features:
 * - Safe dataset closing
 * - Comprehensive error handling
 * - Support for multiple datasets
 * - Null pointer checking
 * - Resource cleanup
 * 
 * @note This module is part of the BigDataStatMeth library.
 * @note This module is part of the BigDataStatMeth library and is not intended to be 
 * used directly by the user. It is used internally by the library.
 */

#ifndef BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP
#define BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP

namespace BigDataStatMeth {

    /**
     * @brief Safely closes a single HDF5 dataset
     * 
     * @param ds1 Pointer to the dataset to close
     * 
     * @throws H5::FileIException on file operation errors
     * @throws H5::DataSetIException on dataset operation errors
     * @throws H5::DataSpaceIException on dataspace operation errors
     * 
     * @note Checks for null pointer before attempting to close
     */
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } 
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"c++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"C++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    
    /**
     * @brief Safely closes two HDF5 datasets
     * 
     * @param ds1 Pointer to the first dataset
     * @param ds2 Pointer to the second dataset
     * 
     * @throws H5::FileIException on file operation errors
     * @throws H5::DataSetIException on dataset operation errors
     * @throws H5::DataSpaceIException on dataspace operation errors
     * 
     * @note Checks each pointer for null before attempting to close
     */
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                         BigDataStatMeth::hdf5Dataset* ds2) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } 
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"c++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"C++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    
    /**
     * @brief Safely closes three HDF5 datasets
     * 
     * @param ds1 Pointer to the first dataset
     * @param ds2 Pointer to the second dataset
     * @param ds3 Pointer to the third dataset
     * 
     * @throws H5::FileIException on file operation errors
     * @throws H5::DataSetIException on dataset operation errors
     * @throws H5::DataSpaceIException on dataspace operation errors
     * 
     * @note Checks each pointer for null before attempting to close
     */
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                         BigDataStatMeth::hdf5Dataset* ds2, 
                         BigDataStatMeth::hdf5Dataset* ds3) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"c++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"C++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    /**
     * @brief Safely closes four HDF5 datasets
     * 
     * @param ds1 Pointer to the first dataset
     * @param ds2 Pointer to the second dataset
     * @param ds3 Pointer to the third dataset
     * @param ds4 Pointer to the fourth dataset
     * 
     * @throws H5::FileIException on file operation errors
     * @throws H5::DataSetIException on dataset operation errors
     * @throws H5::DataSpaceIException on dataspace operation errors
     * 
     * @note Checks each pointer for null before attempting to close
     */
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                         BigDataStatMeth::hdf5Dataset* ds2, 
                         BigDataStatMeth::hdf5Dataset* ds3, 
                         BigDataStatMeth::hdf5Dataset* ds4) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } 
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"c++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"C++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    
    
    
    /**
     * @brief Safely closes five HDF5 datasets
     * 
     * @param ds1 Pointer to the first dataset
     * @param ds2 Pointer to the second dataset
     * @param ds3 Pointer to the third dataset
     * @param ds4 Pointer to the fourth dataset
     * @param ds5 Pointer to the fifth dataset
     * 
     * @throws H5::FileIException on file operation errors
     * @throws H5::DataSetIException on dataset operation errors
     * @throws H5::DataSpaceIException on dataspace operation errors
     * 
     * @note Checks each pointer for null before attempting to close
     */
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                                BigDataStatMeth::hdf5Dataset* ds2, 
                                BigDataStatMeth::hdf5Dataset* ds3, 
                                BigDataStatMeth::hdf5Dataset* ds4, 
                                BigDataStatMeth::hdf5Dataset* ds5) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } else if(ds5->getDatasetptr() !=nullptr){
                ds5->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"c++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"C++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }



    /**
     * @brief Safely closes six HDF5 datasets
     * 
     * @param ds1 Pointer to the first dataset
     * @param ds2 Pointer to the second dataset
     * @param ds3 Pointer to the third dataset
     * @param ds4 Pointer to the fourth dataset
     * @param ds5 Pointer to the fifth dataset
     * @param ds6 Pointer to the sixth dataset
     * 
     * @throws H5::FileIException on file operation errors
     * @throws H5::DataSetIException on dataset operation errors
     * @throws H5::DataSpaceIException on dataspace operation errors
     * 
     * @note Checks each pointer for null before attempting to close
     */
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                                       BigDataStatMeth::hdf5Dataset* ds2, 
                                       BigDataStatMeth::hdf5Dataset* ds3, 
                                       BigDataStatMeth::hdf5Dataset* ds4, 
                                       BigDataStatMeth::hdf5Dataset* ds5,
                                       BigDataStatMeth::hdf5Dataset* ds6) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } else if(ds5->getDatasetptr() !=nullptr){
                ds5->close_file();
            } else if(ds6->getDatasetptr() !=nullptr){
                ds6->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"c++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"C++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    /**
     * @brief Safely closes seven HDF5 datasets
     * 
     * @param ds1 Pointer to the first dataset
     * @param ds2 Pointer to the second dataset
     * @param ds3 Pointer to the third dataset
     * @param ds4 Pointer to the fourth dataset
     * @param ds5 Pointer to the fifth dataset
     * @param ds6 Pointer to the sixth dataset
     * @param ds7 Pointer to the seventh dataset
     * 
     * @throws H5::FileIException on file operation errors
     * @throws H5::DataSetIException on dataset operation errors
     * @throws H5::DataSpaceIException on dataspace operation errors
     * 
     * @note Checks each pointer for null before attempting to close
     */
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                                       BigDataStatMeth::hdf5Dataset* ds2, 
                                       BigDataStatMeth::hdf5Dataset* ds3, 
                                       BigDataStatMeth::hdf5Dataset* ds4, 
                                       BigDataStatMeth::hdf5Dataset* ds5,
                                       BigDataStatMeth::hdf5Dataset* ds6,
                                       BigDataStatMeth::hdf5Dataset* ds7) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } else if(ds5->getDatasetptr() !=nullptr){
                ds5->close_file();
            } else if(ds6->getDatasetptr() !=nullptr){
                ds6->close_file();
            } else if(ds7->getDatasetptr() !=nullptr){
                ds7->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"c++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"C++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }


    /**
     * @brief Safely closes eight HDF5 datasets
     * 
     * @param ds1 Pointer to the first dataset
     * @param ds2 Pointer to the second dataset
     * @param ds3 Pointer to the third dataset
     * @param ds4 Pointer to the fourth dataset
     * @param ds5 Pointer to the fifth dataset
     * @param ds6 Pointer to the sixth dataset
     * @param ds7 Pointer to the seventh dataset
     * @param ds8 Pointer to the eighth dataset
     * 
     * @throws H5::FileIException on file operation errors
     * @throws H5::DataSetIException on dataset operation errors
     * @throws H5::DataSpaceIException on dataspace operation errors
     * 
     * @note Checks each pointer for null before attempting to close
     */
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                                       BigDataStatMeth::hdf5Dataset* ds2, 
                                       BigDataStatMeth::hdf5Dataset* ds3, 
                                       BigDataStatMeth::hdf5Dataset* ds4, 
                                       BigDataStatMeth::hdf5Dataset* ds5,
                                       BigDataStatMeth::hdf5Dataset* ds6,
                                       BigDataStatMeth::hdf5Dataset* ds7,
                                       BigDataStatMeth::hdf5Dataset* ds8) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } else if(ds5->getDatasetptr() !=nullptr){
                ds5->close_file();
            } else if(ds6->getDatasetptr() !=nullptr){
                ds6->close_file();
            } else if(ds7->getDatasetptr() !=nullptr){
                ds7->close_file();
            } else if(ds8->getDatasetptr() !=nullptr){
                ds8->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"c++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"c++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"C++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }

}

#endif // BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP