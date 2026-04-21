#include <BigDataStatMeth.hpp>

/**
 * @file hdf5_getDatasetList.cpp
 * @brief Implementation of HDF5 dataset listing functionality
 * @details This file provides functionality for retrieving lists of datasets
 * from HDF5 files. The implementation supports:
 * - Listing all datasets in a group
 * - Filtering datasets by prefix
 * - Filtering datasets by suffix
 * - Safe HDF5 file operations
 * 
 * Key features:
 * - Flexible dataset filtering
 * - Error handling for HDF5 operations
 * - Memory-safe implementation
 * - Support for read-only operations
 */

/**
 * @brief Lists datasets in an HDF5 group
 * 
 * @details Implements dataset listing functionality with optional filtering
 * by prefix or suffix. The function safely handles HDF5 file operations and
 * provides comprehensive error handling.
 * 
 * Implementation features:
 * - Safe file opening and closing
 * - Optional prefix/suffix filtering
 * - Memory management for HDF5 resources
 * - Comprehensive error handling
 * 
 * @param filename Path to HDF5 file
 * @param group Group path in HDF5 file
 * @param prefix Optional prefix filter
 * 
 * @return Vector of dataset names
 * @throws H5::FileIException if there are HDF5 file operation errors
 * @throws H5::DataSetIException if there are HDF5 dataset operation errors
 * @throws std::exception for other errors
 */

//' List Datasets in HDF5 Group
//'
//' @description
//' Retrieves a list of all datasets within a specified HDF5 group, with optional
//' filtering by prefix or suffix.
//'
//' @details
//' This function provides flexible dataset listing capabilities for HDF5 files.
//' Key features:
//' 
//' * Listing options:
//'   - All datasets in a group
//'   - Datasets matching a prefix
//'   - Datasets matching a suffix
//' 
//' * Implementation features:
//'   - Safe HDF5 file operations
//'   - Memory-efficient implementation
//'   - Comprehensive error handling
//'   - Read-only access to files
//'
//' The function opens the HDF5 file in read-only mode to ensure data safety.
//'
//' @param filename  Character string. Path to the HDF5 file.
//' @param group     Character string or \code{NULL}. Group path within the
//'   HDF5 file. If \code{NULL} (default), the entire file is traversed
//'   recursively and dataset paths are returned relative to the root
//'   (e.g. \code{"INPUT/A"}, \code{"RESULTS/SVD/d"}).
//' @param prefix    Optional character string. Only return datasets whose
//'   name starts with this prefix.
//' @param recursive Logical. If \code{TRUE}, recurse into subgroups and
//'   return full relative paths. Ignored when \code{group = NULL} (always
//'   recursive). Default \code{FALSE}.
//'
//' @return Character vector containing dataset names.
//'
//' @examples
//' \donttest{
//' fn <- tempfile(fileext = ".h5")
//' X  <- hdf5_create_matrix(fn, "INPUT/A",  data = matrix(rnorm(100), 10, 10))
//' Y  <- hdf5_create_matrix(fn, "INPUT/B",  data = matrix(rnorm(100), 10, 10))
//' Z  <- hdf5_create_matrix(fn, "RESULTS/C",data = matrix(rnorm(100), 10, 10))
//'
//' # All datasets in the file (recursive from root)
//' bdgetDatasetsList_hdf5(fn)
//'
//' # Only datasets in INPUT group
//' bdgetDatasetsList_hdf5(fn, group = "INPUT")
//'
//' # INPUT group, recursive (same result here, no subgroups)
//' bdgetDatasetsList_hdf5(fn, group = "INPUT", recursive = TRUE)
//'
//' # Filter by prefix
//' bdgetDatasetsList_hdf5(fn, group = "INPUT", prefix = "A")
//'
//' hdf5_close_all()
//' unlink(fn)
//' }
//'
//' @references
//' * The HDF Group. (2000-2010). HDF5 User's Guide.
//'
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdgetDatasetsList_hdf5(std::string filename,
                             Rcpp::Nullable<std::string> group = R_NilValue,
                             Rcpp::Nullable<std::string> prefix = R_NilValue,
                             bool  recursive = false)
{
    
    // H5File* file = nullptr;
    Rcpp::StringVector groupDatasets;
    BigDataStatMeth::hdf5File* fQuery = nullptr;
    
    try
    {
        
        H5::Exception::dontPrint();
        
        std::string strgroup  = group.isNull()  ? "/"  : Rcpp::as<std::string>(group);
        std::string strprefix = prefix.isNull() ? ""   : Rcpp::as<std::string>(prefix);
        
        // When no group is specified, recurse by default so the user sees everything
        bool brecursive = group.isNull() ? true : recursive;
        
        fQuery = new BigDataStatMeth::hdf5File(filename, false);
        fQuery->openFile("r");
        
        if( fQuery->getFileptr() != nullptr) {
            groupDatasets = fQuery->getAllDatasetNames(strgroup, strprefix, brecursive);
        } else {
            delete fQuery; fQuery = nullptr;
            Rcpp::stop("c++ exception bdgetDatasetsList_hdf5 File does not exist");
            return(R_NilValue);
        }
        
        delete fQuery; fQuery = nullptr;
        
    } catch( H5::FileIException& error ) { 
        Rcpp::stop("c++ exception bdgetDatasetsList_hdf5 (File IException)");
        return(R_NilValue);
    } catch( H5::DataSetIException& error ) { 
        Rcpp::stop("c++ exception bdgetDatasetsList_hdf5 (DataSet IException)");
        return(R_NilValue);
    } catch(std::exception &ex) {
        Rcpp::stop("c++ exception bdgetDatasetsList_hdf5: " + std::string(ex.what()));
        return(R_NilValue);
    }
    
    return(groupDatasets);
    
}
