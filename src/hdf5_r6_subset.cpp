/**
 * @file hdf5_r6_subset.cpp
 * @brief R6 wrapper for subsetting HDF5 datasets
 * @details Implements data reading (subsetting) for HDF5Matrix objects
 */

#include <BigDataStatMeth.hpp>

//' Read block from HDF5 dataset (subsetting)
//'
//' @param ptr_sexp External pointer to hdf5Dataset
//' @param rows Integer vector with row indices (1-based, as in R)
//' @param cols Integer vector with column indices (1-based, as in R)
//' @return Numeric matrix with requested data
//' 
//' @details 
//' This function reads a subset of data from an HDF5 dataset.
//' Indices are 1-based (R convention) and converted internally to 0-based (C++ convention).
//' 
//' The function handles:
//' - Contiguous blocks (e.g., rows 1:10)
//' - Non-contiguous indices (e.g., rows c(1,3,5,7))
//' - Full dimensions (e.g., all rows, specific columns)
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_hdf5dataset_subset(SEXP ptr_sexp,
                                             Rcpp::IntegerVector rows,
                                             Rcpp::IntegerVector cols)
{
    try {
        // Extract pointer
        // Rcpp::XPtr<BigDataStatMeth::hdf5Dataset> ptr(ptr_sexp);
        auto* ptr = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(ptr_sexp));
        
        // if (!ptr) {
        if (ptr == nullptr) {
            Rf_error("Invalid external pointer (NULL)");
        }
        
        if (!ptr->isOpen()) {
            Rf_error("Dataset is closed");
        }
        
        // Get dataset dimensions (as R sees them)
        hsize_t nrows_total = ptr->nrows_r();
        hsize_t ncols_total = ptr->ncols_r();
        
        // Validate indices
        for (int i = 0; i < rows.length(); i++) {
            if (rows[i] < 1 || rows[i] > static_cast<hsize_t>(nrows_total)) {
                std::string msg = "Row index out of bounds: " + std::to_string(rows[i]) + 
                                  " (dataset has " + std::to_string(nrows_total) + " rows)";
                Rf_error("%s", msg.c_str());
            }
        }
        
        for (int j = 0; j < cols.length(); j++) {
            if (cols[j] < 1 || cols[j] > static_cast<hsize_t>(ncols_total)) {
                std::string msg = "Column index out of bounds: " + std::to_string(cols[j]) + 
                    " (dataset has " + std::to_string(ncols_total) + " columns)";
                Rf_error("%s", msg.c_str());
            }
        }
        
        // Check if indices are contiguous (optimized case)
        bool rows_contiguous = true;
        bool cols_contiguous = true;
        
        if (rows.length() > 1) {
            for (int i = 1; i < rows.length(); i++) {
                if (rows[i] != rows[i-1] + 1) {
                    rows_contiguous = false;
                    break;
                }
            }
        }
        
        if (cols.length() > 1) {
            for (int j = 1; j < cols.length(); j++) {
                if (cols[j] != cols[j-1] + 1) {
                    cols_contiguous = false;
                    break;
                }
            }
        }
        
        // CASE 1: Contiguous block (most efficient)
        if (rows_contiguous && cols_contiguous) {
            
            // Convert to 0-based indices
            hsize_t row_start = rows[0] - 1;
            hsize_t col_start = cols[0] - 1;
            hsize_t row_count = rows.length();
            hsize_t col_count = cols.length();
            
            // Allocate result matrix
            Rcpp::NumericMatrix result(row_count, col_count);
            
            // Remember: HDF5 stores transposed relative to R
            // R sees: nrows_r × ncols_r
            // HDF5 has: ncols_r × nrows_r (transposed)
            
            // So when R wants rows [a:b], we need columns [a:b] in HDF5
            // And when R wants cols [c:d], we need rows [c:d] in HDF5
            
            std::vector<hsize_t> offset = {col_start, row_start};  // Swapped!
            std::vector<hsize_t> count = {col_count, row_count};   // Swapped!
            std::vector<hsize_t> stride = {1, 1};
            std::vector<hsize_t> block = {1, 1};
            
            // Read data (stored column-major in HDF5, which matches R)
            ptr->readDatasetBlock(offset, count, stride, block, result.begin());
            
            // Data is already in correct layout for R
            return result;
            
        } else {
            // CASE 2: Non-contiguous indices (read row-by-row or element-by-element)
            
            Rcpp::NumericMatrix result(rows.length(), cols.length());
            
            // For non-contiguous, read each element individually
            // (could be optimized further, but this is simple and works)
            
            for (int i = 0; i < rows.length(); i++) {
                for (int j = 0; j < cols.length(); j++) {
                    
                    // Convert to 0-based
                    hsize_t row_idx = rows[i] - 1;
                    hsize_t col_idx = cols[j] - 1;
                    
                    // Read single element
                    std::vector<hsize_t> offset = {col_idx, row_idx};  // Swapped!
                    std::vector<hsize_t> count = {1, 1};
                    std::vector<hsize_t> stride = {1, 1};
                    std::vector<hsize_t> block = {1, 1};
                    
                    double value;
                    ptr->readDatasetBlock(offset, count, stride, block, &value);
                    
                    result(i, j) = value;
                }
            }
            
            return result;
        }
        
    } catch(H5::Exception& error) {
        std::string msg = "HDF5 error in subset: " + error.getDetailMsg();
        Rf_error("%s", msg.c_str());
    } catch(std::exception& ex) {
        std::string msg = "Error in subset: " + std::string(ex.what());
        Rf_error("%s", msg.c_str());
    }
}


//' Get full dataset as matrix (convenience function)
//'
//' @param ptr_sexp External pointer to hdf5Dataset
//' @return Numeric matrix with all data
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_hdf5dataset_read_all(SEXP ptr_sexp)
{
    try {
        // Rcpp::XPtr<BigDataStatMeth::hdf5Dataset> ptr(ptr_sexp);
        auto* ptr = static_cast<BigDataStatMeth::hdf5Dataset*>(R_ExternalPtrAddr(ptr_sexp));
        
        // if (!ptr) {
        if (ptr == nullptr) {
            Rf_error("Invalid external pointer (NULL)");
        }
        
        if (!ptr->isOpen()) {
            Rf_error("Dataset is closed");
        }
        
        hsize_t nrows = ptr->nrows_r();
        hsize_t ncols = ptr->ncols_r();
        
        // Create indices for all rows and columns
        Rcpp::IntegerVector rows = Rcpp::seq(1, nrows);
        Rcpp::IntegerVector cols = Rcpp::seq(1, ncols);
        
        // Use subset function (will use optimized contiguous path)
        return rcpp_hdf5dataset_subset(ptr_sexp, rows, cols);
        
    } catch(std::exception& ex) {
        std::string msg = "Error reading dataset: " + std::string(ex.what());
        Rf_error("%s", msg.c_str());
    }
}
