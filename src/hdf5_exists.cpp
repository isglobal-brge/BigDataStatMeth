/**
 * @file hdf5_exists.cpp
 * @brief Check if HDF5 element (group or dataset) exists
 * @details Wrapper for exists_HDF5_element() used in tests and validations
 */

#include <BigDataStatMeth.hpp>

//' Check if HDF5 element exists
//'
//' @description
//' Verifies if a group or dataset exists within an HDF5 file.
//' Used for validation and conditional logic in HDF5 operations.
//'
//' @param filename Character. Path to HDF5 file.
//' @param path Character. Full path to element (e.g., "data/matrix" or "group1/group2").
//'   Can be a group path or a full dataset path.
//'
//' @return Logical. \code{TRUE} if element exists, \code{FALSE} otherwise.
//'   Returns \code{FALSE} if file doesn't exist or cannot be opened.
//'
//' @details
//' This function opens the HDF5 file in read-only mode, checks if the
//' specified path exists, then closes the file automatically via destructor.
//' It's safe to call on files that may not exist.
//'
//' **Note:** This checks for **both** groups and datasets. If you need to
//' distinguish between them, use \code{\link{bdExists_hdf5_element}} with
//' additional validation.
//'
//' @examples
//' \donttest{
//' # Create test file
//' bdCreate_hdf5_matrix("test.h5", matrix(1:100, 10, 10), 
//'                      "data", "X", overwriteFile = TRUE)
//'
//' # Check if elements exist
//' hdf5_exists("test.h5", "data")         # TRUE (group exists)
//' hdf5_exists("test.h5", "data/X")       # TRUE (dataset exists)
//' hdf5_exists("test.h5", "data/Y")       # FALSE (doesn't exist)
//' hdf5_exists("test.h5", "missing/path") # FALSE (doesn't exist)
//' hdf5_exists("nonexistent.h5", "data")  # FALSE (file doesn't exist)
//'
//' # Cleanup
//' unlink("test.h5")
//' }
//'
//' @seealso \code{\link{bdExists_hdf5_element}} for the bd* equivalent
//'
//' @export
// [[Rcpp::export]]
bool hdf5_exists(std::string filename, std::string path)
{
    try {
        H5::Exception::dontPrint();
        
        // Check if file exists first (using standard C++)
        std::ifstream file_check(filename);
        if (!file_check.good()) {
            return false;
        }
        file_check.close();
        
        // Create hdf5File object (stack allocation - destructor handles cleanup)
        BigDataStatMeth::hdf5File file_obj(filename, false);
        
        // Open file in read-only mode
        // Returns H5::H5File* pointer, or throws on error
        H5::H5File* pfile = file_obj.openFile("r");
        
        if (pfile == nullptr) {
            return false;
        }
        
        // Check if element exists using the API function
        // exists_HDF5_element(H5::H5File* file, std::string element)
        bool exists = BigDataStatMeth::exists_HDF5_element(pfile, path);
        
        // No need to explicitly close - destructor of file_obj handles it
        return exists;
        
    } catch (const H5::FileIException& e) {
        // File errors → element doesn't exist (or file invalid)
        return false;
    } catch (const H5::Exception& e) {
        // HDF5 errors → element doesn't exist
        return false;
    } catch (const std::runtime_error& e) {
        // exists_HDF5_element throws runtime_error on problems
        return false;
    } catch (const std::exception& e) {
        // Other errors → element doesn't exist
        return false;
    } catch (...) {
        // Unknown errors → element doesn't exist
        return false;
    }
}
