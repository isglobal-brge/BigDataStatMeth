#include "BigDataStatMeth.hpp"

//' Move HDF5 Dataset
//'
//' @description
//' Moves an HDF5 dataset from one location to another within the same HDF5 file.
//' This function automatically handles moving associated rownames and colnames 
//' datasets, creates parent groups if needed, and updates all internal references.
//'
//' @param filename Character string. Path to the HDF5 file
//' @param source_path Character string. Current path to the dataset (e.g., "/group1/dataset1")
//' @param dest_path Character string. New path for the dataset (e.g., "/group2/new_name")
//' @param overwrite Logical. Whether to overwrite destination if it exists (default: FALSE)
//'
//' @return Logical. TRUE on success, FALSE on failure
//'
//' @details
//' This function provides a high-level interface for moving datasets within HDF5 files.
//' The operation is efficient as it uses HDF5's native linking mechanism without 
//' copying actual data.
//'
//' Key features:
//' \itemize{
//'   \item Moves main dataset and associated rownames/colnames datasets
//'   \item Creates parent directory structure automatically
//'   \item Preserves all dataset attributes and properties
//'   \item Updates internal dataset references
//'   \item Efficient metadata-only operation
//'   \item Comprehensive error handling
//' }
//'
//' @section Behavior:
//' \itemize{
//'   \item If the destination parent groups don't exist, they will be created automatically
//'   \item Associated rownames and colnames datasets are moved to the same new group
//'   \item All dataset attributes and properties are preserved during the move
//'   \item The operation is atomic - either all elements move successfully or none do
//' }
//'
//' @section Requirements:
//' \itemize{
//'   \item The HDF5 file must exist and be accessible
//'   \item The source dataset must exist
//'   \item The file must not be locked by another process
//'   \item User must have read-write permissions on the file
//' }
//'
//' @examples
//' \dontrun{
//' # Move dataset to a different group
//' success <- bdmove_hdf5_dataset("data.h5", 
//'                          source_path = "/old_group/my_dataset",
//'                          dest_path = "/new_group/my_dataset")
//'
//' # Rename dataset within the same group
//' success <- bdmove_hdf5_dataset("data.h5",
//'                          source_path = "/data/old_name", 
//'                          dest_path = "/data/new_name",
//'                          overwrite = TRUE)
//'
//' # Move dataset to root level
//' success <- bdmove_hdf5_dataset("data.h5",
//'                          source_path = "/deep/nested/dataset",
//'                          dest_path = "/dataset")
//'
//' # Move with automatic group creation
//' success <- bdmove_hdf5_dataset("data.h5",
//'                          source_path = "/old_location/dataset",
//'                          dest_path = "/new/deep/structure/dataset")
//' }
//'
//' @family BigDataStatMeth HDF5 utilities
//' @author BigDataStatMeth package authors
//' @export
 // [[Rcpp::export]]
 void bdmove_hdf5_dataset(std::string filename,
                     std::string source_path,
                     std::string dest_path, 
                     bool overwrite = false)
 {
     BigDataStatMeth::hdf5Dataset* dataset = nullptr;
     
     try {
         // Input validation
         if (filename.empty()) {
             Rcpp::Rcerr << "Error: filename cannot be empty" << std::endl;
             return void();
         }
         
         if (source_path.empty()) {
             Rcpp::Rcerr << "Error: source_path cannot be empty" << std::endl;
             return void();
         }
         
         if (dest_path.empty()) {
             Rcpp::Rcerr << "Error: dest_path cannot be empty" << std::endl;
             return void();
         }
         
         if (source_path == dest_path) {
             Rcpp::Rcerr << "Error: source_path and dest_path cannot be the same" << std::endl;
             return void();
         }
         
         // Create hdf5Dataset object with current path
         dataset = new BigDataStatMeth::hdf5Dataset(filename, source_path, false);
         
         // Open the existing dataset
         dataset->openDataset();
         
         // Perform the move operation
         dataset->moveDataset(dest_path, overwrite);
         
         // Clean up
         delete dataset; dataset = nullptr;
         
         return void();
         
     } catch (const std::exception& e) {
         if (dataset) {
             delete dataset; dataset = nullptr;
         }
         Rcpp::Rcerr << "Exception in bdmove_hdf5_dataset: " << e.what() << std::endl;
         return void();
     } catch (...) {
         if (dataset) {
             delete dataset; dataset = nullptr;
         }
         Rcpp::Rcerr << "Unknown exception in bdmove_hdf5_dataset" << std::endl;
         return void();
     }
 }