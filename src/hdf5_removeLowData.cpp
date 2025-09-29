#include <BigDataStatMeth.hpp>
#include "hdf5Utilities/hdf5RemoveLowData.hpp"

/**
 * @file hdf5_removeLowData.cpp
 * @brief Implementation of low data removal for HDF5-stored genomic data
 * @details This file provides functionality for removing SNPs (Single Nucleotide
 * Polymorphisms) with low representation in genomic data stored in HDF5 format.
 * The implementation supports:
 * - Row-wise and column-wise filtering
 * - Configurable threshold for removal
 * - Flexible output options
 * - Memory-efficient operations
 * 
 * Key features:
 * - Support for large genomic datasets
 * - Configurable filtering direction
 * - Threshold-based filtering
 * - Memory-efficient implementation
 * - Comprehensive error handling
 */

/**
 * @brief Removes SNPs with low representation from HDF5 dataset
 * 
 * @details Implements efficient removal of SNPs with low representation in
 * genomic data stored in HDF5 format. The function supports both row-wise
 * and column-wise filtering with configurable thresholds.
 * 
 * Implementation features:
 * - Flexible filtering direction (row/column)
 * - Configurable removal threshold
 * - Memory-efficient operations
 * - Safe file operations
 * - Comprehensive error handling
 * 
 * @param filename Path to HDF5 file
 * @param group Input group containing dataset
 * @param dataset Input dataset name
 * @param outgroup Output group for results
 * @param outdataset Output dataset name
 * @param pcent Threshold percentage for removal
 * @param bycols Whether to filter by columns
 * @param overwrite Whether to overwrite existing dataset
 * 
 * @throws H5::FileIException for HDF5 file operation errors
 * @throws H5::DataSetIException for HDF5 dataset operation errors
 * @throws H5::DataSpaceIException for HDF5 dataspace errors
 * @throws H5::DataTypeIException for HDF5 datatype errors
 * @throws std::exception for other errors
 */

//' Remove Low-Representation SNPs from HDF5 Dataset
//'
//' @description
//' Removes SNPs (Single Nucleotide Polymorphisms) with low representation from
//' genomic data stored in HDF5 format.
//'
//' @details
//' This function provides efficient filtering capabilities for genomic data with
//' support for:
//' 
//' * Filtering options:
//'   - Row-wise or column-wise filtering
//'   - Configurable threshold percentage
//'   - Flexible output location
//' 
//' * Implementation features:
//'   - Memory-efficient processing
//'   - Safe file operations
//'   - Comprehensive error handling
//'   - Progress reporting
//'
//' The function supports both in-place modification and creation of new datasets.
//'
//' @param filename Character string. Path to the HDF5 file.
//' @param group Character string. Path to the group containing input dataset.
//' @param dataset Character string. Name of the dataset to filter.
//' @param outgroup Character string. Output group path for filtered data.
//' @param outdataset Character string. Output dataset name for filtered data.
//' @param pcent Numeric (optional). Threshold percentage for removal (0-1).
//'   Default is 0.5. SNPs with representation below this threshold are removed.
//' @param bycols Logical (optional). Whether to filter by columns (TRUE) or
//'   rows (FALSE). Default is TRUE.
//' @param overwrite Logical (optional). Whether to overwrite existing dataset.
//'   Default is FALSE.
//'
//' @return No return value, called for side effects (data filtering).
//'   Prints a warning message indicating the number of rows/columns removed.
//'
//' @examples
//' \dontrun{
//' library(BigDataStatMeth)
//' 
//' # Create test SNP data with missing values
//' snps <- matrix(sample(c(0, 1, 2, NA), 100, replace = TRUE,
//'                      prob = c(0.3, 0.3, 0.3, 0.1)), 10, 10)
//' 
//' # Save to HDF5
//' fn <- "snp_data.hdf5"
//' bdCreate_hdf5_matrix(fn, snps, "genotype", "raw_snps",
//'                      overwriteFile = TRUE)
//' 
//' # Remove SNPs with low representation
//' bdRemovelowdata_hdf5(
//'   filename = fn,
//'   group = "genotype",
//'   dataset = "raw_snps",
//'   outgroup = "genotype_filtered",
//'   outdataset = "filtered_snps",
//'   pcent = 0.3,
//'   bycols = TRUE
//' )
//' 
//' # Cleanup
//' if (file.exists(fn)) {
//'   file.remove(fn)
//' }
//' }
//'
//' @references
//' * The HDF Group. (2000-2010). HDF5 User's Guide.
//' * Marchini, J., & Howie, B. (2010). Genotype imputation for genome-wide
//'   association studies. Nature Reviews Genetics, 11(7), 499-511.
//'
//' @seealso
//' * \code{\link{bdImputeSNPs_hdf5}} for imputing missing SNP values
//' * \code{\link{bdCreate_hdf5_matrix}} for creating HDF5 matrices
//'
//' @export
// [[Rcpp::export]]
void bdRemovelowdata_hdf5( std::string filename, std::string group, std::string dataset, std::string outgroup, std::string outdataset, 
                               Rcpp::Nullable<double> pcent, Rcpp::Nullable<bool> bycols, Rcpp::Nullable<bool> overwrite = R_NilValue)
{
    
    BigDataStatMeth::hdf5Dataset* dsIn = nullptr;
    BigDataStatMeth::hdf5DatasetInternal* dsOut = nullptr;

    try
    {
        
        bool bcols, bforce;
        double dpcent;
        int iremoved = 0;
        
        std::string stroutdata = outgroup +"/" + outdataset;
        std::string strdataset = group +"/" + dataset;
        
        if(bycols.isNull()){  bcols = true ; }
        else{  bcols = Rcpp::as<bool>(bycols); }
        
        if(pcent.isNull()){   dpcent = 0.5 ; }
        else{     dpcent = Rcpp::as<double>(pcent); }
        
        if(overwrite.isNull()) { bforce = false ; }
        else { bforce = Rcpp::as<bool>(overwrite); }
        

        if( strdataset.compare(stroutdata)!= 0)
        {
            dsIn = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
            dsIn->openDataset();
            
            dsOut = new BigDataStatMeth::hdf5DatasetInternal(filename, outgroup, outdataset, bforce);
            
        } else {
            Rcpp::Rcerr << "c++ exception bdRemovelowdata_hdf5: " << "Input and output dataset must be different";
            return void();
        }
        
        if( dsIn->getDatasetptr() != nullptr) {
            iremoved = Rcpp_Remove_Low_Data_hdf5( dsIn, dsOut, bcols, dpcent);
        } else {
            checkClose_file(dsIn, dsOut);
            Rcpp::Rcerr << "c++ exception bdRemovelowdata_hdf5: " << "File does not exist";
            return void();
        }
        
        Rcpp::Function warning("warning");
        if (bycols )
            warning( std::to_string(iremoved) + " Columns have been removed");
        else
            warning( std::to_string(iremoved) + " Rows have been removed");
    
        delete dsIn; dsIn = nullptr;
        delete dsOut; dsOut = nullptr;
        
    } catch( H5::FileIException& error ){  
        checkClose_file(dsIn, dsOut);
        Rcpp::stop ("c++ exception bdRemovelowdata_hdf5 (File IException)");
        return void();
    } catch( H5::DataSetIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::stop ("c++ exception bdRemovelowdata_hdf5 (DataSet IException)");
        return void();
    } catch( H5::DataSpaceIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::stop ("c++ exception bdRemovelowdata_hdf5 (DataSpace IException)");
        return void();
    } catch( H5::DataTypeIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::stop( "c++ exception bdRemovelowdata_hdf5 (DataType IException)");
        return void();
    } catch(std::exception &ex) {
        checkClose_file(dsIn, dsOut);
        Rcpp::stop ("c++ exception bdRemovelowdata_hdf5: %s",ex.what());
        return void();
    } catch (...) {
        checkClose_file(dsIn, dsOut);
        Rcpp::stop ("c++ exception bdRemovelowdata_hdf5 (unknown reason)");
        return void();
    }
    
    return void();
}
