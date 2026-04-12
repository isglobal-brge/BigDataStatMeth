/**
 * @file hdf5_r6_write.cpp
 * @brief R6 wrapper for writing data to HDF5 datasets
 */

#include <BigDataStatMeth.hpp>

//' Write data block to HDF5 dataset (R6 wrapper)
//'
//' @description
//' Writes a block of data to an HDF5 dataset at specified offset.
//' Supports writing scalars, vectors, and matrices.
//'
//' @param ptr_sexp External pointer (SEXP) to hdf5Dataset
//' @param value Data to write (numeric scalar, vector, or matrix)
//' @param row_offset Starting row (0-based in C++, but receives 1-based from R)
//' @param col_offset Starting column (0-based in C++, but receives 1-based from R)
//' @param nrows Number of rows to write
//' @param ncols Number of columns to write
//'
//' @return NULL (invisible)
//'
//' @keywords internal
// [[Rcpp::export]]
void rcpp_hdf5dataset_write_block(SEXP ptr_sexp,
                                   Rcpp::RObject value,
                                   int row_offset,
                                   int col_offset,
                                   int nrows,
                                   int ncols)
{
    try {
        H5::Exception::dontPrint();

        auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(
                       R_ExternalPtrAddr(ptr_sexp));

        if (ds == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!ds->isOpen())
            throw std::runtime_error("Dataset is closed");

        // Convert from R's 1-based to C++'s 0-based indexing
        hsize_t offset_r = static_cast<hsize_t>(row_offset - 1);
        hsize_t offset_c = static_cast<hsize_t>(col_offset - 1);
        hsize_t count_r  = static_cast<hsize_t>(nrows);
        hsize_t count_c  = static_cast<hsize_t>(ncols);

        // Validate bounds
        if (offset_r + count_r > ds->nrows_r() || 
            offset_c + count_c > ds->ncols_r()) {
            throw std::runtime_error("Write indices out of bounds");
        }

        // Prepare offset and count vectors for writeDatasetBlock
        // NOTE: writeDatasetBlock expects HDF5 file order (transposed from R)
        std::vector<hsize_t> vOffset = {offset_c, offset_r};  // swap for HDF5
        std::vector<hsize_t> vCount  = {count_c, count_r};    // swap for HDF5
        std::vector<hsize_t> vStride = {1, 1};
        std::vector<hsize_t> vBlock  = {1, 1};

        // writeDatasetBlock handles the transposition internally
        // bTranspose = false means "data is in R order, please transpose for HDF5"
        ds->writeDatasetBlock(value, vOffset, vCount, vStride, vBlock, false);

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error writing block: %s", e.what());
    }
}


//' Write entire dataset (R6 wrapper)
//'
//' @description
//' Replaces entire HDF5 dataset contents with new data.
//'
//' @param ptr_sexp External pointer (SEXP) to hdf5Dataset
//' @param value Data to write (numeric matrix)
//'
//' @return NULL (invisible)
//'
//' @keywords internal
// [[Rcpp::export]]
void rcpp_hdf5dataset_write_all(SEXP ptr_sexp, Rcpp::RObject value)
{
    try {
        H5::Exception::dontPrint();

        auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(
                       R_ExternalPtrAddr(ptr_sexp));

        if (ds == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!ds->isOpen())
            throw std::runtime_error("Dataset is closed");

        // Validate dimensions match
        Rcpp::IntegerVector dims;
        if (Rcpp::is<Rcpp::NumericMatrix>(value)) {
            Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(value);
            dims = Rcpp::IntegerVector::create(mat.nrow(), mat.ncol());
        } else if (Rcpp::is<Rcpp::IntegerMatrix>(value)) {
            Rcpp::IntegerMatrix mat = Rcpp::as<Rcpp::IntegerMatrix>(value);
            dims = Rcpp::IntegerVector::create(mat.nrow(), mat.ncol());
        } else {
            throw std::runtime_error("Value must be a numeric or integer matrix");
        }

        if (dims[0] != static_cast<int>(ds->nrows_r()) || 
            dims[1] != static_cast<int>(ds->ncols_r())) {
            throw std::runtime_error("Matrix dimensions don't match dataset");
        }

        // Write entire dataset
        ds->writeDataset(value);

    } catch (H5::FileIException& e) {
        Rf_error("HDF5 file error: %s", e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("HDF5 dataset error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("Error writing dataset: %s", e.what());
    }
}


// ---------------------------------------------------------------------------
// rcpp_hdf5_create_matrix
// ---------------------------------------------------------------------------

//' Create an HDF5 dataset with configurable compression (R6 wrapper)
//'
//' @description
//' Creates an HDF5 dataset of size \code{nrows x ncols} and optionally writes
//' data to it. Replaces \code{bdCreate_hdf5_matrix()} /
//' \code{bdCreate_hdf5_emptyDataset()} in the R6+S3 interface so that
//' compression can be controlled from R.
//'
//' @param filename Character. Path to the HDF5 file.
//' @param group Character. Group path inside the file.
//' @param dataset Character. Dataset name.
//' @param nrows Integer. Number of rows (>= 1).
//' @param ncols Integer. Number of columns (>= 1).
//' @param data Optional numeric/integer matrix or data.frame; NULL creates
//'   an empty (zero-filled) dataset.
//' @param dtype Character. Element type: "real" (default), "int", "logical".
//' @param overwrite_file Logical. Recreate file if it already exists.
//' @param overwrite_dataset Logical. Replace dataset if it already exists.
//' @param compression Integer 0-9. gzip compression level (0 = no compression,
//'   6 = balanced default). Applied to the new dataset only.
//'
//' @return Named list with \code{filename} and \code{path} of the created dataset.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5_create_matrix(std::string filename,
                                    std::string group,
                                    std::string dataset,
                                    int nrows,
                                    int ncols,
                                    Rcpp::Nullable<Rcpp::RObject> data = R_NilValue,
                                    std::string dtype                   = "real",
                                    bool overwrite_file                 = false,
                                    bool overwrite_dataset              = false,
                                    int compression                     = 6)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "",
        Rcpp::Named("path")     = "");

    BigDataStatMeth::hdf5File*    objFile    = nullptr;
    BigDataStatMeth::hdf5Dataset* objDataset = nullptr;
    BigDataStatMeth::hdf5Dims*    dsdims     = nullptr;

    try {
        H5::Exception::dontPrint();

        if (nrows < 1 || ncols < 1)
            Rf_error("rcpp_hdf5_create_matrix: nrows and ncols must be >= 1");

        objFile = new BigDataStatMeth::hdf5File(filename, overwrite_file);
        int iRes = objFile->createFile();

        if (iRes != EXEC_OK && iRes != EXEC_WARNING)
            Rf_error("rcpp_hdf5_create_matrix: could not create/open file");

        if (iRes == EXEC_WARNING)
            objFile->openFile("rw");

        objDataset = new BigDataStatMeth::hdf5Dataset(
            objFile, group, dataset, overwrite_dataset);

        // Apply compression before dataset creation
        objDataset->setCompressionLevel(compression);
        objDataset->createDataset(nrows, ncols, dtype);

        if (data.isNotNull()) {
            Rcpp::RObject robj = Rcpp::as<Rcpp::RObject>(data);

            if (Rf_inherits(robj, "data.frame")) {
                SEXP mat = Rcpp::Language("as.matrix", robj).eval();
                if (Rf_isMatrix(mat))
                    objDataset->writeDataset(Rcpp::as<Rcpp::NumericMatrix>(mat));
                else
                    Rf_error("rcpp_hdf5_create_matrix: cannot coerce data.frame to matrix");
            } else {
                objDataset->writeDataset(robj);
            }

            // Write dimnames if present
            SEXP dn = robj.attr("dimnames");
            if (!Rf_isNull(dn)) {
                Rcpp::List dimnames(dn);
                if (dimnames.size() >= 2) {
                    dsdims = new BigDataStatMeth::hdf5Dims(objDataset);
                    Rcpp::CharacterVector svrows(1), svcols(1);
                    if (!Rf_isNull(dimnames[0]))
                        svrows = Rcpp::as<Rcpp::CharacterVector>(dimnames[0]);
                    if (!Rf_isNull(dimnames[1]))
                        svcols = Rcpp::as<Rcpp::CharacterVector>(dimnames[1]);
                    dsdims->writeDimnames(
                        (svrows.size() < nrows) ? Rcpp::CharacterVector(1) : svrows,
                        (svcols.size()  < ncols) ? Rcpp::CharacterVector(1) : svcols);
                    delete dsdims; dsdims = nullptr;
                }
            }
        }

        delete objDataset; objDataset = nullptr;
        delete objFile;    objFile    = nullptr;

        lst["filename"] = filename;
        lst["path"]     = group + "/" + dataset;

    } catch (H5::FileIException& e) {
        if (dsdims)   { delete dsdims;   dsdims   = nullptr; }
        if (objFile)  { delete objFile;  objFile  = nullptr; }
        BigDataStatMeth::checkClose_file(objDataset);
        Rf_error("HDF5 file error (rcpp_hdf5_create_matrix): %s",
                 e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        if (dsdims)   { delete dsdims;   dsdims   = nullptr; }
        if (objFile)  { delete objFile;  objFile  = nullptr; }
        BigDataStatMeth::checkClose_file(objDataset);
        Rf_error("HDF5 dataset error (rcpp_hdf5_create_matrix): %s",
                 e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        if (dsdims)   { delete dsdims;   dsdims   = nullptr; }
        if (objFile)  { delete objFile;  objFile  = nullptr; }
        BigDataStatMeth::checkClose_file(objDataset);
        Rf_error("Error (rcpp_hdf5_create_matrix): %s", e.what());
    }

    return lst;
}
