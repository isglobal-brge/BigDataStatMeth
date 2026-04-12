/**
 * @file hdf5_r6_normalize.cpp
 * @brief Rcpp wrapper for block-wise HDF5 matrix normalization (R6 ptr API).
 *
 * Follows the same pattern as hdf5_r6_arithmetic.cpp: receives filename /
 * group / dataset strings instead of external pointers, and opens all
 * datasets internally with unique_ptr.  This avoids HDF5 metadata-cache
 * conflicts entirely — no open-handle checks are needed.
 */

// #include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(Rhdf5lib)]]

#include "BigDataStatMeth.hpp"

/**
 * @brief Block-wise normalization of an HDF5 matrix.
 *
 * Opens \p in_file / \p in_group / \p in_dataset internally and writes the
 * normalized result to \p out_file / \p out_group / \p out_dataset.
 * All datasets are managed via \c unique_ptr — no external pointer is held
 * during the operation, so any number of R-side \c HDF5Matrix objects may
 * remain open without risk of metadata-cache conflicts.
 *
 * @param in_file    Path to the HDF5 file containing the input dataset.
 * @param in_group   HDF5 group of the input dataset.
 * @param in_dataset Name of the input dataset.
 * @param out_file   Path to the HDF5 file for the output dataset.
 * @param out_group  HDF5 group for the output dataset.
 * @param out_dataset Name of the output dataset.
 * @param center     Subtract column/row means (default TRUE).
 * @param scale      Divide by column/row standard deviations (default TRUE).
 * @param byrows     Normalize by rows if TRUE; by columns (default) if FALSE.
 * @param wsize      Block size for HDF5 reads (NULL = auto).
 * @return Named list: \code{center} and \code{scale} numeric vectors so the
 *   R layer can attach them as \code{scaled:center} / \code{scaled:scale}.
 *
 * @export
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_normalize(std::string in_file,
                                       std::string in_group,
                                       std::string in_dataset,
                                       std::string out_file,
                                       std::string out_group,
                                       std::string out_dataset,
                                       bool center = true,
                                       bool scale  = true,
                                       bool byrows = false,
                                       Rcpp::Nullable<int> wsize       = R_NilValue,
                                       Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        // ── Open input dataset ────────────────────────────────────────────
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(in_file, in_group, in_dataset, false)
        );
        dsA->openDataset();

        if (dsA->getDatasetptr() == nullptr)
            Rf_error("rcpp_hdf5dataset_normalize: cannot open input dataset "
                     "'%s/%s' in '%s'.",
                     in_group.c_str(), in_dataset.c_str(), in_file.c_str());

        const hsize_t nrows = dsA->nrows();
        const hsize_t ncols = dsA->ncols();

        // ── Compute mean and SD in one pass ───────────────────────────────
        Eigen::MatrixXd datanormal;
        if (!byrows) {
            datanormal = Eigen::MatrixXd::Zero(2, nrows);
            BigDataStatMeth::get_HDF5_mean_sd_by_column(dsA.get(), datanormal,
                                                        true, true, wsize);
        } else {
            datanormal = Eigen::MatrixXd::Zero(2, ncols);
            BigDataStatMeth::get_HDF5_mean_sd_by_row(dsA.get(), datanormal,
                                                     true, true, wsize);
        }

        // ── Adjust scale factor to match base R scale() contract ──────────
        // base R scale(center=FALSE, scale=TRUE) uses
        //   RMS = sqrt( sum(x^2) / (n-1) )
        // Our stat function computes
        //   SD  = sqrt( sum((x-mu)^2) / (n-1) )
        // Relationship: RMS^2 = SD^2 + (n/(n-1)) * mu^2
        // nrows() = HDF5 first dim = R ncols;  ncols() = HDF5 second dim = R nrows.
        if (!center && scale) {
            // n = number of observations per variable:
            //   column-wise (!byrows): each R-column has ncols (= R nrows) observations
            //   row-wise    ( byrows): each R-row    has nrows (= R ncols) observations
            const double n  = static_cast<double>(!byrows ? ncols : nrows);
            const double nf = n / (n - 1.0);
            datanormal.row(1) = (datanormal.row(1).array().square() +
                                 nf * datanormal.row(0).array().square()).sqrt();
        }

        // ── Pre-create output file if it is different from the input file ────
        // hdf5Dataset constructor calls checkHDF5File() -> H5::H5File::isHdf5()
        // which throws when the file does not exist yet.  Creating it here with
        // H5F_ACC_TRUNC is safe: if the file already exists it is simply opened.
        if (out_file != in_file) {
            std::ifstream probe(out_file);
            if (!probe.good()) {
                H5::H5File newf(out_file, H5F_ACC_TRUNC);
                newf.close();
            }
        }

        // ── Create output dataset ─────────────────────────────────────────
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsOut(
            new BigDataStatMeth::hdf5Dataset(out_file, out_group, out_dataset, true)
        );
        dsOut->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());
        dsOut->createDataset(dsA.get(), "real");

        if (dsOut->getDatasetptr() == nullptr)
            Rf_error("rcpp_hdf5dataset_normalize: failed to create output "
                     "dataset '%s/%s' in '%s'.",
                     out_group.c_str(), out_dataset.c_str(), out_file.c_str());

        // ── Normalize block-wise ──────────────────────────────────────────
        bool bcorrected = false;
        BigDataStatMeth::RcppNormalizeHdf5(dsA.get(), dsOut.get(),
                                            datanormal, wsize,
                                            center, scale, byrows, bcorrected);

        // ── Return center and scale vectors to R ──────────────────────────
        Rcpp::NumericVector r_center = Rcpp::wrap(datanormal.row(0));
        Rcpp::NumericVector r_scale  = Rcpp::wrap(datanormal.row(1));

        return Rcpp::List::create(
            Rcpp::Named("center") = center ? r_center : Rcpp::NumericVector(0),
            Rcpp::Named("scale")  = scale  ? r_scale  : Rcpp::NumericVector(0)
        );

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_normalize (File IException): %s",
                 e.getDetailMsg().c_str());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_normalize (DataSet IException): %s",
                 e.getDetailMsg().c_str());
    } catch (H5::DataSpaceIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_normalize (DataSpace IException): %s",
                 e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_normalize: %s", e.what());
    }
    return R_NilValue;
}
