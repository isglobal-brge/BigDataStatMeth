/**
 * @file hdf5_r6_omics.cpp
 * @brief Thin Rcpp wrappers for omics-specific HDF5 operations.
 *
 * Exposes three header-level functions as clean R-callable wrappers,
 * following the same string-based pattern used by all other R6 wrappers:
 *   - rcpp_hdf5dataset_impute_snps()
 *   - rcpp_hdf5dataset_filter_low_coverage()
 *   - rcpp_hdf5dataset_filter_maf()
 *
 * All functions open their own unique_ptr handles internally; callers
 * must close their own HDF5Matrix handle before calling (Rule 1.7).
 *
 * Return value: named list with keys file / group / dataset / n_removed
 * (n_removed = -1 for impute, which does not remove features).  The two
 * filtering wrappers additionally return "kept": the 1-based indices, along
 * the filtered axis, of the elements that survived the filter.  The R6 layer
 * uses them to carry the dimension names of the input over to the result.
 */

#include "BigDataStatMeth.hpp"
#include "hdf5Utilities/hdf5ImputeData.hpp"
#include "hdf5Utilities/hdf5RemoveLowData.hpp"
#include "hdf5Omics/hdf5RemoveMAF.hpp"


/**
 * @brief Impute missing SNP values in an HDF5 dataset.
 *
 * Fills the entries flagged as missing (the value 3 in 0/1/2/3-coded data) by
 * computing column or row means of the non-missing values. The result can
 * overwrite the input dataset (same group/dataset) or be written to a new
 * location.
 *
 * @param in_file      Path to the HDF5 file.
 * @param in_group     Group of the input dataset.
 * @param in_dataset   Name of the input dataset.
 * @param out_group    Group for the output dataset.
 * @param out_dataset  Name for the output dataset.
 * @param by_cols      If true, impute by columns; otherwise by rows.
 * @param threads      Number of OpenMP threads (-1 = auto).
 * @param overwrite    Overwrite existing output dataset.
 * @return Named list: file, group, dataset, n_removed (always -1).
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_impute_snps(
        std::string in_file,
        std::string in_group,
        std::string in_dataset,
        std::string out_group,
        std::string out_dataset,
        bool        by_cols   = true,
        int         threads   = -1,
        bool        overwrite = false,
        Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsIn(
            new BigDataStatMeth::hdf5Dataset(in_file, in_group, in_dataset, false));
        dsIn->openDataset();
        if (!dsIn->getDatasetptr())
            throw std::runtime_error("cannot open '" + in_group + "/" +
                                     in_dataset + "'");

        std::unique_ptr<BigDataStatMeth::hdf5DatasetInternal> dsOut(
            new BigDataStatMeth::hdf5DatasetInternal(
                in_file, out_group, out_dataset, overwrite));
        dsOut->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsIn->getCompressionLevel());

        // Rcpp::Nullable<int> nth = (threads < 0) ?
        //     Rcpp::Nullable<int>(R_NilValue) :
        //     Rcpp::Nullable<int>(threads);

        Rcpp::Nullable<int> nth = R_NilValue;
        if (threads >= 0) {
            nth = Rcpp::wrap(threads);
        }
        
        BigDataStatMeth::Rcpp_Impute_snps_hdf5(
            dsIn.get(), dsOut.get(), by_cols, out_dataset, nth);

        return Rcpp::List::create(
            Rcpp::Named("file")      = in_file,
            Rcpp::Named("group")     = out_group,
            Rcpp::Named("dataset")   = out_dataset,
            Rcpp::Named("n_removed") = -1
        );

    } catch (H5::FileIException& e) {
        Rf_error("rcpp_hdf5dataset_impute_snps (File): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("rcpp_hdf5dataset_impute_snps (DataSet): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_impute_snps: %s", e.what());
    }
    return R_NilValue;
}


/**
 * @brief Remove columns or rows with too many missing values.
 *
 * Missing values are the entries equal to 3, the missing-data code of the
 * 0/1/2/3 genotype encoding used across the omics functions; R \c NA / \c NaN
 * entries are not counted. A column (by_cols = true) or row (by_cols = false)
 * is removed when its proportion of 3s is greater than or equal to \p pcent.
 * Input and output datasets must be different.
 *
 * @param in_file      Path to the HDF5 file.
 * @param in_group     Group of the input dataset.
 * @param in_dataset   Name of the input dataset.
 * @param out_group    Group for the output dataset.
 * @param out_dataset  Name for the output dataset.
 * @param pcent        Missing-data threshold [0,1]. Default 0.05.
 * @param by_cols      If true, filter columns; otherwise rows.
 * @param overwrite    Overwrite existing output dataset.
 * @return Named list: file, group, dataset, n_removed, kept.
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_filter_low_coverage(
        std::string in_file,
        std::string in_group,
        std::string in_dataset,
        std::string out_group,
        std::string out_dataset,
        double      pcent     = 0.05,
        bool        by_cols   = true,
        bool        overwrite = false,
        Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        if (in_group + "/" + in_dataset == out_group + "/" + out_dataset)
            throw std::runtime_error("input and output dataset must be different");

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsIn(
            new BigDataStatMeth::hdf5Dataset(in_file, in_group, in_dataset, false));
        dsIn->openDataset();
        if (!dsIn->getDatasetptr())
            throw std::runtime_error("cannot open '" + in_group + "/" +
                                     in_dataset + "'");

        std::unique_ptr<BigDataStatMeth::hdf5DatasetInternal> dsOut(
            new BigDataStatMeth::hdf5DatasetInternal(
                in_file, out_group, out_dataset, overwrite));
        dsOut->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsIn->getCompressionLevel());

        std::vector<hsize_t> kept;

        int n_removed = BigDataStatMeth::Rcpp_Remove_Low_Data_hdf5(
            dsIn.get(), dsOut.get(), by_cols, pcent, &kept);

        // Surviving positions, 1-based, along the filtered axis: R columns when
        // by_cols is true (they are HDF5 dimension 0), R rows otherwise.
        Rcpp::IntegerVector keptR(kept.size());
        for (size_t k = 0; k < kept.size(); k++)
            keptR[k] = static_cast<int>(kept[k]) + 1;

        return Rcpp::List::create(
            Rcpp::Named("file")      = in_file,
            Rcpp::Named("group")     = out_group,
            Rcpp::Named("dataset")   = out_dataset,
            Rcpp::Named("n_removed") = std::abs(n_removed),
            Rcpp::Named("kept")      = keptR
        );

    } catch (H5::FileIException& e) {
        Rf_error("rcpp_hdf5dataset_filter_low_coverage (File): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("rcpp_hdf5dataset_filter_low_coverage (DataSet): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_filter_low_coverage: %s",
                 e.what());
    }
    return R_NilValue;
}


/**
 * @brief Remove SNPs whose Minor Allele Frequency is at or below a threshold.
 *
 * MAF formula (for 0/1/2-coded diploid data):
 *   maf = (n0/n) + 0.5*(n1/n);  if maf > 0.5: maf = 1 - maf
 * SNPs with maf <= maf_threshold are removed, i.e. only the SNPs whose MAF is
 * strictly above the threshold are kept.
 *
 * @param in_file       Path to the HDF5 file.
 * @param in_group      Group of the input dataset.
 * @param in_dataset    Name of the input dataset.
 * @param out_group     Group for the output dataset.
 * @param out_dataset   Name for the output dataset.
 * @param maf_threshold MAF threshold [0, 0.5]. Default 0.05.
 * @param by_cols       If true, process columns; otherwise rows.
 * @param block_size    I/O block size. Default 100.
 * @param overwrite     Overwrite existing output dataset.
 * @return Named list: file, group, dataset, n_removed, kept.
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_filter_maf(
        std::string in_file,
        std::string in_group,
        std::string in_dataset,
        std::string out_group,
        std::string out_dataset,
        double      maf_threshold = 0.05,
        bool        by_cols       = false,
        int         block_size    = 100,
        bool        overwrite     = false,
        Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        if (in_group + "/" + in_dataset == out_group + "/" + out_dataset)
            throw std::runtime_error("input and output dataset must be different");

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsIn(
            new BigDataStatMeth::hdf5Dataset(in_file, in_group, in_dataset, false));
        dsIn->openDataset();
        if (!dsIn->getDatasetptr())
            throw std::runtime_error("cannot open '" + in_group + "/" +
                                     in_dataset + "'");

        std::unique_ptr<BigDataStatMeth::hdf5DatasetInternal> dsOut(
            new BigDataStatMeth::hdf5DatasetInternal(
                in_file, out_group, out_dataset, overwrite));
        dsOut->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsIn->getCompressionLevel());

        std::vector<hsize_t> kept;

        int n_removed = BigDataStatMeth::Rcpp_Remove_MAF_hdf5(
            dsIn.get(), dsOut.get(), by_cols, maf_threshold, block_size, &kept);

        // Surviving positions, 1-based, along the filtered axis: R columns when
        // by_cols is true (they are HDF5 dimension 0), R rows otherwise.
        Rcpp::IntegerVector keptR(kept.size());
        for (size_t k = 0; k < kept.size(); k++)
            keptR[k] = static_cast<int>(kept[k]) + 1;

        return Rcpp::List::create(
            Rcpp::Named("file")      = in_file,
            Rcpp::Named("group")     = out_group,
            Rcpp::Named("dataset")   = out_dataset,
            Rcpp::Named("n_removed") = std::abs(n_removed),
            Rcpp::Named("kept")      = keptR
        );

    } catch (H5::FileIException& e) {
        Rf_error("rcpp_hdf5dataset_filter_maf (File): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("rcpp_hdf5dataset_filter_maf (DataSet): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_filter_maf: %s", e.what());
    }
    return R_NilValue;
}
