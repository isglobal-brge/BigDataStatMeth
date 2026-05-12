/**
 * @file hdf5_r6_aggregations.cpp
 * @brief R6 wrappers for block-wise aggregation operations on HDF5 datasets
 *
 * @details
 * Provides Rcpp-exported entry points used by the HDF5Matrix R6/S3 interface
 * for standard R aggregation generics:
 *
 *  Column-wise (result length = ncols_R):
 *   - rcpp_hdf5dataset_colSums   → colSums(X)
 *   - rcpp_hdf5dataset_colMeans  → colMeans(X)
 *   - rcpp_hdf5dataset_colMins   → apply(X, 2, min)
 *   - rcpp_hdf5dataset_colMaxs   → apply(X, 2, max)
 *   - rcpp_hdf5dataset_colVars   → apply(X, 2, var)
 *   - rcpp_hdf5dataset_colSds    → apply(X, 2, sd)
 *
 *  Row-wise (result length = nrows_R):
 *   - rcpp_hdf5dataset_rowSums   → rowSums(X)
 *   - rcpp_hdf5dataset_rowMeans  → rowMeans(X)
 *   - rcpp_hdf5dataset_rowMins   → apply(X, 1, min)
 *   - rcpp_hdf5dataset_rowMaxs   → apply(X, 1, max)
 *   - rcpp_hdf5dataset_rowVars   → apply(X, 1, var)
 *   - rcpp_hdf5dataset_rowSds    → apply(X, 1, sd)
 *
 *  Scalar (whole matrix as flat vector):
 *   - rcpp_hdf5dataset_scalar_sum   → sum(X)
 *   - rcpp_hdf5dataset_scalar_mean  → mean(X)
 *   - rcpp_hdf5dataset_scalar_min   → min(X)
 *   - rcpp_hdf5dataset_scalar_max   → max(X)
 *   - rcpp_hdf5dataset_scalar_var   → var(as.vector(X))
 *   - rcpp_hdf5dataset_scalar_sd    → sd(as.vector(X))
 *
 * All functions follow the same conventions as hdf5_r6_arithmetic.cpp:
 *  - Input is an external pointer (SEXP) to an open hdf5Dataset object.
 *  - ALL logic inside try{} block.
 *  - Validation throws std::runtime_error (never Rcpp::stop inside try{}).
 *  - A fresh hdf5Dataset is opened from the same file for the actual work,
 *    leaving the caller's handle untouched.
 *  - Errors use Rf_error() in catch{} to avoid SIGABRT on macOS ARM64.
 *
 * @note
 * CRITICAL: Rcpp::stop() must only be called inside catch{} handlers.
 * Use Rf_error() (not Rcpp::stop()) there to avoid memory corruption
 * (double-free / SIGABRT) on macOS ARM64 with Rcpp >= 1.1.0.
 */

#include <BigDataStatMeth.hpp>

// ---------------------------------------------------------------------------
// Internal helper: extract raw pointer from SEXP and open a fresh dataset
// ---------------------------------------------------------------------------

/** @brief Open a fresh read-only hdf5Dataset from the metadata in ptr. */
static std::unique_ptr<BigDataStatMeth::hdf5Dataset>
agg_open_dataset(SEXP ptr)
{
    auto* raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
        R_ExternalPtrAddr(ptr));

    if (raw == nullptr)
        throw std::runtime_error("Invalid external pointer (NULL)");
    if (!raw->isOpen())
        throw std::runtime_error("Dataset is closed");

    //.. 20260304 ..// std::unique_ptr<BigDataStatMeth::hdf5Dataset> ds( new BigDataStatMeth::hdf5Dataset( raw->getFileName(), raw->getGroup(), raw->getDatasetName(), false));
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> ds( new BigDataStatMeth::hdf5Dataset( raw->getFullPath(), raw->getGroup(), raw->getDatasetName(), false));
    ds->openDataset();

    if (ds->getDatasetptr() == nullptr)
        throw std::runtime_error("Failed to open dataset: " +
                                 raw->getGroup() + "/" + raw->getDatasetName());
    return ds;
}


// ===========================================================================
// Column-wise wrappers
// ===========================================================================

//' Column sums of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{colSums(X)} for an HDF5
//' matrix referenced by an external pointer.
//'
//' @param ptr     External pointer (SEXP) to an open hdf5Dataset.
//' @param paral   Logical or NULL; enable OpenMP parallelisation.
//' @param wsize   Integer or NULL; block size (NULL = auto).
//' @param threads Integer or NULL; thread count (NULL = auto).
//'
//' @return Numeric vector of length ncols_R.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_colSums(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_colSums(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_colSums: %s", e.what());
    }
}

//' Column means of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{colMeans(X)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length ncols_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_colMeans(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_colMeans(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_colMeans: %s", e.what());
    }
}

//' Column minimums of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{apply(X, 2, min)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length ncols_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_colMins(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_colMins(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_colMins: %s", e.what());
    }
}

//' Column maximums of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{apply(X, 2, max)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length ncols_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_colMaxs(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_colMaxs(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_colMaxs: %s", e.what());
    }
}

//' Column variances of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{apply(X, 2, var)}.
//' Uses Bessel's correction (n-1).
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length ncols_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_colVars(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_colVars(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_colVars: %s", e.what());
    }
}

//' Column standard deviations of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{apply(X, 2, sd)}.
//' Uses Bessel's correction (n-1).
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length ncols_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_colSds(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_colSds(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_colSds: %s", e.what());
    }
}


// ===========================================================================
// Row-wise wrappers
// ===========================================================================

//' Row sums of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{rowSums(X)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length nrows_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_rowSums(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_rowSums(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_rowSums: %s", e.what());
    }
}

//' Row means of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{rowMeans(X)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length nrows_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_rowMeans(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_rowMeans(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_rowMeans: %s", e.what());
    }
}

//' Row minimums of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{apply(X, 1, min)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length nrows_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_rowMins(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_rowMins(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_rowMins: %s", e.what());
    }
}

//' Row maximums of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{apply(X, 1, max)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length nrows_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_rowMaxs(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_rowMaxs(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_rowMaxs: %s", e.what());
    }
}

//' Row variances of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{apply(X, 1, var)}.
//' Uses Bessel's correction (n-1).
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length nrows_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_rowVars(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_rowVars(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_rowVars: %s", e.what());
    }
}

//' Row standard deviations of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise, OpenMP-parallel computation of \code{apply(X, 1, sd)}.
//' Uses Bessel's correction (n-1).
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Numeric vector of length nrows_R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_rowSds(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return Rcpp::wrap(
            BigDataStatMeth::get_HDF5_rowSds(ds.get(), bparal, wsize, threads));
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_rowSds: %s", e.what());
    }
}


// ===========================================================================
// Scalar wrappers
// ===========================================================================

//' Sum of all elements of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise computation of \code{sum(X)}.  Equivalent to
//' \code{sum(as.matrix(X))} but without loading the full matrix into RAM.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Scalar numeric.
//' @keywords internal
// [[Rcpp::export]]
double rcpp_hdf5dataset_scalar_sum(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return BigDataStatMeth::get_HDF5_scalar_sum(ds.get(), bparal, wsize, threads);
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_scalar_sum: %s", e.what());
    }
}

//' Mean of all elements of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise computation of \code{mean(X)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Scalar numeric.
//' @keywords internal
// [[Rcpp::export]]
double rcpp_hdf5dataset_scalar_mean(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return BigDataStatMeth::get_HDF5_scalar_mean(ds.get(), bparal, wsize, threads);
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_scalar_mean: %s", e.what());
    }
}

//' Minimum of all elements of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise computation of \code{min(X)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Scalar numeric.
//' @keywords internal
// [[Rcpp::export]]
double rcpp_hdf5dataset_scalar_min(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return BigDataStatMeth::get_HDF5_scalar_min(ds.get(), bparal, wsize, threads);
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_scalar_min: %s", e.what());
    }
}

//' Maximum of all elements of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise computation of \code{max(X)}.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Scalar numeric.
//' @keywords internal
// [[Rcpp::export]]
double rcpp_hdf5dataset_scalar_max(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return BigDataStatMeth::get_HDF5_scalar_max(ds.get(), bparal, wsize, threads);
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_scalar_max: %s", e.what());
    }
}

//' Variance of all elements of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise computation of \code{var(as.vector(X))}.
//' Uses Bessel's correction (N-1) where N is the total number of elements.
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Scalar numeric.
//' @keywords internal
// [[Rcpp::export]]
double rcpp_hdf5dataset_scalar_var(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return BigDataStatMeth::get_HDF5_scalar_var(ds.get(), bparal, wsize, threads);
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_scalar_var: %s", e.what());
    }
}

//' Standard deviation of all elements of an HDF5 dataset (R6 wrapper)
//'
//' @description
//' Block-wise computation of \code{sd(as.vector(X))}.
//' Uses Bessel's correction (N-1).
//'
//' @inheritParams rcpp_hdf5dataset_colSums
//' @return Scalar numeric.
//' @keywords internal
// [[Rcpp::export]]
double rcpp_hdf5dataset_scalar_sd(
        SEXP                 ptr,
        Rcpp::Nullable<bool> paral   = R_NilValue,
        Rcpp::Nullable<int>  wsize   = R_NilValue,
        Rcpp::Nullable<int>  threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        auto ds = agg_open_dataset(ptr);
        const bool bparal = paral.isNotNull() && Rcpp::as<bool>(paral);
        return BigDataStatMeth::get_HDF5_scalar_sd(ds.get(), bparal, wsize, threads);
    } catch (std::exception& e) {
        Rf_error("rcpp_hdf5dataset_scalar_sd: %s", e.what());
    }
}
