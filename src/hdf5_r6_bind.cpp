/**
 * @file hdf5_r6_bind.cpp
 * @brief Rcpp wrappers for block-wise cbind / rbind of HDF5Matrix objects.
 *
 * HDF5 coordinate system note (critical):
 *   R stores matrices in column-major; BigDataStatMeth stores them transposed
 *   in HDF5 (row-major).  Therefore:
 *     - HDF5 dim[0]  =  R ncols   (hdf5Dataset::dim()[0])
 *     - HDF5 dim[1]  =  R nrows   (hdf5Dataset::dim()[1])
 *   All readDatasetBlock / writeDatasetBlock calls use HDF5 internal coords,
 *   NOT R-visible nrows_r() / ncols_r().
 *
 * Block strategy:
 *   Iterate in chunks along HDF5 dim[1] (= R rows), keeping RAM bounded at
 *   block_rows * max(dim[0]_A, dim[0]_B) * 8 bytes regardless of matrix size.
 *
 * cbind  (R: same nrows, append ncols => HDF5: same dim[1], append dim[0])
 * rbind  (R: same ncols, append nrows => HDF5: same dim[0], append dim[1])
 */

#include "BigDataStatMeth.hpp"

// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_bind(
        std::string file_a,
        std::string group_a,
        std::string dataset_a,
        std::string file_b,
        std::string group_b,
        std::string dataset_b,
        std::string out_file,
        std::string out_group,
        std::string out_dataset,
        std::string func,        // "cbind" or "rbind"
        bool        overwrite  = false,
        int         block_rows = 1000,
        Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        if (func != "cbind" && func != "rbind")
            Rcpp::stop("rcpp_hdf5dataset_bind: func must be 'cbind' or 'rbind'");

        const bool is_cbind = (func == "cbind");

        // ── Open inputs ──────────────────────────────────────────────────
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(file_a, group_a, dataset_a, false));
        dsA->openDataset();
        if (!dsA->getDatasetptr())
            Rcpp::stop("rcpp_hdf5dataset_bind: cannot open A '%s/%s'",
                       group_a.c_str(), dataset_a.c_str());

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(file_b, group_b, dataset_b, false));
        dsB->openDataset();
        if (!dsB->getDatasetptr())
            Rcpp::stop("rcpp_hdf5dataset_bind: cannot open B '%s/%s'",
                       group_b.c_str(), dataset_b.c_str());

        // ── HDF5 internal dimensions (NOT R-visible dimensions) ───────────
        //   dim()[0] = R ncols,  dim()[1] = R nrows
        const hsize_t dA0 = dsA->dim()[0];   // R ncols of A
        const hsize_t dA1 = dsA->dim()[1];   // R nrows of A
        const hsize_t dB0 = dsB->dim()[0];   // R ncols of B
        const hsize_t dB1 = dsB->dim()[1];   // R nrows of B

        // ── Dimension validation ─────────────────────────────────────────
        if (is_cbind && dA1 != dB1)
            Rcpp::stop("cbind: nrow(A)=%llu != nrow(B)=%llu",
                       (unsigned long long)dA1, (unsigned long long)dB1);
        if (!is_cbind && dA0 != dB0)
            Rcpp::stop("rbind: ncol(A)=%llu != ncol(B)=%llu",
                       (unsigned long long)dA0, (unsigned long long)dB0);

        // ── Output dimensions ────────────────────────────────────────────
        //   cbind: createDataset(nrows_R = dA1,       ncols_R = dA0+dB0)
        //          => HDF5 dim[0] = dA0+dB0, dim[1] = dA1
        //   rbind: createDataset(nrows_R = dA1+dB1,   ncols_R = dA0)
        //          => HDF5 dim[0] = dA0,     dim[1] = dA1+dB1
        const hsize_t out_nrows_R = is_cbind ? dA1       : dA1 + dB1;
        const hsize_t out_ncols_R = is_cbind ? dA0 + dB0 : dA0;

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsOut(
            new BigDataStatMeth::hdf5Dataset(out_file, out_group, out_dataset, overwrite));
        dsOut->setCompressionLevel(compression.isNotNull() ? Rcpp::as<int>(compression) : dsA->getCompressionLevel());
        dsOut->createDataset(out_nrows_R, out_ncols_R, "real");

        // ── Block copy helper ─────────────────────────────────────────────
        //   Copies src (HDF5: d0 rows x d1 cols) into dsOut block-wise.
        //   hdf5_row_off : offset along HDF5 dim[0] in destination (cbind)
        //   hdf5_col_off : offset along HDF5 dim[1] in destination (rbind)
        //
        //   Follows the same read/write pattern as RcppBind_datasets_hdf5:
        //     readDatasetBlock({0, col}, {d0, chunk})
        //     Eigen::Map<RowMajor>(buf, d0, chunk) -> mat
        //     writeDatasetBlock(wrap(mat), offset, {chunk, d0}, ..., bTranspose=true)
        //
        //   With bTranspose=true, writeDatasetBlock internally swaps count to
        //   {d0, chunk} and transposes data -> correct HDF5 hyperslab write.
        const hsize_t bsz = static_cast<hsize_t>(block_rows);
        std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};

        auto copy_blocks = [&](BigDataStatMeth::hdf5Dataset* src,
                                hsize_t d0, hsize_t d1,
                                hsize_t hdf5_row_off, hsize_t hdf5_col_off)
        {
            hsize_t col = 0;
            while (col < d1) {
                const hsize_t chunk = std::min(bsz, d1 - col);

                std::vector<double> buf(d0 * chunk);
                src->readDatasetBlock({0, col}, {d0, chunk},
                                      stride, blk, buf.data());

                Eigen::Map<Eigen::Matrix<double,
                    Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
                    mat(buf.data(), static_cast<int>(d0),
                                    static_cast<int>(chunk));

                // count = {mat.cols(), mat.rows()} = {chunk, d0}
                // bTranspose=true -> hsCount swapped to {d0, chunk}
                // -> HDF5 hyperslab at {hdf5_row_off, col+hdf5_col_off}
                //    of size {d0, chunk}  ✓
                dsOut->writeDatasetBlock(
                    Rcpp::wrap(mat),
                    {hdf5_row_off, col + hdf5_col_off},
                    {chunk, d0},
                    stride, blk, true);

                col += chunk;
            }
        };

        if (is_cbind) {
            // A -> output at HDF5 row offset 0   (R cols [0 .. dA0-1])
            copy_blocks(dsA.get(), dA0, dA1, 0,   0);
            // B -> output at HDF5 row offset dA0 (R cols [dA0 .. dA0+dB0-1])
            copy_blocks(dsB.get(), dB0, dB1, dA0, 0);
        } else {
            // A -> output at HDF5 col offset 0   (R rows [0 .. dA1-1])
            copy_blocks(dsA.get(), dA0, dA1, 0, 0);
            // B -> output at HDF5 col offset dA1 (R rows [dA1 .. dA1+dB1-1])
            copy_blocks(dsB.get(), dB0, dB1, 0, dA1);
        }

        return Rcpp::List::create(
            Rcpp::Named("file")    = out_file,
            Rcpp::Named("group")   = out_group,
            Rcpp::Named("dataset") = out_dataset
        );

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_bind (File): %s", e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_bind (DataSet): %s", e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_bind: %s", e.what());
    }
    return R_NilValue;
}
