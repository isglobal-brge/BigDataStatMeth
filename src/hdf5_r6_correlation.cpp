/**
 * @file hdf5_r6_correlation.cpp
 * @brief Rcpp wrapper for HDF5 matrix correlation with cross-file output support.
 *
 * Delegates computation entirely to the existing bdCorr_hdf5() engine
 * (which always writes results into the input file). When the caller requests
 * a different output file, this wrapper uses hdf5Dataset::copyToFile() to
 * transfer the correlation (and optional p-values) dataset block-by-block —
 * no full matrix is ever loaded into RAM.
 *
 * All dataset handles are opened with unique_ptr; no external R pointer is
 * received, so HDF5 metadata-cache conflicts cannot occur.
 */

#include "BigDataStatMeth.hpp"

/**
 * @brief Block-wise HDF5 matrix correlation with optional cross-file output.
 *
 * Runs bdCorr_hdf5() to compute the correlation (always written to in_file_x),
 * then copies the result block-by-block to out_file via hdf5Dataset::copyToFile()
 * when out_file differs from in_file_x.
 *
 * @param in_file_x          Path to the HDF5 file containing matrix X.
 * @param in_group_x         Group of dataset X.
 * @param in_dataset_x       Name of dataset X.
 * @param in_file_y          Path to file containing matrix Y ("" = single-matrix).
 * @param in_group_y         Group of dataset Y ("" = single-matrix).
 * @param in_dataset_y       Name of dataset Y ("" = single-matrix).
 * @param out_file           Output HDF5 file. "" = same as in_file_x.
 * @param out_group          Output group. "" = auto ("CORR/<dataset>").
 * @param trans_x            Transpose X before correlating.
 * @param trans_y            Transpose Y before correlating.
 * @param method             "pearson" or "spearman".
 * @param use_complete_obs   Use only complete (non-NA) observations.
 * @param compute_pvalues    Also compute and store p-values.
 * @param block_size         Block size for bdCorr_hdf5 computation.
 * @param copy_blockrows     Row-block size for cross-file copy (default 500).
 * @return Named list with keys: filename, group, correlation, method,
 *   correlation_type, n_variables, n_observations, use_complete_obs,
 *   pvalues, has_pvalues. All paths refer to the output file.
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_cor(std::string in_file_x,
                                 std::string in_group_x,
                                 std::string in_dataset_x,
                                 std::string in_file_y       = "",
                                 std::string in_group_y      = "",
                                 std::string in_dataset_y    = "",
                                 std::string out_file        = "",
                                 std::string out_group       = "",
                                 bool        trans_x         = false,
                                 bool        trans_y         = false,
                                 std::string method          = "pearson",
                                 bool        use_complete_obs = false,
                                 bool        compute_pvalues  = true,
                                 int         block_size       = 1000,
                                 int         copy_blockrows   = 500,
                                 Rcpp::Nullable<int> threads  = R_NilValue,
                                 Rcpp::Nullable<int> compression = R_NilValue)
{
    try {
        H5::Exception::dontPrint();

        // ── Resolve output file ───────────────────────────────────────────
        // bdCorr_hdf5 always writes to in_file_x when output_filename == "".
        // We pass "" so it always writes there; we copy afterwards if needed.
        const bool cross_file = (!out_file.empty() && out_file != in_file_x);

        // // ── Run correlation via existing bdCorr_hdf5 engine ───────────────
        // Rcpp::List res = bdCorr_hdf5(
        //     in_file_x, in_group_x, in_dataset_x,
        //     in_file_y, in_group_y, in_dataset_y,
        //     trans_x, trans_y,
        //     method, use_complete_obs, compute_pvalues,
        //     block_size,
        //     /*overwrite=*/true,
        //     /*output_filename=*/"",
        //     out_group,
        //     /*output_dataset_corr=*/"",
        //     /*output_dataset_pval=*/"",
        //     /*threads=*/-1
        // );

        
        // ── Nullable wrappers para los parámetros opcionales ─────────────────
        Rcpp::Nullable<Rcpp::CharacterVector> n_out_group =
            out_group.empty() ? R_NilValue : Rcpp::wrap(out_group);
        Rcpp::Nullable<Rcpp::CharacterVector> n_out_corr  = R_NilValue;
        Rcpp::Nullable<Rcpp::CharacterVector> n_out_pval  = R_NilValue;
        Rcpp::Nullable<int> n_threads = threads;
        
        // ── Single o cross según si hay dataset Y ────────────────────────────
        Rcpp::List res;
        const bool is_cross = !in_dataset_y.empty();
        
        int icompression = compression.isNotNull() ? Rcpp::as<int>(compression) : 6;
        
        if (!is_cross) {
            res = BigDataStatMeth::RcppbdCorr_hdf5_single(
                in_file_x, in_group_x, in_dataset_x,
                method, use_complete_obs, compute_pvalues,
                block_size, /*overwrite=*/true,
                n_out_group, n_out_corr, n_out_pval,
                !trans_x, n_threads, icompression
            );
        } else {
            res = BigDataStatMeth::RcppbdCorr_hdf5_cross(
                in_file_x, in_group_x, in_dataset_x,
                in_file_y, in_group_y, in_dataset_y,
                method, use_complete_obs, compute_pvalues,
                block_size, /*overwrite=*/true,
                in_file_x,   // output_filename = in_file_x (siempre)
                n_out_group, n_out_corr, n_out_pval,
                !trans_x, !trans_y, n_threads, icompression
            );
        }
        
        // ── If same file (or no output file specified) — return as-is ─────
        if (!cross_file) {
            return res;
        }

        // ── Cross-file: copy correlation (and pvalues) to out_file ────────
        const std::string src_file  = Rcpp::as<std::string>(res["filename"]);
        const std::string src_group = Rcpp::as<std::string>(res["group"]);
        const std::string corr_name = Rcpp::as<std::string>(res["correlation"]);
        const bool        has_pval  = Rcpp::as<bool>       (res["has_pvalues"]);
        const std::string pval_name = Rcpp::as<std::string>(res["pvalues"]);

        // Determine destination group: use out_group if set, else mirror src_group
        const std::string dst_group = out_group.empty() ? src_group : out_group;

        // Copy correlation dataset
        {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsSrc(
                new BigDataStatMeth::hdf5Dataset(src_file, src_group, corr_name, false)
            );
            dsSrc->openDataset();
            if (dsSrc->getDatasetptr() == nullptr)
                Rf_error("rcpp_hdf5dataset_cor: cannot open correlation result "
                         "'%s/%s' in '%s'",
                         src_group.c_str(), corr_name.c_str(), src_file.c_str());

            dsSrc->copyToFile(out_file, dst_group, corr_name,
                              static_cast<hsize_t>(copy_blockrows), true);
        }

        // Copy p-values dataset (if computed)
        if (has_pval && !pval_name.empty()) {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsPval(
                new BigDataStatMeth::hdf5Dataset(src_file, src_group, pval_name, false)
            );
            dsPval->openDataset();
            if (dsPval->getDatasetptr() != nullptr) {
                dsPval->copyToFile(out_file, dst_group, pval_name,
                                   static_cast<hsize_t>(copy_blockrows), true);
            }
        }

        // ── Return list pointing to the output file ───────────────────────
        return Rcpp::List::create(
            Rcpp::Named("filename")         = out_file,
            Rcpp::Named("group")            = dst_group,
            Rcpp::Named("correlation")      = corr_name,
            Rcpp::Named("method")           = res["method"],
            Rcpp::Named("correlation_type") = res["correlation_type"],
            Rcpp::Named("n_variables")      = res["n_variables"],
            Rcpp::Named("n_observations")   = res["n_observations"],
            Rcpp::Named("use_complete_obs") = res["use_complete_obs"],
            Rcpp::Named("pvalues")          = res["pvalues"],
            Rcpp::Named("has_pvalues")      = res["has_pvalues"]
        );

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_cor (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_cor (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_cor: %s", e.what());
    }
    return R_NilValue;
}
