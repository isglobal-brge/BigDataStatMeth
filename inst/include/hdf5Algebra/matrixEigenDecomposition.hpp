/**
 * @file matrixEigenDecomposition.hpp
 * @brief Eigenvalue decomposition algorithms for big data matrices stored in HDF5 format
 * 
 * This header provides efficient eigenvalue and eigenvector computation functions
 * for large matrices stored in HDF5 format using the Spectra library (same version
 * as used in BigDataStatMeth SVD) for consistent results with RSpectra.
 * 
 * @author BigDataStatMeth Development Team
 * @date 2025
 * @version 1.0
 * @note 2026-03-07 Output datasets now inherit compression level from input datasets
 *         via setCompressionLevel() called before every createDataset() invocation.
 * @note Updated for Spectra 1.0.1 compatibility
 * 
 * @details
 * The implementation uses Spectra library for eigenvalue computation, which provides:
 * - Support for both symmetric and non-symmetric matrices
 * - Consistent results with RSpectra package
 * - Memory-efficient algorithms for large matrices
 * - Integration with existing BigDataStatMeth infrastructure
 */

#ifndef BIGDATASTATMETH_HDF5_MATRIXEIGEN_HPP
#define BIGDATASTATMETH_HDF5_MATRIXEIGEN_HPP

#include "Spectra/SymEigsSolver.h"
#include "Spectra/GenEigsSolver.h"
#include "Spectra/MatOp/DenseSymMatProd.h"
#include "Spectra/MatOp/DenseGenMatProd.h"

#include <random>


namespace BigDataStatMeth {
    
    /**
     * @brief Convert 'which' string to Spectra SortRule for symmetric matrices
     */
    inline Spectra::SortRule getSymmetricSortRule(const std::string& which) {
        if (which == "LA") return Spectra::SortRule::LargestAlge;
        if (which == "SA") return Spectra::SortRule::SmallestAlge;
        if (which == "LM") return Spectra::SortRule::LargestMagn;
        if (which == "SM") return Spectra::SortRule::SmallestMagn;
        return Spectra::SortRule::LargestMagn;
    }
    
    /**
     * @brief Convert 'which' string to Spectra SortRule for general matrices
     */
    inline Spectra::SortRule getGeneralSortRule(const std::string& which) {
        if (which == "LM") return Spectra::SortRule::LargestMagn;
        if (which == "SM") return Spectra::SortRule::SmallestMagn;
        if (which == "LR") return Spectra::SortRule::LargestReal;
        if (which == "SR") return Spectra::SortRule::SmallestReal;
        if (which == "LI") return Spectra::SortRule::LargestImag;
        if (which == "SI") return Spectra::SortRule::SmallestImag;
        return Spectra::SortRule::LargestMagn;
    }
    
    /**
     * @brief Validate and adjust Spectra parameters for convergence
     * @note Only call this for partial decomposition (k < n). Full decomposition
     *       (k >= n) must use SelfAdjointEigenSolver, not Spectra.
     */
    inline std::tuple<int, int> validateSpectraParams(int n, int k, int ncv) {
        if (k <= 0) k = std::min(n, 6);
        if (k >= n) k = n - 1;
        
        if (ncv <= 0)
            ncv = std::min(n, std::max(2 * k + 1, k + 2));
        
        ncv = std::max(ncv, k + 2);
        ncv = std::min(ncv, n);
        if (ncv - k < 2)
            ncv = std::min(n, k + 2);
        
        return std::make_tuple(k, ncv);
    }
    
    /**
     * @brief Matrix symmetry detection
     */
    inline bool isMatrixSymmetric(const Eigen::MatrixXd& X, int sample_size = 100) {
        if (X.rows() != X.cols()) return false;
        int n = X.rows();
        if (n <= 3) return true;
        if (n <= 50) {
            double max_diff = (X - X.transpose()).cwiseAbs().maxCoeff();
            return max_diff <= 1e-12 * std::max(1.0, X.cwiseAbs().maxCoeff());
        }
        sample_size = std::min(sample_size, n * n / 4);
        Eigen::VectorXd diffs(sample_size);
        int count = 0;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, n - 1);
        for (int s = 0; s < sample_size && count < sample_size; ++s) {
            int i = dis(gen), j = dis(gen);
            if (i != j) diffs(count++) = std::abs(X(i,j) - X(j,i));
        }
        if (count == 0) return true;
        return diffs.head(count).maxCoeff() <= 1e-12 * std::max(1.0, X.cwiseAbs().maxCoeff());
    }
    
    /**
     * @brief Structure to hold eigendecomposition results
     */
    struct eigdecomp {
        Eigen::VectorXd eigenvalues_real;
        Eigen::VectorXd eigenvalues_imag;
        Eigen::MatrixXd eigenvectors_real;
        Eigen::MatrixXd eigenvectors_imag;
        bool bcomputevectors = true;
        bool bconv           = false;
        bool is_symmetric    = false;
    };
    
    /**
     * @brief Eigenvalue decomposition using Spectra (partial, k < n)
     */
    inline eigdecomp RcppbdEigen_spectra(const Eigen::MatrixXd& X, int k,
                                          const std::string& which = "LM",
                                          int ncv = 0,
                                          bool bcenter = false, bool bscale = false,
                                          double tol = 1e-10, int max_iter = 1000)
    {
        eigdecomp reteig;
        Eigen::MatrixXd nX;
        
        try {
            int n = X.rows();
            if (X.rows() != X.cols())
                throw std::runtime_error("Matrix must be square for eigendecomposition");
            
            std::tie(k, ncv) = validateSpectraParams(n, k, ncv);
            reteig.is_symmetric = isMatrixSymmetric(X);
            
            if (reteig.is_symmetric) {
                Eigen::MatrixXd Xcp = (bcenter || bscale)
                    ? RcppNormalize_Data(X, bcenter, bscale, false) : X;
                Spectra::DenseSymMatProd<double> op(Xcp);
                Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, k, ncv);
                eigs.init();
                (void)eigs.compute(getSymmetricSortRule(which), max_iter, tol);
                if (eigs.info() == Spectra::CompInfo::Successful) {
                    reteig.eigenvalues_real  = eigs.eigenvalues();
                    reteig.eigenvalues_imag  = Eigen::VectorXd::Zero(k);
                    reteig.eigenvectors_real = eigs.eigenvectors();
                    reteig.eigenvectors_imag = Eigen::MatrixXd::Zero(n, k);
                    reteig.bconv = true;
                }
            } else {
                Eigen::MatrixXd Xcp = (bcenter || bscale)
                    ? RcppNormalize_Data(X, bcenter, bscale, false) : X;
                Spectra::DenseGenMatProd<double> op(Xcp);
                Spectra::GenEigsSolver<Spectra::DenseGenMatProd<double>> eigs(op, k, ncv);
                eigs.init();
                (void)eigs.compute(getGeneralSortRule(which), max_iter, tol);
                if (eigs.info() == Spectra::CompInfo::Successful) {
                    Eigen::VectorXcd ev = eigs.eigenvalues();
                    Eigen::MatrixXcd evec = eigs.eigenvectors();
                    reteig.eigenvalues_real  = ev.real();
                    reteig.eigenvalues_imag  = ev.imag();
                    reteig.eigenvectors_real = evec.real();
                    reteig.eigenvectors_imag = evec.imag();
                    reteig.bconv = true;
                }
            }
            reteig.bcomputevectors = true;
            
        } catch(std::exception &ex) {
            throw std::runtime_error(std::string("C++ exception RcppbdEigen_spectra: ") + ex.what());
        } catch (...) {
            throw std::runtime_error("C++ exception RcppbdEigen_spectra (unknown reason)");
        }
        return reteig;
    }
    
    /**
     * @brief Block-wise eigendecomposition for large HDF5 matrices using Spectra
     *
     * @param dsA Input dataset (square matrix)
     * @param dsd Output dataset for eigenvalues (1 x k)
     * @param dsu Output dataset for eigenvectors (n x k), may be nullptr
     * @param k   Number of eigenvalues (must be < n, caller ensures this)
     */
    template <class T>
    inline void RcppbdEigen_hdf5_Block(T* dsA,
                                        BigDataStatMeth::hdf5Dataset* dsd,
                                        BigDataStatMeth::hdf5Dataset* dsu,
                                        int k, const std::string& which, int ncv,
                                        bool bcenter, bool bscale,
                                        double tol, int max_iter,
                                        bool compute_vectors,
                                        Rcpp::Nullable<int> threads = R_NilValue)
    {
        static_assert(std::is_same<T*, BigDataStatMeth::hdf5Dataset*>::value ||
                      std::is_same<T*, BigDataStatMeth::hdf5DatasetInternal*>::value,
                      "Error - type not allowed");
        try {
            std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
            
            hsize_t n = dsA->nrows();
            if (n != dsA->ncols())
                throw std::runtime_error("Matrix must be square for eigendecomposition");
            
            std::tie(k, ncv) = validateSpectraParams((int)n, k, ncv);
            
            std::vector<double> vdA(n * n);
            dsA->readDatasetBlock({0, 0}, {n, n}, stride, block, vdA.data());
            Eigen::MatrixXd X = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic,
                                                          Eigen::Dynamic, Eigen::RowMajor>>(
                vdA.data(), n, n);
            if (bcenter || bscale)
                X = RcppNormalize_Data(X, bcenter, bscale, false);
            
            eigdecomp reteig = RcppbdEigen_spectra(X, k, which, ncv, false, false, tol, max_iter);
            if (!reteig.bconv)
                throw std::runtime_error("Eigendecomposition failed to converge");
            
            int comp = (int)dsA->getCompressionLevel();
            
            dsd->inheritCompressionLevel(comp);
            dsd->createDataset(1, reteig.eigenvalues_real.size(), "real");
            dsd->writeDataset(Rcpp::wrap(reteig.eigenvalues_real));
            
            if (compute_vectors && reteig.bcomputevectors && dsu != nullptr) {
                dsu->inheritCompressionLevel(comp);
                dsu->createDataset(reteig.eigenvectors_real.rows(),
                                   reteig.eigenvectors_real.cols(), "real");
                dsu->writeDataset(Rcpp::wrap(reteig.eigenvectors_real));
            }
            
            if (!reteig.is_symmetric && reteig.eigenvalues_imag.cwiseAbs().maxCoeff() > 1e-14) {
                BigDataStatMeth::hdf5Dataset* dsd_imag = new BigDataStatMeth::hdf5Dataset(
                    dsA->getFullPath(), "EIGEN/" + dsA->getDatasetName(), "values_imag", true);
                dsd_imag->inheritCompressionLevel(comp);
                dsd_imag->createDataset(1, reteig.eigenvalues_imag.size(), "real");
                dsd_imag->writeDataset(Rcpp::wrap(reteig.eigenvalues_imag));
                delete dsd_imag;
                
                if (compute_vectors && reteig.bcomputevectors && dsu != nullptr &&
                    reteig.eigenvectors_imag.cwiseAbs().maxCoeff() > 1e-14) {
                    BigDataStatMeth::hdf5Dataset* dsu_imag = new BigDataStatMeth::hdf5Dataset(
                        dsA->getFullPath(), "EIGEN/" + dsA->getDatasetName(), "vectors_imag", true);
                    dsu_imag->inheritCompressionLevel(comp);
                    dsu_imag->createDataset(reteig.eigenvectors_imag.rows(),
                                            reteig.eigenvectors_imag.cols(), "real");
                    dsu_imag->writeDataset(Rcpp::wrap(reteig.eigenvectors_imag));
                    delete dsu_imag;
                }
            }
        } catch(H5::FileIException& error) {
            throw std::runtime_error("c++ exception RcppbdEigen_hdf5_Block (File IException)");
        } catch(H5::DataSetIException& error) {
            throw std::runtime_error("c++ exception RcppbdEigen_hdf5_Block (DataSet IException)");
        } catch(H5::GroupIException& error) {
            throw std::runtime_error("c++ exception RcppbdEigen_hdf5_Block (Group IException)");
        } catch(std::exception &ex) {
            throw std::runtime_error(std::string("C++ exception RcppbdEigen_hdf5_Block: ") + ex.what());
        } catch (...) {
            throw std::runtime_error("C++ exception RcppbdEigen_hdf5_Block (unknown reason)");
        }
    }
    
    /**
     * @brief Main eigenvalue decomposition for HDF5 matrices
     *
     * @param filename    Path to HDF5 file
     * @param strsubgroup Group path
     * @param strdataset  Dataset name (square matrix)
     * @param k           Number of eigenvalues. 0 or >= n -> full decomposition
     *                    via SelfAdjointEigenSolver. 0 < k < n -> Spectra (partial).
     * @param which       Which eigenvalues: LM, SM, LR, SR, LI, SI, LA, SA
     * @param ncv         Lanczos vectors (0 = auto)
     * @param bcenter     Center columns
     * @param bscale      Scale columns
     * @param tolerance   Convergence tolerance
     * @param max_iter    Maximum iterations
     * @param compute_vectors Compute eigenvectors
     * @param bforce      Overwrite existing results
     * @param threads     OpenMP threads
     *
     * @note Results written under "EIGEN/<dataset>/values" and "EIGEN/<dataset>/vectors"
     */
    inline void RcppbdEigen_hdf5(std::string filename,
                                  std::string strsubgroup,
                                  std::string strdataset,
                                  int k = 0,
                                  const std::string& which = "LM",
                                  int ncv = 0,
                                  bool bcenter = false, bool bscale = false,
                                  double tolerance = 1e-10, int max_iter = 1000,
                                  bool compute_vectors = true, bool bforce = false,
                                  Rcpp::Nullable<int> threads = R_NilValue)
    {
        try {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(nullptr);
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsd(nullptr);
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsu(nullptr);
            
            std::vector<hsize_t> stride = {1, 1}, block = {1, 1},
                                  offset = {0, 0}, count  = {0, 0};
            
            dsA.reset(new BigDataStatMeth::hdf5Dataset(filename, strsubgroup, strdataset, false));
            dsA->openDataset();
            
            if (dsA->getDatasetptr() != nullptr) {
                
                std::string stroutgroup = "EIGEN/" + strdataset;
                
                std::vector<hsize_t> dims_out = {dsA->nrows(), dsA->ncols()};
                count = {dims_out[0], dims_out[1]};
                
                if (dims_out[0] != dims_out[1])
                    throw std::runtime_error("Matrix must be square for eigendecomposition");
                
                hsize_t n = dims_out[0];
                
                // Detect full decomposition BEFORE validateSpectraParams clamps k.
                // k <= 0 means "all" (R layer converts 0 -> n but guard here as well).
                bool full_decomp = (k <= 0 || k >= (int)n);
                
                // For Spectra path: validate params now. For full path: skip (not needed).
                if (!full_decomp)
                    std::tie(k, ncv) = validateSpectraParams((int)n, k, ncv);
                
                dsd.reset(new BigDataStatMeth::hdf5Dataset(filename, stroutgroup, "values", bforce));
                if (compute_vectors)
                    dsu.reset(new BigDataStatMeth::hdf5Dataset(filename, stroutgroup, "vectors", bforce));
                
                // ── Small matrices: load entirely into memory ──────────────────────
                if (n * n < MAXELEMSINBLOCK / 20) {
                    
                    std::vector<double> vdA(count[0] * count[1]);
                    dsA->readDatasetBlock({offset[0], offset[1]}, {count[0], count[1]},
                                          stride, block, vdA.data());
                    Eigen::MatrixXd X = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic,
                                                                   Eigen::Dynamic,
                                                                   Eigen::RowMajor>>(
                        vdA.data(), count[0], count[1]);
                    
                    eigdecomp reteig;
                    
                    if (full_decomp) {
                        // All eigenvalues: use direct solver (Spectra requires k < n)
                        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(X);
                        if (solver.info() != Eigen::Success)
                            throw std::runtime_error("SelfAdjointEigenSolver failed to converge");
                        
                        int nf = (int)n;
                        reteig.eigenvalues_real.resize(nf);
                        for (int i = 0; i < nf; ++i)
                            reteig.eigenvalues_real(i) = solver.eigenvalues()(nf - 1 - i);
                        
                        if (compute_vectors) {
                            reteig.eigenvectors_real.resize(nf, nf);
                            for (int j = 0; j < nf; ++j)
                                reteig.eigenvectors_real.col(j) =
                                    solver.eigenvectors().col(nf - 1 - j);
                        }
                        reteig.bconv           = true;
                        reteig.is_symmetric    = true;
                        reteig.bcomputevectors = compute_vectors;
                        
                    } else {
                        reteig = RcppbdEigen_spectra(X, k, which, ncv,
                                                      bcenter, bscale, tolerance, max_iter);
                        if (!reteig.bconv)
                            throw std::runtime_error("Eigendecomposition failed to converge");
                    }
                    
                    // Write eigenvalues (1 x k dataset)
                    dsd->inheritCompressionLevel(dsA->getCompressionLevel());
                    dsd->createDataset(1, reteig.eigenvalues_real.size(), "real");
                    dsd->writeDataset(Rcpp::wrap(reteig.eigenvalues_real));
                    
                    // Write eigenvectors (n x k dataset)
                    if (compute_vectors && reteig.bcomputevectors) {
                        dsu->inheritCompressionLevel(dsA->getCompressionLevel());
                        dsu->createDataset(reteig.eigenvectors_real.rows(),
                                           reteig.eigenvectors_real.cols(), "real");
                        dsu->writeDataset(Rcpp::wrap(reteig.eigenvectors_real));
                    }
                    
                    // Imaginary parts (non-symmetric matrices only)
                    if (!reteig.is_symmetric &&
                        reteig.eigenvalues_imag.cwiseAbs().maxCoeff() > 1e-14) {
                        
                        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsd_imag(
                            new BigDataStatMeth::hdf5Dataset(filename, stroutgroup,
                                                              "values_imag", bforce));
                        dsd_imag->inheritCompressionLevel(dsA->getCompressionLevel());
                        dsd_imag->createDataset(1, reteig.eigenvalues_imag.size(), "real");
                        dsd_imag->writeDataset(Rcpp::wrap(reteig.eigenvalues_imag));
                        
                        if (compute_vectors && reteig.bcomputevectors &&
                            reteig.eigenvectors_imag.cwiseAbs().maxCoeff() > 1e-14) {
                            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsu_imag(
                                new BigDataStatMeth::hdf5Dataset(filename, stroutgroup,
                                                                  "vectors_imag", bforce));
                            dsu_imag->inheritCompressionLevel(dsA->getCompressionLevel());
                            dsu_imag->createDataset(reteig.eigenvectors_imag.rows(),
                                                    reteig.eigenvectors_imag.cols(), "real");
                            dsu_imag->writeDataset(Rcpp::wrap(reteig.eigenvectors_imag));
                        }
                    }
                    
                } else {
                    // ── Large matrices: block-wise path ───────────────────────────
                    // For full_decomp on large matrices: use n-1 (best we can do with Spectra)
                    int k_block = full_decomp ? (int)(n - 1) : k;
                    RcppbdEigen_hdf5_Block(dsA.get(), dsd.get(),
                                            compute_vectors ? dsu.get() : nullptr,
                                            k_block, which, ncv,
                                            bcenter, bscale, tolerance, max_iter,
                                            compute_vectors, threads);
                }
            }
            
        } catch(H5::FileIException& error) {
            throw std::runtime_error("c++ exception RcppbdEigen_hdf5 (File IException)");
        } catch(H5::DataSetIException& error) {
            throw std::runtime_error(
                std::string("c++ exception RcppbdEigen_hdf5 (DataSet IException): ")
                + error.getDetailMsg());
        } catch(std::exception &ex) {
            throw std::runtime_error(
                std::string("c++ exception RcppbdEigen_hdf5: ") + ex.what());
        } catch (...) {
            throw std::runtime_error("C++ exception RcppbdEigen_hdf5 (unknown reason)");
        }
    }

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_HDF5_MATRIXEIGEN_HPP
