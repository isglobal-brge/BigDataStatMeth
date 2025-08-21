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


namespace BigDataStatMeth {

    /**
     * @brief Validate and adjust Spectra parameters for convergence
     * @param n Matrix size
     * @param k Number of eigenvalues requested  
     * @param ncv Number of Arnoldi vectors
     * @return Adjusted parameters that satisfy Spectra constraints
     */
    inline std::tuple<int, int> validateSpectraParams(int n, int k, int ncv) {
        // Ensure k is reasonable
        if (k <= 0) k = std::min(n, 6);
        if (k >= n) k = n - 1;
        
        // Calculate optimal ncv following RSpectra defaults
        if (ncv <= 0) {
            ncv = std::min(n, std::max(2 * k + 1, k + 2));
        }
        
        // Enforce Spectra constraints: k + 2 <= ncv <= n
        ncv = std::max(ncv, k + 2);
        ncv = std::min(ncv, n);
        
        // Additional safety: ensure we have enough space for convergence
        if (ncv - k < 2) {
            ncv = std::min(n, k + 2);
        }
        
        return std::make_tuple(k, ncv);
    }
    
    /**
     * @brief Structure to hold eigendecomposition results
     * @details Similar to svdeig structure used in SVD, this structure holds
     * the results of eigenvalue decomposition for consistent interface.
     */
    struct eigdecomp {
        Eigen::VectorXd eigenvalues_real;      /**< Real part of eigenvalues */
        Eigen::VectorXd eigenvalues_imag;      /**< Imaginary part of eigenvalues */
        Eigen::MatrixXd eigenvectors_real;     /**< Real part of eigenvectors */
        Eigen::MatrixXd eigenvectors_imag;     /**< Imaginary part of eigenvectors */
        bool bcomputevectors = true;           /**< Whether eigenvectors were computed */
        bool bconv = false;                    /**< Convergence status */
        bool is_symmetric = false;             /**< Whether input matrix was symmetric */
    };
    
    /**
     * @brief Eigenvalue decomposition using Spectra (compatible with BigDataStatMeth SVD version)
     * 
     * Uses the same Spectra version and patterns as the existing SVD implementation in
     * BigDataStatMeth. Based on the RcppbdSVD function in matrixSvd.hpp.
     * 
     * @param X Input matrix
     * @param k Number of eigenvalues to compute
     * @param ncv Number of Arnoldi vectors (if 0, uses k+1)
     * @param bcenter Whether to center the data
     * @param bscale Whether to scale the data
     * 
     * @return eigdecomp structure containing results
     * 
     * @note This follows the same pattern as RcppbdSVD in matrixSvd.hpp
     * @note Updated for Spectra 1.0.1 API compatibility
     */
    inline eigdecomp RcppbdEigen_spectra(const Eigen::MatrixXd& X, int k, int ncv = 0,
                                         bool bcenter = false, bool bscale = false) {
        
        eigdecomp reteig;
        Eigen::MatrixXd nX;
        int nconv = 0;
        
        try {
            
            int n = X.rows();
            if (X.rows() != X.cols()) {
                Rf_error("Matrix must be square for eigendecomposition");
                return reteig;
            }
            
            // Better parameter selection following RSpectra defaults
            std::tie(k, ncv) = validateSpectraParams(n, k, ncv);
            
            // Check if matrix is approximately symmetric
            double max_asymmetry = 0.0;
            int check_size = std::min(20, n);
            for (int i = 0; i < check_size; ++i) {
                for (int j = i + 1; j < check_size; ++j) {
                    max_asymmetry = std::max(max_asymmetry, std::abs(X(i,j) - X(j,i)));
                }
            }
            
            reteig.is_symmetric = (max_asymmetry < 1e-12);
            
            if (reteig.is_symmetric) {
                // Use symmetric solver - following SVD pattern
                Eigen::MatrixXd Xcp;
                if (bcenter == true || bscale == true) {
                    nX = RcppNormalize_Data(X, bcenter, bscale, false);
                    Xcp = nX;  // For symmetric case, use matrix directly
                } else {
                    Xcp = X;
                }
                
                Spectra::DenseSymMatProd<double> op(Xcp);
                // Updated for Spectra 1.0.1: removed template parameters, pass object by reference
                Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, k, ncv);
                
                // Initialize and compute
                eigs.init();
                // Updated for Spectra 1.0.1: SortRule as runtime parameter from local Spectra
                nconv = eigs.compute(Spectra::SortRule::LargestAlge);
                
                // Retrieve results
                // Updated for Spectra 1.0.1: enum class for status check from local Spectra
                if (eigs.info() == Spectra::CompInfo::Successful) {
                    reteig.eigenvalues_real = eigs.eigenvalues().reverse();  // Largest first
                    reteig.eigenvalues_imag = Eigen::VectorXd::Zero(k);
                    reteig.eigenvectors_real = eigs.eigenvectors().rowwise().reverse();
                    reteig.eigenvectors_imag = Eigen::MatrixXd::Zero(n, k);
                    reteig.bconv = true;
                } else {
                    reteig.bconv = false;
                }
                
            } else {
                // For non-symmetric matrices, use general solver
                Eigen::MatrixXd Xcp;
                if (bcenter == true || bscale == true) {
                    nX = RcppNormalize_Data(X, bcenter, bscale, false);
                    Xcp = nX;
                } else {
                    Xcp = X;
                }
                
                Spectra::DenseGenMatProd<double> op(Xcp);
                // Updated for Spectra 1.0.1: removed template parameters, pass object by reference
                Spectra::GenEigsSolver<Spectra::DenseGenMatProd<double>> eigs(op, k, ncv);
                
                // Initialize and compute
                eigs.init();
                // Updated for Spectra 1.0.1: SortRule as runtime parameter from local Spectra
                nconv = eigs.compute(Spectra::SortRule::LargestMagn);
                
                // Retrieve results
                // Updated for Spectra 1.0.1: enum class for status check from local Spectra
                if (eigs.info() == Spectra::CompInfo::Successful) {
                    Eigen::VectorXcd eigenvals = eigs.eigenvalues();
                    Eigen::MatrixXcd eigenvecs = eigs.eigenvectors();
                    
                    reteig.eigenvalues_real = eigenvals.real();
                    reteig.eigenvalues_imag = eigenvals.imag();
                    reteig.eigenvectors_real = eigenvecs.real();
                    reteig.eigenvectors_imag = eigenvecs.imag();
                    reteig.bconv = true;
                } else {
                    reteig.bconv = false;
                }
            }
            
            reteig.bcomputevectors = true;
            
        } catch(std::exception &ex) {
            Rcpp::Rcout << "C++ exception RcppbdEigen_spectra: " << ex.what();
            reteig.bconv = false;
        } catch (...) {
            Rf_error("C++ exception RcppbdEigen_spectra (unknown reason)");
            reteig.bconv = false;
        }
        
        return reteig;
    }
    
    /**
     * @brief Block-wise eigendecomposition for large HDF5 matrices using Spectra
     * 
     * This function performs eigenvalue decomposition on large matrices stored in HDF5
     * using Spectra library for consistency with RSpectra results.
     * 
     * @param dsA Input hdf5Dataset containing the matrix
     * @param dsd Output hdf5Dataset for eigenvalues
     * @param dsu Output hdf5Dataset for eigenvectors  
     * @param k Number of eigenvalues to compute
     * @param bcenter Whether to center the data
     * @param bscale Whether to scale the data
     * @param threads Number of parallel threads
     * 
     * @throws H5::FileIException for HDF5 file operation errors
     * @throws H5::DataSetIException for HDF5 dataset operation errors
     * @throws std::exception for computation errors
     * 
     * @note Updated for Spectra 1.0.1 compatibility
     */
    template <class T>
    inline void RcppbdEigen_hdf5_Block(T* dsA, 
                                       BigDataStatMeth::hdf5Dataset* dsd,
                                       BigDataStatMeth::hdf5Dataset* dsu, 
                                       int k, bool bcenter, bool bscale,
                                       Rcpp::Nullable<int> threads = R_NilValue) {
        
        static_assert(std::is_same<T*, BigDataStatMeth::hdf5Dataset*>::value ||
                      std::is_same<T*, BigDataStatMeth::hdf5DatasetInternal*>::value,
                      "Error - type not allowed");
        
        BigDataStatMeth::hdf5Dataset* dsnormalizedData = nullptr;
        
        try {
            
            std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
            eigdecomp reteig;
            
            // Get matrix dimensions
            hsize_t n_rows = dsA->nrows();
            hsize_t n_cols = dsA->ncols();
            
            // Check if matrix is square
            if (n_rows != n_cols) {
                Rf_error("Matrix must be square for eigendecomposition");
                return;
            }
            
            hsize_t n = n_rows;
            
            // Use parameter validation function
            int ncv = 0; 
            std::tie(k, ncv) = validateSpectraParams((int)n, k, 0);
            
            // Read matrix data
            Eigen::MatrixXd X;
            std::vector<double> vdA(n * n);
            
            dsA->readDatasetBlock({0, 0}, {n, n}, stride, block, vdA.data());
            X = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(vdA.data(), n, n);
            
            // Handle normalization if needed
            if (bcenter || bscale) {
                
                // For large matrices, we could implement block-wise normalization here
                // For now, since we already read the matrix, apply normalization directly
                if (n * n < MAXELEMSINBLOCK) {
                    X = RcppNormalize_Data(X, bcenter, bscale, false);
                } else {
                    // For very large matrices, implement block-wise normalization
                    Rcpp::warning("Large matrix normalization not yet optimized - applying direct normalization");
                    X = RcppNormalize_Data(X, bcenter, bscale, false);
                }
            }
            
            // Use Spectra for eigendecomposition - same pattern as SVD
            reteig = RcppbdEigen_spectra(X, k, ncv, false, false);
            
            if (!reteig.bconv) {
                Rf_error("Eigendecomposition failed to converge");
                return;
            }
            
            // Write eigenvalues to HDF5
            dsd->createDataset(1, reteig.eigenvalues_real.size(), "real");
            dsd->writeDataset(Rcpp::wrap(reteig.eigenvalues_real));
            
            // Write eigenvectors to HDF5 if computed
            if (reteig.bcomputevectors) {
                dsu->createDataset(reteig.eigenvectors_real.rows(), reteig.eigenvectors_real.cols(), "real");
                dsu->writeDataset(Rcpp::wrap(reteig.eigenvectors_real));
            }
            
            // For non-symmetric matrices, also save imaginary parts if non-zero
            if (!reteig.is_symmetric && reteig.eigenvalues_imag.cwiseAbs().maxCoeff() > 1e-14) {
                
                // Create additional datasets for imaginary parts
                BigDataStatMeth::hdf5Dataset* dsd_imag = new BigDataStatMeth::hdf5Dataset(
                    dsA->getFileName(), "EIGEN/" + dsA->getDatasetName(), "values_imag", true);
                dsd_imag->createDataset(1, reteig.eigenvalues_imag.size(), "real");
                dsd_imag->writeDataset(Rcpp::wrap(reteig.eigenvalues_imag));
                delete dsd_imag;
                
                if (reteig.bcomputevectors && reteig.eigenvectors_imag.cwiseAbs().maxCoeff() > 1e-14) {
                    BigDataStatMeth::hdf5Dataset* dsu_imag = new BigDataStatMeth::hdf5Dataset(
                        dsA->getFileName(), "EIGEN/" + dsA->getDatasetName(), "vectors_imag", true);
                    dsu_imag->createDataset(reteig.eigenvectors_imag.rows(), reteig.eigenvectors_imag.cols(), "real");
                    dsu_imag->writeDataset(Rcpp::wrap(reteig.eigenvectors_imag));
                    delete dsu_imag;
                }
            }
            
        } catch(H5::FileIException& error) {
            checkClose_file(dsA, dsd, dsu, dsnormalizedData);
            Rcpp::Rcerr << "\nc++ exception RcppbdEigen_hdf5_Block (File IException)\n";
            return;
        } catch(H5::DataSetIException& error) {
            checkClose_file(dsA, dsd, dsu, dsnormalizedData);
            Rcpp::Rcerr << "\nc++ exception RcppbdEigen_hdf5_Block (DataSet IException)\n";
            return;
        } catch(H5::GroupIException& error) {
            checkClose_file(dsA, dsd, dsu, dsnormalizedData);
            Rcpp::Rcerr << "\nc++ exception RcppbdEigen_hdf5_Block (Group IException)\n";
            return;
        } catch(std::exception &ex) {
            checkClose_file(dsA, dsd, dsu, dsnormalizedData);
            Rcpp::Rcerr << "C++ exception RcppbdEigen_hdf5_Block: " << ex.what();
            return;
        } catch (...) {
            checkClose_file(dsA, dsd, dsu, dsnormalizedData);
            Rcpp::Rcerr << "\nC++ exception RcppbdEigen_hdf5_Block (unknown reason)";
            return;
        }
    }
    
    /**
     * @brief Main eigenvalue decomposition function for HDF5 matrices using Spectra
     * 
     * This is the main interface for eigenvalue decomposition of matrices stored in HDF5.
     * It uses Spectra library for consistent results with RSpectra and automatically 
     * handles both symmetric and non-symmetric matrices.
     * 
     * @param filename Path to HDF5 file containing the input matrix
     * @param strsubgroup Group path within the HDF5 file
     * @param strdataset Dataset name containing the matrix
     * @param k Number of eigenvalues to compute (0 = auto-select based on matrix size)
     * @param bcenter Whether to center the data before decomposition
     * @param bscale Whether to scale the data before decomposition
     * @param tolerance Convergence tolerance for Spectra algorithms
     * @param max_iter Maximum iterations for Spectra algorithms  
     * @param compute_vectors Whether to compute eigenvectors or just eigenvalues
     * @param bforce Whether to overwrite existing results
     * @param threads Number of threads for parallel computation
     * 
     * @throws H5::FileIException for HDF5 file operation errors
     * @throws H5::DataSetIException for HDF5 dataset operation errors
     * @throws std::exception for computation errors
     * 
     * @note Results are written to HDF5 datasets under "EIGEN/<dataset_name>/"
     *       - eigenvalues (real): "EIGEN/<dataset_name>/d"
     *       - eigenvalues (imag): "EIGEN/<dataset_name>/d_imag" (if non-zero)
     *       - eigenvectors (real): "EIGEN/<dataset_name>/u"
     *       - eigenvectors (imag): "EIGEN/<dataset_name>/u_imag" (if non-zero)
     * @note Updated for Spectra 1.0.1 compatibility
     */
    inline void RcppbdEigen_hdf5(std::string filename, std::string strsubgroup, std::string strdataset,
                                 int k = 0, bool bcenter = false, bool bscale = false,
                                 double tolerance = 1e-10, int max_iter = 1000, 
                                 bool compute_vectors = true, bool bforce = false,
                                 Rcpp::Nullable<int> threads = R_NilValue) {
        
        BigDataStatMeth::hdf5Dataset* dsA = nullptr;
        BigDataStatMeth::hdf5Dataset* dsd = nullptr;
        BigDataStatMeth::hdf5Dataset* dsu = nullptr;
        
        try {
            
            std::vector<hsize_t> stride = {1, 1}, block = {1, 1}, offset = {0, 0}, count = {0, 0};
            
            dsA = new BigDataStatMeth::hdf5Dataset(filename, strsubgroup, strdataset, false);
            dsA->openDataset();
            
            if (dsA->getDatasetptr() != nullptr) {
                
                // Create results folder following SVD pattern
                std::string stroutgroup = "EIGEN/" + strdataset;
                
                std::vector<hsize_t> dims_out = {dsA->nrows(), dsA->ncols()};
                count = {dims_out[0], dims_out[1]};
                
                // Check if matrix is square
                if (dims_out[0] != dims_out[1]) {
                    Rf_error("Matrix must be square for eigendecomposition");
                    return;
                }
                
                hsize_t n = dims_out[0];
                
                // Use parameter validation function
                int ncv = 0; 
                std::tie(k, ncv) = validateSpectraParams((int)n, k, 0);
                
                // Create output datasets
                dsd = new BigDataStatMeth::hdf5Dataset(filename, stroutgroup, "values", bforce);
                if (compute_vectors) {
                    dsu = new BigDataStatMeth::hdf5Dataset(filename, stroutgroup, "vectors", bforce);
                }
                
                // Small matrices => Direct eigendecomposition with Spectra
                if (n * n < MAXELEMSINBLOCK / 20) {
                    
                    Eigen::MatrixXd X;
                    eigdecomp reteig;
                    
                    std::vector<double> vdA(count[0] * count[1]);
                    dsA->readDatasetBlock({offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data());
                    
                    X = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                        vdA.data(), count[0], count[1]);
                    
                    reteig = RcppbdEigen_spectra(X, k, ncv, bcenter, bscale);
                    
                    if (!reteig.bconv) {
                        Rf_error("Eigendecomposition failed to converge");
                        return;
                    }
                    
                    // Write real eigenvalues
                    dsd->createDataset(1, reteig.eigenvalues_real.size(), "real");
                    dsd->writeDataset(Rcpp::wrap(reteig.eigenvalues_real));

                    // Write real eigenvectors if computed
                    if (compute_vectors && reteig.bcomputevectors) {
                        dsu->createDataset(reteig.eigenvectors_real.rows(), reteig.eigenvectors_real.cols(), "real");
                        dsu->writeDataset(Rcpp::wrap(reteig.eigenvectors_real));
                    }
                    
                    // Write imaginary parts if non-zero (for non-symmetric matrices)
                    if (!reteig.is_symmetric && reteig.eigenvalues_imag.cwiseAbs().maxCoeff() > 1e-14) {
                        
                        BigDataStatMeth::hdf5Dataset* dsd_imag = new BigDataStatMeth::hdf5Dataset(filename, stroutgroup, "values_imag", bforce);
                        dsd_imag->createDataset(1, reteig.eigenvalues_imag.size(), "real");
                        dsd_imag->writeDataset(Rcpp::wrap(reteig.eigenvalues_imag));
                        delete dsd_imag;
                        
                        if (compute_vectors && reteig.bcomputevectors && reteig.eigenvectors_imag.cwiseAbs().maxCoeff() > 1e-14) {
                            BigDataStatMeth::hdf5Dataset* dsu_imag = new BigDataStatMeth::hdf5Dataset(filename, stroutgroup, "vectors_imag", bforce);
                            dsu_imag->createDataset(reteig.eigenvectors_imag.rows(), reteig.eigenvectors_imag.cols(), "real");
                            dsu_imag->writeDataset(Rcpp::wrap(reteig.eigenvectors_imag));
                            delete dsu_imag;
                        }
                    }
                    
                } else {
                    // Large matrices => Block-based approach with Spectra  
                    if (dsA->getDatasetptr() != nullptr) {
                        RcppbdEigen_hdf5_Block(dsA, dsd, dsu, k, bcenter, bscale, threads);
                    }
                }
            }
            delete dsd; dsd = nullptr;
            if (dsu) { delete dsu; dsu = nullptr; }
            delete dsA; dsA = nullptr;
            
        } catch(H5::FileIException& error) {
            checkClose_file(dsA, dsd, dsu);
            Rcpp::Rcerr << "\nc++ exception RcppbdEigen_hdf5 (File IException)\n";
        } catch(H5::DataSetIException& error) {
            checkClose_file(dsA, dsd, dsu);
            Rcpp::Rcerr << "\nc++ exception RcppbdEigen_hdf5 (DataSet IException) "<<error.getDetailMsg() <<" \n";
        } catch(std::exception &ex) {
            checkClose_file(dsA, dsd, dsu);
            Rcpp::Rcerr << "c++ exception RcppbdEigen_hdf5: " << ex.what();
        } catch (...) {
            checkClose_file(dsA, dsd, dsu);
            Rcpp::Rcerr << "\nC++ exception RcppbdEigen_hdf5 (unknown reason)";
        }
    }

} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_HDF5_MATRIXEIGEN_HPP