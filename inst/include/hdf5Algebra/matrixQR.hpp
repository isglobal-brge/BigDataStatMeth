/**
 * @file matrixQR.hpp
 * @brief QR decomposition for HDF5 matrices
 * @note 2026-03-07 Output datasets now inherit compression level from input datasets
 *         via setCompressionLevel() called before every createDataset() invocation.
 * @details This header file provides implementations for computing the QR
 * decomposition of matrices stored in HDF5 format. The implementation includes:
 * 
 * Key features:
 * - Full and thin QR decomposition
 * - Memory-efficient algorithms
 * - Parallel processing support
 * - Rank computation
 * - Householder transformations
 * 
 * Supported operations:
 * - Full QR decomposition
 * - Thin QR decomposition
 * - Rank-revealing QR
 * - In-memory computation
 * - HDF5 storage computation
 * 
 * Performance features:
 * - Eigen optimizations
 * - Cache-friendly algorithms
 * - Multi-threaded processing (via Eigen::setNbThreads for LAPACK method;
 *   via OpenMP block parallelism for TSQR method)
 * - I/O optimization
 * - Memory management
 * 
 * The implementation uses:
 * - Householder QR algorithm  (method = "lapack")
 * - TSQR (Tall-Skinny QR)     (method = "tsqr")
 * - Block-based computation
 * - HDF5 chunked storage
 * - Parallel processing
 * - Vectorized operations
 *
 * Algorithm selection (method = "auto"):
 *   - TSQR  when m > 5*n  AND  m > 1000  (tall-skinny, benefits from parallelism)
 *   - LAPACK otherwise                     (square / small / wide matrices)
 *
 * @note 20260304: Added TSQR algorithm and method parameter.
 *   - "lapack"  : existing Eigen HouseholderQR, single-threaded Householder.
 *                 Eigen internal BLAS threads controlled via Eigen::setNbThreads.
 *   - "tsqr"    : parallel row-block QR + R-factor reduction.
 *                 Q produced is always thin (m×k); for thin=FALSE a basis
 *                 completion step produces the full m×m orthogonal factor
 *                 (expensive: O(m² k), use only when needed).
 *   - "auto"    : selects TSQR for tall-skinny matrices, LAPACK otherwise.
 *   The R factor is identical (up to sign) for both methods when full rank.
 */

#ifndef BIGDATASTATMETH_HDF5_MATRIXQR_HPP
#define BIGDATASTATMETH_HDF5_MATRIXQR_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

namespace BigDataStatMeth {

// ── Forward declarations ─────────────────────────────────────────────────────

/**
 * @brief Template function for QR decomposition
 * @details Computes the QR decomposition of a matrix using Householder
 * transformations. Supports both full and thin decomposition.
 * 
 * @tparam M Matrix type (Eigen::MatrixXd or mapped matrix)
 * @param X Input matrix to decompose
 * @param bthin Whether to compute thin QR
 * @return Structure containing Q and R matrices
 */
template< typename M> inline strQR RcppQR ( M X, bool bthin );

/**
 * @brief QR decomposition for HDF5 matrices (Householder / LAPACK backend)
 * @details Computes the QR decomposition of a matrix stored in HDF5 format.
 * Supports both full and thin decomposition with optional Eigen thread control.
 * 
 * @param dsA        Input matrix dataset
 * @param dsQ        Output Q matrix dataset
 * @param dsR        Output R matrix dataset
 * @param bthin      Whether to compute thin QR
 * @param block_size Block size for processing (optional, currently unused for LAPACK)
 * @param threads    Number of Eigen internal threads (-1 = auto)
 */
inline void RcppQRHdf5( BigDataStatMeth::hdf5Dataset* dsA,
                        BigDataStatMeth::hdf5Dataset* dsQ,
                        BigDataStatMeth::hdf5Dataset* dsR, 
                        bool bthin,
                        Rcpp::Nullable<int> block_size,
                        Rcpp::Nullable<int> threads );

/**
 * @brief TSQR (Tall-Skinny QR) decomposition for HDF5 matrices
 * @details Parallel QR for tall-skinny matrices (m >> n).
 * Divides A into row blocks, factorises each block in parallel via OpenMP,
 * then reduces the local R factors with a sequential QR on their stack.
 * Assembles the global thin Q from local Q factors and the reduction Q.
 *
 * For thin=FALSE the thin Q is completed to a full m×m orthogonal factor
 * via a Householder basis-completion step (O(m² k), expensive for large m).
 *
 * @param dsA        Input matrix dataset (m × n, m >> n recommended)
 * @param dsQ        Output Q matrix dataset
 * @param dsR        Output R matrix dataset
 * @param bthin      If true, Q is m×k (thin); if false, Q is m×m (full)
 * @param block_size Minimum rows per block (NULL = auto based on threads)
 * @param threads    OpenMP threads for the local QR phase (NULL = auto)
 *
 * @note The R factor from TSQR is mathematically equivalent to that from
 *   Householder QR (up to sign flips on rows). The Q factor is generally
 *   different but equally valid orthogonal matrix.
 */
inline void RcppTSQRHdf5( BigDataStatMeth::hdf5Dataset* dsA,
                           BigDataStatMeth::hdf5Dataset* dsQ,
                           BigDataStatMeth::hdf5Dataset* dsR,
                           bool bthin,
                           Rcpp::Nullable<int> block_size,
                           Rcpp::Nullable<int> threads );

/**
 * @brief Unified QR entry point — dispatches to LAPACK or TSQR based on method.
 * @details Wrapper that selects the appropriate QR algorithm:
 *   - "lapack" : Eigen HouseholderQR (existing implementation)
 *   - "tsqr"   : Parallel TSQR (new, efficient for tall-skinny)
 *   - "auto"   : TSQR if m > 5*n AND m > 1000, otherwise LAPACK
 *
 * @param dsA        Input matrix dataset
 * @param dsQ        Output Q matrix dataset
 * @param dsR        Output R matrix dataset
 * @param bthin      Thin (true) or full (false) QR
 * @param block_size Block size hint (NULL = auto)
 * @param threads    Thread count (NULL = auto)
 * @param method     "auto", "lapack", or "tsqr"
 */
inline void RcppQRHdf5_dispatch( BigDataStatMeth::hdf5Dataset* dsA,
                                  BigDataStatMeth::hdf5Dataset* dsQ,
                                  BigDataStatMeth::hdf5Dataset* dsR,
                                  bool bthin,
                                  Rcpp::Nullable<int> block_size,
                                  Rcpp::Nullable<int> threads,
                                  const std::string& method );


// ── RcppQR (in-memory template) ──────────────────────────────────────────────

template< typename M>
inline strQR RcppQR ( M X, bool bthin )
{
    
    static_assert(std::is_same<M, Eigen::MatrixXd >::value || 
                  std::is_same<M, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                  std::is_same<M, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
                  "Error - type not allowed");
    
    Eigen::MatrixXd A = X;
    int m = A.rows(), n = A.cols();
    int irank;
    strQR vQR;
    
    Eigen::MatrixXd R;
    Eigen::MatrixXd Q;
    
    Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    
    qr.compute(A);
    irank = lu_decomp.rank();
    
    if (irank == m + 1 || irank == n + 1 )
    {
        vQR.R = qr.matrixQR().template triangularView<Eigen::Upper>();
    } else {
        vQR.R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>(); 
    }
    
    if (bthin == false)
    {
        vQR.Q =  qr.householderQ();       // Full decomposition
    } else {
        
        vQR.Q = Eigen::MatrixXd::Identity(m,n);
        vQR.Q = qr.householderQ() * vQR.Q;    // Thin decomposition
    }
    
    return(vQR);
}



// extern strQR RcppQR( Eigen::MatrixXd & A, bool bthin)
// {
//     
//     int m = A.rows(), n = A.cols();
//     int irank;
//     strQR vQR;
//     
//     Eigen::MatrixXd R;
//     Eigen::MatrixXd Q;
//     
//     Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
//     Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
//     
//     qr.compute(A);
//     irank = lu_decomp.rank();
//     
//     if (irank == m + 1 || irank == n + 1 )
//     {
//         vQR.R = qr.matrixQR().template triangularView<Eigen::Upper>();
//     } else {
//         vQR.R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>(); 
//     }
//     
//     if (bthin == false)
//     {
//         vQR.Q =  qr.householderQ();       // Full decomposition
//     } else {
//         
//         vQR.Q = Eigen::MatrixXd::Identity(m,n);
//         vQR.Q = qr.householderQ() * vQR.Q;    // Thin decomposition
//     }
//     
//     return(vQR);
//     
// }



// ── RcppQRHdf5 (LAPACK / Householder backend) ────────────────────────────────

inline void RcppQRHdf5( BigDataStatMeth::hdf5Dataset* dsA, 
                        BigDataStatMeth::hdf5Dataset* dsQ, 
                        BigDataStatMeth::hdf5Dataset* dsR, 
                        bool bthin,
                        Rcpp::Nullable<int> block_size = R_NilValue,
                        Rcpp::Nullable<int> threads    = R_NilValue )
{
    
    try {
        // ── Optionally set Eigen BLAS thread count ────────────────────────────
        // Eigen::setNbThreads uses its own internal pool (separate from OpenMP),
        // so it is safe to call here; it respects the getDTthreads limit via the
        // caller capping threads before passing them in.
        if (threads.isNotNull()) {
            int nth = std::max(1, Rcpp::as<int>(threads));
            Eigen::setNbThreads(nth);
        }

        // ── Read A from HDF5 (row-major → transpose to get R column-major) ───
        std::vector<hsize_t> offset = {0,0},
            count  = {dsA->nrows(), dsA->ncols()},
            stride = {1,1},
            block  = {1,1};
        
        std::vector<double> vdA( count[0] * count[1] );
        dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]},
                               stride, block, vdA.data() );
        
        // count[0]=ncols_R, count[1]=nrows_R  (HDF5 is transposed vs R)
        Eigen::MatrixXd A = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdA.data(), count[0], count[1]);
        A.transposeInPlace();   // now A is m×n  (R dimensions)
        
        const int m = static_cast<int>(A.rows());
        const int n = static_cast<int>(A.cols());
        const int k = std::min(m, n);   // number of significant components
        
        // ── QR decomposition ─────────────────────────────────────────────────
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
        
        if (!bthin) {
            // ── Full QR: Q is m×m, R is m×n (upper trapezoidal) ─────────────
            // Q
            Eigen::MatrixXd Qfull = qr.householderQ();
            dsQ->inheritCompressionLevel(dsA->getCompressionLevel());
            dsQ->createDataset( Qfull.rows(), Qfull.cols(), "real" );
            dsQ->writeDataset( Rcpp::wrap(Qfull) );
            
            // R: m×n upper trapezoidal (zeros below the first k rows diagonal)
            Eigen::MatrixXd Rfull = Eigen::MatrixXd::Zero(m, n);
            Rfull.topRows(k) = qr.matrixQR().topRows(k)
                 .template triangularView<Eigen::Upper>();
            dsR->inheritCompressionLevel(dsA->getCompressionLevel());
            dsR->createDataset( Rfull.rows(), Rfull.cols(), "real" );
            dsR->writeDataset( Rcpp::wrap(Rfull) );
            
        } else {
            // ── Thin QR: Q is m×k, R is k×n ─────────────────────────────────
            // Q thin
            Eigen::MatrixXd Qthin = qr.householderQ()
            * Eigen::MatrixXd::Identity(m, k);
            dsQ->inheritCompressionLevel(dsA->getCompressionLevel());
            dsQ->createDataset( Qthin.rows(), Qthin.cols(), "real" );
            dsQ->writeDataset( Rcpp::wrap(Qthin) );
            
            // R thin: k×n upper triangular
            Eigen::MatrixXd Rthin = qr.matrixQR().topRows(k)
                                      .template triangularView<Eigen::Upper>();
            dsR->inheritCompressionLevel(dsA->getCompressionLevel());
            dsR->createDataset( Rthin.rows(), Rthin.cols(), "real" );
            dsR->writeDataset( Rcpp::wrap(Rthin) );
        }
        
    } catch( H5::FileIException& error ) {
        throw std::runtime_error("c++ exception RcppQRHdf5 (File IException)");
    } catch( H5::GroupIException & error ) {
        throw std::runtime_error("c++ exception RcppQRHdf5 (Group IException)");
    } catch( H5::DataSetIException& error ) {
        throw std::runtime_error("c++ exception RcppQRHdf5 (DataSet IException)");
    } catch(std::exception& ex) {
        throw std::runtime_error(std::string("c++ exception RcppQRHdf5: ") + ex.what());
    }
    /**
    
    try {

        int irank; //,
            // iblockfactor = 1;
        std::vector<hsize_t> offset = {0,0},
            count = {dsA->nrows(), dsA->ncols()},
            stride = {1,1},
            block = {1,1};
        Eigen::MatrixXd Q;
        
        
        std::vector<double> vdA( count[0] * count[1] ); 
        dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
        Eigen::MatrixXd A = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdA.data(), count[0], count[1] );
        A.transposeInPlace();
        
        
        Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(A);
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
        
        qr.compute(A);
        irank = lu_decomp.rank();
        
        
        //..//if (irank == count[0] + 1 || irank == count[1] + 1 )
        if (static_cast<unsigned long long>(irank) == count[0] + 1ULL ||
            static_cast<unsigned long long>(irank) == count[1] + 1ULL)
        {
            Eigen::MatrixXd R = qr.matrixQR().template triangularView<Eigen::Upper>();
            dsR->inheritCompressionLevel(dsA->getCompressionLevel());
            dsR->createDataset( R.rows(), R.cols(), "real" );
            dsR->writeDataset(Rcpp::wrap(R));
        } else {
            Eigen::MatrixXd R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>();
            dsR->inheritCompressionLevel(dsA->getCompressionLevel());
            dsR->createDataset( R.rows(), R.cols(), "real" );
            dsR->writeDataset(Rcpp::wrap(R));
        }
        
        
        if (bthin == false)
        {
            Eigen::MatrixXd Q = qr.householderQ();
            dsQ->inheritCompressionLevel(dsA->getCompressionLevel());
            dsQ->createDataset( Q.rows(), Q.cols(), "real" );
            dsQ->writeDataset( Rcpp::wrap(Q) );
        } else {
            if (bthin) {
                Eigen::MatrixXd R = qr.matrixQR()
                                      .topLeftCorner(irank, A.cols())
                                      .triangularView<Eigen::Upper>();
                dsR->inheritCompressionLevel(dsA->getCompressionLevel());
                dsR->createDataset( R.rows(), R.cols(), "real" );
                dsR->writeDataset(Rcpp::wrap(R));
            } else {
                // Full QR: R is m×n upper trapezoidal
                Eigen::MatrixXd R = Eigen::MatrixXd::Zero(A.rows(), A.cols());
                R.topRows(irank) = qr.matrixQR()
                 .topRows(irank)
                 .triangularView<Eigen::Upper>();
                dsR->inheritCompressionLevel(dsA->getCompressionLevel());
                dsR->createDataset( R.rows(), R.cols(), "real" );
                dsR->writeDataset(Rcpp::wrap(R));
            }
        }
        
        //. 2026/03/02 .//  else {
            
            //. 2026/03/02 .//      //. 2025/01/15 error with thin calculus.// Eigen::MatrixXd Qthin = qr.householderQ() * Eigen::MatrixXd::Identity(count[0], count[1]);
            //. 2026/03/02 .//  //. 2025/01/15 .// dsQ->createDataset( Qthin.rows(), Qthin.cols(), "real" );
            //. 2026/03/02 .//  //. 2025/01/15 .// dsQ->writeDataset( Rcpp::wrap(Qthin));
            
            //. 2026/03/02 .//  Eigen::MatrixXd Qthin = qr.householderQ() * Eigen::MatrixXd::Identity(count[1], count[0]);
            //. 2026/03/02 .//  dsQ->createDataset( Qthin.rows(), Qthin.cols(), "real" );
            //. 2026/03/02 .//  dsQ->writeDataset( Rcpp::wrap(Qthin));
            
        //. 2026/03/02 .//  }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        throw std::runtime_error("c++ exception RcppQRHdf5 (File IException)");
    } catch( H5::GroupIException & error ) { // catch failure caused by the DataSet operations
        throw std::runtime_error("c++ exception RcppQRHdf5 (Group IException)");
    } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
        throw std::runtime_error("c++ exception RcppQRHdf5 (DataSet IException)");
    } catch(std::exception& ex) {
        throw std::runtime_error(std::string("c++ exception RcppQRHdf5: ") + ex.what());
    }
     **/
    
    return void();
    
}


// ── RcppTSQRHdf5 (Tall-Skinny QR — parallel row-block algorithm) ─────────────
//
// Algorithm outline (one-level TSQR):
//   1. Partition A (m×n) into p row-blocks A_1 … A_p  (each ≥ n rows).
//   2. Compute local QR in parallel:  A_b = Q_b R_b  (OpenMP for loop).
//   3. Stack local R factors:  R_stack = [R_1; R_2; …; R_p]  ((p·n)×n).
//   4. Compute QR of R_stack:  R_stack = Q_top · R_final.
//   5. Assemble global thin Q:  Q[block_b] = Q_b · Q_top[b·n:(b+1)·n, :].
//   6. Optionally expand thin Q (m×n) to full Q (m×m) via basis completion.
//
// References:
//   Demmel et al. (2012) "Communication-optimal parallel and sequential QR
//   and LU factorizations", SIAM J. Sci. Comput., 34(1):A206–A239.

inline void RcppTSQRHdf5( BigDataStatMeth::hdf5Dataset* dsA,
                           BigDataStatMeth::hdf5Dataset* dsQ,
                           BigDataStatMeth::hdf5Dataset* dsR,
                           bool bthin,
                           Rcpp::Nullable<int> block_size = R_NilValue,
                           Rcpp::Nullable<int> threads    = R_NilValue )
{
    try {

        // ── HDF5 dimension convention (package-wide) ──────────────────────────
        // nrows()   = HDF5 first dim  = R ncols  (n)
        // ncols()   = HDF5 second dim = R nrows  (m)
        const int n = static_cast<int>(dsA->nrows());   // R columns
        const int m = static_cast<int>(dsA->ncols());   // R rows

        if (n <= 0 || m <= 0)
            throw std::runtime_error("RcppTSQRHdf5: empty matrix");
        if (m < n)
            throw std::runtime_error(
                "RcppTSQRHdf5: matrix must have m >= n (use method='lapack' for wide/square matrices)");
        
        if (m < n)
            throw std::runtime_error("RcppTSQRHdf5: matrix must have m >= n");
        
        // Warn if not genuinely tall-skinny — LAPACK will be more efficient
        if (m < 2 * n)
            Rf_warning("method='tsqr' is inefficient for nearly-square matrices "
                           "(m < 2*n); consider method='lapack' or method='auto'");
        

        // ── Thread count (CRAN-compliant via get_number_threads) ──────────────
        const int nthreads = static_cast<int>( get_number_threads(threads, R_NilValue));
        
        // Guard: nthreads must be >= 1
        // get_number_threads returns unsigned; if cast gives <= 0, clamp to 1
        const int safe_nthreads = (nthreads > 0) ? nthreads : 1;
        
        // ── Block size: each block needs at least n rows for a rank-n local QR ─
        int bsize;
        if (block_size.isNotNull()) {
            bsize = std::max(n, Rcpp::as<int>(block_size));
        } else {
            // Auto: spread m rows evenly across threads, min n rows per block
            //.. 20260428 ..// bsize = std::max(n, static_cast<int>((m + nthreads - 1) / nthreads));
            bsize = std::max(n, static_cast<int>((m + safe_nthreads - 1) / safe_nthreads));
        }

        const int nblocks = (m + bsize - 1) / bsize;

        // ── Storage for local Q and R factors ─────────────────────────────────
        std::vector<Eigen::MatrixXd> local_R(nblocks);
        std::vector<Eigen::MatrixXd> local_Q(nblocks);

        const std::vector<hsize_t> stride = {1, 1}, blk = {1, 1};

        // ── Step 1: parallel local QR per row block ───────────────────────────
        //.. 20260428 ..// #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
        #pragma omp parallel for num_threads(safe_nthreads) schedule(dynamic)
        for (int b = 0; b < nblocks; ++b) {
            const int row_start = b * bsize;
            const int brows     = std::min(bsize, m - row_start);

            // Read block A[row_start : row_start+brows, :]  from HDF5.
            // In HDF5 storage (transposed): columns are R rows, rows are R cols.
            //   offset = {0, row_start}   (HDF5 dim0=R-cols, HDF5 dim1=R-rows)
            //   count  = {n, brows}
            std::vector<double> vd(static_cast<std::size_t>(n) * brows);
            #pragma omp critical(hdf5_tsqr_read)
            {
                dsA->readDatasetBlock(
                    {0, static_cast<hsize_t>(row_start)},
                    {static_cast<hsize_t>(n), static_cast<hsize_t>(brows)},
                    stride, blk, vd.data());
            }

            // Map to Eigen row-major (n × brows) then transpose → brows × n
            Eigen::MatrixXd Ablock =
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic,
                                         Eigen::Dynamic, Eigen::RowMajor>>(
                    vd.data(), n, brows);
            Ablock.transposeInPlace();   // now brows × n  (R dimensions)

            // Local thin QR:  Ablock (brows×n) = Q_b (brows×n) · R_b (n×n)
            Eigen::HouseholderQR<Eigen::MatrixXd> qr_local(Ablock);
            local_R[b] = qr_local.matrixQR().topRows(n)
                                  .template triangularView<Eigen::Upper>();
            // Thin Q: brows × n
            local_Q[b] = qr_local.householderQ() *
                         Eigen::MatrixXd::Identity(brows, n);
        }

        // ── Step 2: stack local R factors and compute top-level QR ───────────
        // R_stack is (nblocks·n) × n
        Eigen::MatrixXd R_stack(nblocks * n, n);
        for (int b = 0; b < nblocks; ++b)
            R_stack.block(b * n, 0, n, n) = local_R[b];

        Eigen::HouseholderQR<Eigen::MatrixXd> qr_top(R_stack);

        // Final R (n×n), upper triangular
        const Eigen::MatrixXd R_final =
            qr_top.matrixQR().topRows(n)
                              .template triangularView<Eigen::Upper>();

        // Q_top: (nblocks·n) × n  (thin factor of R_stack)
        const Eigen::MatrixXd Q_top =
            qr_top.householderQ() *
            Eigen::MatrixXd::Identity(nblocks * n, n);

        // ── Step 3: assemble global thin Q (m×n) ─────────────────────────────
        // Q[block_b rows] = local_Q[b]  ·  Q_top[b·n : (b+1)·n, :]
        Eigen::MatrixXd Q_thin(m, n);
        for (int b = 0; b < nblocks; ++b) {
            const int row_start = b * bsize;
            const int brows     = std::min(bsize, m - row_start);
            Q_thin.block(row_start, 0, brows, n) =
                local_Q[b] * Q_top.block(b * n, 0, n, n);
        }
        
        // ── Step 4 + 5: write Q and R ─────────────────────────────────────────────
        // R is written ONCE with the correct final dimensions — never twice.
        if (bthin) {
            // Thin R: n × n
            dsR->inheritCompressionLevel(dsA->getCompressionLevel());
            dsR->createDataset(n, n, "real");
            dsR->writeDataset(Rcpp::wrap(R_final));
            
            // Thin Q: m × n
            dsQ->inheritCompressionLevel(dsA->getCompressionLevel());
            dsQ->createDataset(m, n, "real");
            dsQ->writeDataset(Rcpp::wrap(Q_thin));
            
        } else {
            // Full Q: m × m via basis completion
            Eigen::HouseholderQR<Eigen::MatrixXd> qr_full(Q_thin);
            Eigen::MatrixXd Q_full = qr_full.householderQ()
                * Eigen::MatrixXd::Identity(m, m);
            const Eigen::VectorXd diag_signs =
            qr_full.matrixQR().diagonal().head(n).array().sign().matrix();
            for (int i = 0; i < n; ++i)
                Q_full.col(i) *= diag_signs(i);
            
            dsQ->inheritCompressionLevel(dsA->getCompressionLevel());
            dsQ->createDataset(m, m, "real");
            dsQ->writeDataset(Rcpp::wrap(Q_full));
            
            // Full R: m × n (padded with zeros below row n) — created ONCE
            Eigen::MatrixXd R_padded = Eigen::MatrixXd::Zero(m, n);
            R_padded.topRows(n) = R_final;
            dsR->inheritCompressionLevel(dsA->getCompressionLevel());
            dsR->createDataset(m, n, "real");
            dsR->writeDataset(Rcpp::wrap(R_padded));
        }

        // // ── Step 4: write R ───────────────────────────────────────────────────
        // dsR->inheritCompressionLevel(dsA->getCompressionLevel());
        // dsR->createDataset(n, n, "real");
        // dsR->writeDataset(Rcpp::wrap(R_final));
        // 
        // // ── Step 5: write Q ───────────────────────────────────────────────────
        // if (bthin) {
        //     // Thin Q: m × n
        //     dsQ->inheritCompressionLevel(dsA->getCompressionLevel());
        //     dsQ->createDataset(m, n, "real");
        //     dsQ->writeDataset(Rcpp::wrap(Q_thin));
        // } else {
        //     // Full Q: m × m
        //     // Basis completion: compute QR of Q_thin directly.
        //     // Since Q_thin already has orthonormal columns, its Householder QR
        //     // has R = I_n, so householderQ() returns a full m×m orthogonal matrix
        //     // whose first n columns ARE Q_thin.  This guarantees:
        //     //   Q_full[:, 1:n]  = Q_thin
        //     //   Q_full %*% R_padded = Q_thin %*% R_final = A  (exact reconstruction)
        //     //
        //     // Previous approach (QR of [Q_thin | random]) was INCORRECT because
        //     // Householder reorders column pivots and Q_thin is not preserved
        //     // as the first n columns of the result. (20260304)
        //     //
        //     // Cost: O(m² n) — expensive for large m; use thin=TRUE when possible.
        //     
        //     //.. 20260304 ..// Eigen::HouseholderQR<Eigen::MatrixXd> qr_full(Q_thin);
        //     //.. 20260304 ..// const Eigen::MatrixXd Q_full = qr_full.householderQ();
        // 
        //     // Gram-Schmidt: Q_full[:,1:n] = Q_thin por construcción → garantiza
        //     // Q_full * R_padded = Q_thin * R_final = A  (20260304)
        //     
        //     Eigen::HouseholderQR<Eigen::MatrixXd> qr_full(Q_thin);
        //     Eigen::MatrixXd Q_full = qr_full.householderQ()
        //         * Eigen::MatrixXd::Identity(m, m);
        //     // Correct sign flips: Q_full[:,i] *= sign(R_ii) for i=0..n-1
        //     // This makes Q_full[:,1:n] == Q_thin exactly (R diagonal is ±1 since Q_thin is orthonormal)
        //     const Eigen::VectorXd diag_signs =
        //     qr_full.matrixQR().diagonal().head(n).array().sign().matrix();
        //     for (int i = 0; i < n; ++i)
        //         Q_full.col(i) *= diag_signs(i);
        //     
        //     dsQ->inheritCompressionLevel(dsA->getCompressionLevel());
        //     dsQ->createDataset(m, m, "real");
        //     dsQ->writeDataset(Rcpp::wrap(Q_full));
        // 
        //     // Pad R to m × n (zeros below row n) to match full-QR convention
        //     Eigen::MatrixXd R_padded = Eigen::MatrixXd::Zero(m, n);
        //     R_padded.topRows(n) = R_final;
        //     // Re-create R dataset with updated dimensions
        //     dsR->inheritCompressionLevel(dsA->getCompressionLevel());
        //     dsR->createDataset(m, n, "real");
        //     dsR->writeDataset(Rcpp::wrap(R_padded));
        // }

    } catch( H5::FileIException& error ) {
        throw std::runtime_error("c++ exception RcppTSQRHdf5 (File IException)");
    } catch( H5::GroupIException& error ) {
        throw std::runtime_error("c++ exception RcppTSQRHdf5 (Group IException)");
    } catch( H5::DataSetIException& error ) {
        throw std::runtime_error("c++ exception RcppTSQRHdf5 (DataSet IException)");
    } catch(std::exception& ex) {
        throw std::runtime_error(std::string("c++ exception RcppTSQRHdf5: ") + ex.what());
    }

    return void();
}


// ── RcppQRHdf5_dispatch ───────────────────────────────────────────────────────

inline void RcppQRHdf5_dispatch( BigDataStatMeth::hdf5Dataset* dsA,
                                  BigDataStatMeth::hdf5Dataset* dsQ,
                                  BigDataStatMeth::hdf5Dataset* dsR,
                                  bool bthin,
                                  Rcpp::Nullable<int> block_size,
                                  Rcpp::Nullable<int> threads,
                                  const std::string& method )
{
    // R dimensions (corrected for HDF5 transposition)
    const int m = static_cast<int>(dsA->ncols());   // R rows
    const int n = static_cast<int>(dsA->nrows());   // R cols

    // Resolve effective method
    std::string eff_method = method;
    if (eff_method == "auto") {
        // TSQR is beneficial for tall-skinny matrices and when matrix is large
        // enough to amortise the parallel overhead
        eff_method = (m > 5 * n && m > 1000) ? "tsqr" : "lapack";
    }

    if (eff_method == "tsqr") {
        if (m < n)
            throw std::runtime_error(
                "method='tsqr' requires m >= n; use method='lapack' for wide matrices");
        RcppTSQRHdf5(dsA, dsQ, dsR, bthin, block_size, threads);
    } else if (eff_method == "lapack") {
        RcppQRHdf5(dsA, dsQ, dsR, bthin, block_size, threads);
    } else {
        throw std::runtime_error(
            std::string("RcppQRHdf5_dispatch: unknown method '") + eff_method +
            "'. Use 'auto', 'lapack', or 'tsqr'.");
    }
}


}   // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_HDF5_MATRIXQR_HPP
