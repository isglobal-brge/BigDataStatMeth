/**
 * @file hdf5Algebra.hpp
 * @brief Header file for HDF5-based matrix algebra operations
 * @details This header file provides a comprehensive set of matrix algebra operations
 * for large-scale data stored in HDF5 format. The implementation includes:
 * 
 * Key features:
 * - Matrix decompositions (SVD, PCA, QR)
 * - Matrix operations (multiplication, addition, subtraction)
 * - Statistical computations
 * - Linear algebra solvers
 * - Memory-efficient implementations
 * 
 * Components:
 * - matrixNormalization: Matrix normalization operations
 * - matrixSdMean: Standard deviation and mean calculations
 * - multiplication: Dense matrix multiplication
 * - matrixPCA: Principal Component Analysis
 * - matrixEquationSolver: Linear equation solvers
 * - matrixSvd: Singular Value Decomposition
 * - matrixSvdBlock: Block-based SVD for large matrices
 * - multiplicationSparse: Sparse matrix multiplication
 * - matrixQR: QR decomposition
 * - tcrossprod: Transposed cross-product
 * - matrixInvCholesky: Cholesky-based matrix inversion
 * - matrixSum: Matrix summation operations
 * - matrixDiagonal: Diagonal matrix operations
 * - vectormatrix: Vector-matrix operations
 * - matrixSubstract: Matrix subtraction
 * - matrixPseudoinverse: Moore-Penrose pseudoinverse
 * - matrixTriangular: Triangular matrix operations
 * - crossprod: Cross-product operations
 * 
 * The module is designed for efficient processing of large matrices
 * that cannot fit in memory, utilizing HDF5's chunked storage and
 * block-based algorithms.
 */

#ifndef BIGDATASTATMETH_HDF5ALGEBRA_HPP
#define BIGDATASTATMETH_HDF5ALGEBRA_HPP

    #include "hdf5Algebra/matrixNormalization.hpp"
    #include "hdf5Algebra/matrixSdMean.hpp"
    #include "hdf5Algebra/multiplication.hpp"
    #include "hdf5Algebra/matrixPCA.hpp"
    #include "hdf5Algebra/matrixEquationSolver.hpp"
    #include "hdf5Algebra/matrixSvd.hpp"
    #include "hdf5Algebra/matrixNormalization.hpp"
    #include "hdf5Algebra/matrixSvdBlock.hpp"
    #include "hdf5Algebra/multiplicationSparse.hpp"
    #include "hdf5Algebra/matrixQR.hpp"
    #include "hdf5Algebra/tcrossprod.hpp"
    #include "hdf5Algebra/matrixInvCholesky.hpp"
    #include "hdf5Algebra/matrixSum.hpp"
    #include "hdf5Algebra/matrixDiagonal.hpp"
    #include "hdf5Algebra/vectormatrix.hpp"
    #include "hdf5Algebra/matrixSubstract.hpp"
    #include "hdf5Algebra/matrixPseudoinverse.hpp"
    #include "hdf5Algebra/matrixTriangular.hpp"
    #include "hdf5Algebra/crossprod.hpp"


#endif // BIGDATASTATMETH_HD5ALGEBRA_HPP

