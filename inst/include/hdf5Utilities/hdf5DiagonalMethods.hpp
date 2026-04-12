/**
 * @file diagonalOperations.hpp
 * @brief Namespace for diagonal operations between HDF5 datasets
 * @note 2026-03-07 Output datasets now inherit compression level from input datasets
 *         via setCompressionLevel() called before every createDataset() invocation.
 * @details This file provides functionality for performing operations between
 * diagonal elements of matrices or vectors stored in HDF5 format.
 * 
 * Key features:
 * - Operations between matrix diagonals or vectors
 * - Automatic type detection (matrix vs vector)
 * - Flexible target specification (in-place or new dataset)
 * - Memory-efficient processing
 * - Parallel processing support
 */

#ifndef BIGDATASTATMETH_DIAGONAL_OPERATIONS_HPP
#define BIGDATASTATMETH_DIAGONAL_OPERATIONS_HPP

#include "hdf5CheckClose.hpp"
namespace BigDataStatMeth {
    namespace DiagonalOps {
    
        /**
         * @brief Validate vector dataset and return its size
         * @details Checks dataset dimensions and returns the number of elements in the vector.
         * Used for dimension validation before performing vector operations.
         * 
         * Size calculation:
         * - Row vector (1×N): returns N
         * - Column vector (N×1): returns N
         * - Scalar (1×1): returns 1
         * - Matrix (N×M): returns 0 (invalid)
         * 
         * @param ds Dataset to validate
         * @return Vector size if valid vector, 0 if not a vector
         * 
         * @note Return value 0 indicates invalid vector (matrix with N×M where N,M > 1)
         * @note Used for dimension compatibility checking before operations
         */
        inline hsize_t validateVectorDataset(BigDataStatMeth::hdf5Dataset* ds)
        {
            hsize_t rows = ds->nrows();
            hsize_t cols = ds->ncols();
            
            if (rows == 1 && cols > 1) {
                return cols;  // Row vector 1×N
            } else if (cols == 1 && rows > 1) {
                return rows;  // Column vector N×1
            } else if (rows == 1 && cols == 1) {
                return 1;     // Scalar 1×1
            }
            return 0;  // Not a vector (matrix N×M)
        }
        
        
        /**
         * @brief Check if dataset represents a diagonal vector
         * @details Validates dataset dimensions to determine if it's a vector (1×N or N×1).
         * A dataset is considered a diagonal vector if it has vector dimensions.
         * Future enhancements could include .diag suffix validation.
         * 
         * Validation criteria:
         * - Row vector: 1×N where N > 1
         * - Column vector: N×1 where N > 1  
         * - Scalar: 1×1 (edge case)
         * - Invalid: N×M where both N,M > 1
         * 
         * @param ds Dataset to check
         * @return true if dataset is a valid vector, false if matrix
         * 
         * @note Future versions may include .diag suffix validation
         * @note Scalar datasets (1×1) are considered vectors
         */
        inline bool isDiagonalVector(BigDataStatMeth::hdf5Dataset* ds)
        {
            hsize_t rows = ds->nrows();
            hsize_t cols = ds->ncols();
            
            return (rows == 1 && cols > 1) ||     // Row vector
                (cols == 1 && rows > 1) ||     // Column vector
                (rows == 1 && cols == 1);      // Scalar
        }
        
        
        /**
         * @brief Clean up temporary datasets created during operations
         * @details Safely deletes temporary datasets and releases memory.
         * Used internally to ensure proper resource cleanup even when exceptions occur.
         * 
         * Safety features:
         * - Null-pointer safe (handles nullptr gracefully)
         * - Exception safe (won't throw during cleanup)
         * - Memory leak prevention
         * - Used in exception handling paths
         * 
         * @param tempA Temporary dataset A (can be nullptr)
         * @param tempB Temporary dataset B (can be nullptr)
         * 
         * @note Silent cleanup - doesn't throw exceptions to avoid masking original errors
         * @note Called automatically by performMatrixDiagonalOperation()
         */
        inline void cleanup_temp_datasets(BigDataStatMeth::hdf5Dataset* tempA, 
                                          BigDataStatMeth::hdf5Dataset* tempB)
        {
            try {
                if (tempA != nullptr) {
                    delete tempA;
                    tempA = nullptr;
                }
                if (tempB != nullptr) {
                    delete tempB;
                    tempB = nullptr;
                }
            } catch(...) {
                // Silent cleanup
            }
        }
        
        
        
        /**
         * @brief Extract diagonal from matrix and save as vector dataset
         * @details Extracts diagonal elements from a matrix and creates a new vector dataset.
         * Uses existing getDiagonalfromMatrix() function for optimized diagonal reading.
         * The result vector is stored in 1×N format for compatibility with vector operations.
         * 
         * Extraction process:
         * - Validates square matrix requirement
         * - Uses getDiagonalfromMatrix() for efficient extraction
         * - Creates 1×N vector dataset for storage efficiency
         * - Preserves data precision and handles large matrices
         * 
         * @param dsMatrix Input matrix dataset (must be square N×N)
         * @param dsVector Output vector dataset (will be created as 1×N)
         * 
         * @throws std::exception if matrix is not square or extraction fails
         * @note Matrix must be square (N×N) for diagonal extraction
         * @note Output vector is stored as 1×N for compatibility with vector operations
         * @see getDiagonalfromMatrix() for the underlying extraction algorithm
         */
        inline void extractDiagonalToVector(BigDataStatMeth::hdf5Dataset* dsMatrix, 
                                            BigDataStatMeth::hdf5Dataset* dsVector)
        {
            try {
                if (dsMatrix->nrows() != dsMatrix->ncols()) {
                    throw std::runtime_error("extractDiagonalToVector: Matrix must be square");
                }
                
                hsize_t matrix_size = dsMatrix->nrows();
                
                // Extract diagonal using existing function
                Rcpp::NumericVector diagonal = getDiagonalfromMatrix(dsMatrix);
                
                // Create vector dataset (1×N format)
                dsVector->inheritCompressionLevel(dsMatrix->getCompressionLevel());
                //.. 17/03/2026 ..// dsVector->createDataset(matrix_size, 1, "numeric");
                dsVector->createDataset(1, matrix_size, "numeric");
                
                // Write diagonal data
                std::vector<double> diag_vector = Rcpp::as<std::vector<double>>(diagonal);
                std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
                //.. 17/03/2026 ..// dsVector->writeDatasetBlock(diag_vector, {0, 0}, {1, matrix_size}, stride, block);
                dsVector->writeDatasetBlock(diag_vector, {0, 0}, {matrix_size, 1}, stride, block);
                
            } catch(std::exception& ex) {
                throw std::runtime_error(std::string("Error in extractDiagonalToVector: ") + ex.what());
            }
        }
        
        
        
        /**
         * @brief Write diagonal vector to matrix diagonal
         * @details Takes a vector dataset and writes its values to the diagonal of a matrix.
         * Uses existing setDiagonalMatrix() function for optimized diagonal writing.
         * Matrix must already exist and be square. Vector must match matrix diagonal size.
         * 
         * Writing process:
         * - Validates vector and matrix dimensions compatibility
         * - Reads vector data efficiently
         * - Uses setDiagonalMatrix() for optimized diagonal writing
         * - Preserves all non-diagonal matrix elements
         * 
         * @param dsVector Input vector dataset containing diagonal values (1×N or N×1)
         * @param dsMatrix Target matrix dataset (must exist and be square N×N)
         * 
         * @throws std::exception if dimensions don't match or writing fails
         * @note Matrix must already exist with proper dimensions
         * @note Vector size must match matrix diagonal size (N elements for N×N matrix)
         * @see setDiagonalMatrix() for the underlying writing algorithm
         */
        inline void writeDiagonalFromVector(BigDataStatMeth::hdf5Dataset* dsVector,
                                            BigDataStatMeth::hdf5Dataset* dsMatrix)
        {
            try {
                if (dsMatrix->nrows() != dsMatrix->ncols()) {
                    throw std::runtime_error("writeDiagonalFromVector: Matrix must be square");
                }
                
                hsize_t vector_size = validateVectorDataset(dsVector);
                hsize_t matrix_size = dsMatrix->nrows();
                
                if (vector_size == 0) {
                    throw std::runtime_error("writeDiagonalFromVector: Input is not a valid vector");
                }
                
                if (vector_size != matrix_size) {
                    throw std::runtime_error("writeDiagonalFromVector: Vector size must match matrix diagonal size");
                }
                
                // Read vector data
                std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
                std::vector<double> vector_data(vector_size);
                
                if (dsVector->nrows() == 1) {
                    dsVector->readDatasetBlock({0, 0}, {1, vector_size}, stride, block, vector_data.data());
                } else {
                    dsVector->readDatasetBlock({0, 0}, {vector_size, 1}, stride, block, vector_data.data());
                }
                
                // Write using existing function
                Rcpp::NumericVector diagonal_values = Rcpp::wrap(vector_data);
                setDiagonalMatrix(dsMatrix, diagonal_values);
                
            } catch(std::exception& ex) {
                throw std::runtime_error(std::string("Error in writeDiagonalFromVector: ") + ex.what());
            }
        }
        
        /**
         * @brief Perform diagonal operations on matrices using extract-operate-write strategy
         * @details Implements the extract-operate-write pattern for matrix diagonal operations.
         * It processes only diagonal elements instead of full matrices.
         * 
         * Strategy implementation:
         * 1. Extract diagonals from matrices (if needed) using extractDiagonalToVector()
         * 2. Perform vector operation using optimized functions from vectorOperations.hpp
         * 3. Result is stored as vector (most efficient) or written back to matrix diagonal
         * 
         * Operation codes:
         * - 0: Addition (A_diag + B_diag)
         * - 1: Subtraction (A_diag - B_diag)  
         * - 2: Multiplication (A_diag * B_diag)
         * - 3: Division (A_diag / B_diag)
         * 
         * @param dsA First input dataset
         * @param dsB Second input dataset  
         * @param dsResult Result dataset (only used if target="new")
         * @param operation Operation type: 0=add, 1=subtract, 2=multiply, 3=divide
         * @param target Where to write result: "A", "B", or "new"
         * @param bparal Whether to use parallel processing
         * @param threads Number of threads
         * 
         * @throws std::exception if matrices are not square or dimensions don't match
         * @note Creates temporary vector datasets for matrix inputs (cleaned up automatically)
         * @note Validates all dimensions before performing operations
         * @note Uses existing getDiagonalfromMatrix() and setDiagonalMatrix() functions
         */
        inline void performMatrixDiagonalOperation(BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB,
                                                   BigDataStatMeth::hdf5Dataset* dsResult, int operation, std::string target,
                                                   bool bparal, Rcpp::Nullable<int> threads)
        {
            
            BigDataStatMeth::hdf5Dataset* operationTarget = nullptr;
            
            try {

                // BigDataStatMeth::hdf5Dataset* tempA = nullptr;
                // BigDataStatMeth::hdf5Dataset* tempB = nullptr;
                // BigDataStatMeth::hdf5Dataset* tempResult = nullptr;
                std::unique_ptr<BigDataStatMeth::hdf5Dataset> tempA(nullptr);
                std::unique_ptr<BigDataStatMeth::hdf5Dataset> tempB(nullptr);
                std::unique_ptr<BigDataStatMeth::hdf5Dataset> tempResult(nullptr);


                bool isVectorA = isDiagonalVector(dsA);
                bool isVectorB = isDiagonalVector(dsB);
                
                BigDataStatMeth::hdf5Dataset* finalA = dsA;
                BigDataStatMeth::hdf5Dataset* finalB = dsB;
                
                // Extract diagonal from A if it's a matrix
                if (!isVectorA) {
                    if (dsA->nrows() != dsA->ncols()) {
                        // checkClose_file(tempA, tempB, tempResult);
                        throw std::runtime_error("Matrix A must be square for diagonal operations");
                    }
                    std::string tempNameA = dsA->getDatasetName() + "_temp_diag_A";
                    // tempA = new BigDataStatMeth::hdf5Dataset(dsA->getFileptr(), dsA->getGroup(), tempNameA, true);
                    tempA.reset(new BigDataStatMeth::hdf5Dataset(dsA->getFileptr(), dsA->getGroup(), tempNameA, true));
                    extractDiagonalToVector(dsA, tempA.get());
                    finalA = tempA.get();
                }
                
                // Extract diagonal from B if it's a matrix
                if (!isVectorB) {
                    if (dsB->nrows() != dsB->ncols()) {
                        // checkClose_file(tempA, tempB, tempResult);
                        throw std::runtime_error("Matrix B must be square for diagonal operations");
                        // cleanup_temp_datasets(tempA, tempB);
                        return;
                    }
                    std::string tempNameB = dsB->getDatasetName() + "_temp_diag_B";
                    // tempB = new BigDataStatMeth::hdf5Dataset(dsB->getFileptr(), dsB->getGroup(), tempNameB, true);
                    tempB.reset(new BigDataStatMeth::hdf5Dataset(dsB->getFileptr(), dsB->getGroup(), tempNameB, true));
                    extractDiagonalToVector(dsB, tempB.get());
                    finalB = tempB.get();
                }
                
                // Validate dimensions
                hsize_t sizeA = validateVectorDataset(finalA);
                hsize_t sizeB = validateVectorDataset(finalB);
                
                if (sizeA == 0 || sizeB == 0 || sizeA != sizeB) {
                    // checkClose_file(tempA, tempB, tempResult);
                    throw std::runtime_error("Invalid or incompatible diagonal dimensions: vs");
                    // cleanup_temp_datasets(tempA, tempB);
                    return;
                }
                
                // Determine target for operation result
                
                // std::unique_ptr<BigDataStatMeth::hdf5Dataset> operationTarget(nullptr);
                
                
                if (target == "new") {
                    operationTarget = dsResult;
                } else if (target == "A") {
                    if (isVectorA) {
                        operationTarget = dsA;
                    } else {
                        std::string tempNameResult = dsA->getDatasetName() + "_temp_result";
                        // tempResult = new BigDataStatMeth::hdf5Dataset(dsA->getFileptr(), dsA->getGroup(), tempNameResult, true);
                        tempResult.reset(new BigDataStatMeth::hdf5Dataset(dsA->getFileptr(), dsA->getGroup(), tempNameResult, true));
                        operationTarget = tempResult.get();
                    }
                } else if (target == "B") {
                    if (isVectorB) {
                        operationTarget = dsB;
                    } else {
                        std::string tempNameResult = dsB->getDatasetName() + "_temp_result";
                        // tempResult = new BigDataStatMeth::hdf5Dataset(dsB->getFileptr(), dsB->getGroup(), tempNameResult, true);
                        tempResult.reset(new BigDataStatMeth::hdf5Dataset(dsB->getFileptr(), dsB->getGroup(), tempNameResult, true));
                        operationTarget = tempResult.get();
                    }
                }
                
                // Perform vector operation
                switch (operation) {
                    case 0: Rcpp_vector_add_hdf5(finalA, finalB, operationTarget, bparal, threads); break;
                    case 1: Rcpp_vector_subtract_hdf5(finalA, finalB, operationTarget, bparal, threads); break;
                    case 2: Rcpp_vector_multiply_hdf5(finalA, finalB, operationTarget, bparal, threads); break;
                    case 3: Rcpp_vector_divide_hdf5(finalA, finalB, operationTarget, bparal, threads); break;
                    throw std::runtime_error("Unknown diagonal operation");
                }
                
                // Write result back to matrix diagonal if needed
                if (target == "A" && !isVectorA) {
                    writeDiagonalFromVector(tempResult.get(), dsA);
                } else if (target == "B" && !isVectorB) {
                    writeDiagonalFromVector(tempResult.get(), dsB);
                }
                
                //......// delete operationTarget; operationTarget = nullptr;
                // Cleanup
                // cleanup_temp_datasets(tempA, tempB);
                // if (tempResult) { delete tempResult; tempResult = nullptr; }
                
            } catch(std::exception& ex) {
                checkClose_file(operationTarget);
                // checkClose_file(tempA, tempB, tempResult);
                // cleanup_temp_datasets(tempA, tempB);
                // if (tempResult) { delete tempResult; tempResult = nullptr; }
                throw std::runtime_error(std::string("Error in performMatrixDiagonalOperation: ") + ex.what());
            }
        }
        
        
        
        /**
         * @brief Add diagonal elements from two matrices or vectors
         * @details Performs optimized diagonal addition C_diag = A_diag + B_diag.
         * Automatically detects whether inputs are matrices (extracts diagonals) or 
         * vectors (direct vector operation). Uses vectorOperations.hpp for maximum efficiency.
         * 
         * Operation modes:
         * - Matrix + Matrix: Extract diagonals → vector addition → write diagonal
         * - Matrix + Vector: Extract diagonal → vector addition → write diagonal  
         * - Vector + Vector: Direct vector addition (most efficient)
         * 
         * @param dsA First input dataset (matrix or vector)
         * @param dsB Second input dataset (matrix or vector)  
         * @param dsResult Result dataset (only used if target="new")
         * @param target Where to write result: "A", "B", or "new"
         * @param bparal Whether to use parallel processing
         * @param threads Number of threads for parallel processing
         * 
         * @throws std::exception if dimension mismatch or invalid inputs
         * @note Result dataset name automatically gets .diag suffix if inputs are vectors
         */
        inline void addDiagonals(BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB, 
                                 BigDataStatMeth::hdf5Dataset* dsResult, std::string target = "new",
                                 bool bparal = false, Rcpp::Nullable<int> threads = R_NilValue)
        {
            BigDataStatMeth::hdf5Dataset* targetDataset = nullptr;

            try {
                bool isVectorA = isDiagonalVector(dsA);
                bool isVectorB = isDiagonalVector(dsB);
                
                if (isVectorA && isVectorB && (target == "A" || target == "B")) {
                    targetDataset = (target == "A") ? dsA : dsB;
                    Rcpp_vector_add_hdf5(dsA, dsB, targetDataset, bparal, threads);
                    // delete targetDataset; targetDataset = nullptr;
                } else if (isVectorA && isVectorB && target == "new") {
                    Rcpp_vector_add_hdf5(dsA, dsB, dsResult, bparal, threads);
                } else {
                    performMatrixDiagonalOperation(dsA, dsB, dsResult, 0, target, bparal, threads);
                }
            } catch(std::exception& ex) {
                checkClose_file(targetDataset);
                throw std::runtime_error(std::string("Error in addDiagonals: ") + ex.what());
            }
        }
        
        
        /**
         * @brief Subtract diagonal elements from two matrices or vectors
         * @details Performs optimized diagonal subtraction C_diag = A_diag - B_diag.
         * Same optimization strategy as addDiagonals but for subtraction operation.
         * 
         * @param dsA First input dataset (minuend)
         * @param dsB Second input dataset (subtrahend)
         * @param dsResult Result dataset (will be created)
         * @param target Where to write result: "A", "B", or "new"
         * @param bparal Whether to use parallel processing
         * @param threads Number of threads for parallel processing
         * 
         * @throws std::exception if dimension mismatch or invalid inputs
         * @see addDiagonals() for detailed operation modes
         */
        inline void subtractDiagonals(BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB,
                                      BigDataStatMeth::hdf5Dataset* dsResult, std::string target = "new",
                                      bool bparal = false, Rcpp::Nullable<int> threads = R_NilValue)
        {
            BigDataStatMeth::hdf5Dataset* targetDataset = nullptr;

            try {
    
                bool isVectorA = isDiagonalVector(dsA);
                bool isVectorB = isDiagonalVector(dsB);
                
                if (isVectorA && isVectorB && (target == "A" || target == "B")) {
                    targetDataset = (target == "A") ? dsA : dsB;
                    Rcpp_vector_subtract_hdf5(dsA, dsB, targetDataset, bparal, threads);
                    // delete targetDataset; targetDataset = nullptr;
                } else if (isVectorA && isVectorB && target == "new") {
                    Rcpp_vector_subtract_hdf5(dsA, dsB, dsResult, bparal, threads);
                } else {
                    performMatrixDiagonalOperation(dsA, dsB, dsResult, 1, target, bparal, threads);
                }
            } catch(std::exception& ex) {
                checkClose_file(targetDataset);
                throw std::runtime_error(std::string("Error in subtractDiagonals: ") + ex.what());
            }
        }
        
        
        
        
        /**
         * @brief Multiply diagonal elements from two matrices or vectors
         * @details Performs optimized diagonal multiplication C_diag = A_diag * B_diag.
         * Same optimization strategy as addDiagonals but for element-wise multiplication.
         * 
         * @param dsA First input dataset
         * @param dsB Second input dataset
         * @param dsResult Result dataset (will be created)
         * @param target Where to write result: "A", "B", or "new"
         * @param bparal Whether to use parallel processing
         * @param threads Number of threads for parallel processing
         * 
         * @throws std::exception if dimension mismatch or invalid inputs
         * @note This is element-wise multiplication, not matrix multiplication
         */
        inline void multiplyDiagonals(BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB,
                                      BigDataStatMeth::hdf5Dataset* dsResult, std::string target = "new",
                                      bool bparal = false, Rcpp::Nullable<int> threads = R_NilValue)
        {
            BigDataStatMeth::hdf5Dataset* targetDataset = nullptr;

            try {
                
                bool isVectorA = isDiagonalVector(dsA);
                bool isVectorB = isDiagonalVector(dsB);
                
                if (isVectorA && isVectorB && (target == "A" || target == "B")) {
                    targetDataset = (target == "A") ? dsA : dsB;
                    Rcpp_vector_multiply_hdf5(dsA, dsB, targetDataset, bparal, threads);
                    // delete targetDataset; targetDataset = nullptr;
                } else if (isVectorA && isVectorB && target == "new") {
                    Rcpp_vector_multiply_hdf5(dsA, dsB, dsResult, bparal, threads);
                } else {
                    performMatrixDiagonalOperation(dsA, dsB, dsResult, 2, target, bparal, threads);
                }
            } catch(std::exception& ex) {
                checkClose_file(targetDataset);
                throw std::runtime_error(std::string("Error in multiplyDiagonals: ") + ex.what());
            }
        }
        
        
        
        
        /**
         * @brief Divide diagonal elements from two matrices or vectors
         * @details Performs optimized diagonal division C_diag = A_diag / B_diag.
         * Same optimization strategy as addDiagonals but for element-wise division.
         * Division by zero follows IEEE 754 standard (results in infinity).
         * 
         * @param dsA First input dataset (dividend)
         * @param dsB Second input dataset (divisor)
         * @param dsResult Result dataset (will be created)
         * @param target Where to write result: "A", "B", or "new"
         * @param bparal Whether to use parallel processing
         * @param threads Number of threads for parallel processing
         * 
         * @throws std::exception if dimension mismatch or invalid inputs
         * @note Division by zero results in infinity (IEEE 754 standard behavior)
         */
        inline void divideDiagonals(BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB,
                                    BigDataStatMeth::hdf5Dataset* dsResult, std::string target = "new",
                                    bool bparal = false, Rcpp::Nullable<int> threads = R_NilValue)
        {
            BigDataStatMeth::hdf5Dataset* targetDataset = nullptr;

            try {
                
                bool isVectorA = isDiagonalVector(dsA);
                bool isVectorB = isDiagonalVector(dsB);
                
                if (isVectorA && isVectorB && (target == "A" || target == "B")) {
                    BigDataStatMeth::hdf5Dataset* targetDataset = (target == "A") ? dsA : dsB;
                    Rcpp_vector_divide_hdf5(dsA, dsB, targetDataset, bparal, threads);
                    // delete targetDataset; targetDataset = nullptr;
                } else if (isVectorA && isVectorB && target == "new") {
                    Rcpp_vector_divide_hdf5(dsA, dsB, dsResult, bparal, threads);
                } else {
                    performMatrixDiagonalOperation(dsA, dsB, dsResult, 3, target, bparal, threads);
                }
            } catch(std::exception& ex) {
                checkClose_file(targetDataset);
                throw std::runtime_error(std::string("Error in divideDiagonals: ") + ex.what());
            }
        }
        
        
        /**
         * @brief Divide diagonal elements from two matrices or vectors
         * @details Performs optimized diagonal power C_diag = A_diag ^ B_diag.
         * Same optimization strategy as addDiagonals but for element-wise division.
         * 
         * @param dsA First input dataset (dividend)
         * @param dsB Second input dataset (divisor)
         * @param dsResult Result dataset (will be created)
         * @param target Where to write result: "A", "B", or "new"
         * @param bparal Whether to use parallel processing
         * @param threads Number of threads for parallel processing
         * 
         * @throws std::exception if dimension mismatch or invalid inputs
         */
        inline void powerDiagonals(BigDataStatMeth::hdf5Dataset* dsA, BigDataStatMeth::hdf5Dataset* dsB,
                                    BigDataStatMeth::hdf5Dataset* dsResult, std::string target = "new",
                                    bool bparal = false, Rcpp::Nullable<int> threads = R_NilValue)
        {
            BigDataStatMeth::hdf5Dataset* targetDataset = nullptr;

            try {
                
                bool isVectorA = isDiagonalVector(dsA);
                bool isVectorB = isDiagonalVector(dsB);
                
                if (isVectorA && isVectorB && (target == "A" || target == "B")) {
                    targetDataset = (target == "A") ? dsA : dsB;
                    Rcpp_vector_power_hdf5(dsA, dsB, targetDataset, bparal, threads);
                    // delete targetDataset; targetDataset = nullptr;
                } else if (isVectorA && isVectorB && target == "new") {
                    Rcpp_vector_power_hdf5(dsA, dsB, dsResult, bparal, threads);
                } else {
                    performMatrixDiagonalOperation(dsA, dsB, dsResult, 3, target, bparal, threads);
                }
            } catch(std::exception& ex) {
                checkClose_file(targetDataset);
                throw std::runtime_error(std::string("Error in powerDiagonals: ") + ex.what());
            }
        }
    
        /**
         * @brief Perform scalar operations on diagonal elements
         * @details Applies scalar operations (add, subtract, multiply, divide) to diagonal
         * elements. Automatically detects if input is matrix (extracts diagonal) or vector (direct).
         * 
         * @param dsInput Input dataset (matrix or vector)
         * @param dsResult Result dataset (only used if target="new")
         * @param scalar Scalar value for operation
         * @param operation Operation type: 0=add, 1=subtract, 2=multiply, 3=divide
         * @param target Where to write result: "input", "new"
         * @param bparal Whether to use parallel processing
         * @param threads Number of threads
         */
        inline void scalarOperation(BigDataStatMeth::hdf5Dataset* dsInput,
                                    BigDataStatMeth::hdf5Dataset* dsResult,
                                    double scalar,
                                    int operation,
                                    std::string target = "new",
                                    bool bparal = false,
                                    Rcpp::Nullable<int> threads = R_NilValue)
        {
            
            try {
                
                H5::Exception::dontPrint();
                
                // BigDataStatMeth::hdf5Dataset* tempInput = nullptr;
                // BigDataStatMeth::hdf5Dataset* tempResult = nullptr;
                std::unique_ptr<BigDataStatMeth::hdf5Dataset> tempInput(nullptr);
                std::unique_ptr<BigDataStatMeth::hdf5Dataset> tempResult(nullptr);

                BigDataStatMeth::hdf5Dataset* finalInput = dsInput;


                bool isVectorInput = isDiagonalVector(dsInput);

                // Extract diagonal from input if it's a matrix
                if (!isVectorInput) {

                    if (dsInput->nrows() != dsInput->ncols()) {
                        throw std::runtime_error("Input matrix must be square for diagonal operations");
                    }

                    std::string tempNameInput = dsInput->getDatasetName() + "_temp_scalar_input";
                    // tempInput = new BigDataStatMeth::hdf5Dataset(dsInput->getFileptr(), dsInput->getGroup(), tempNameInput, true);
                    tempInput.reset(new BigDataStatMeth::hdf5Dataset(dsInput->getFileptr(), dsInput->getGroup(), tempNameInput, true));
                    extractDiagonalToVector(dsInput, tempInput.get());
                    finalInput = tempInput.get();
                }
                
                // Validate input
                hsize_t sizeInput = validateVectorDataset(finalInput);
                if (sizeInput == 0) {
                    throw std::runtime_error("Invalid input dimensions for scalar operation");
                    // cleanup_temp_datasets(tempInput, nullptr);
                    return;
                }
                
                // Determine target for operation result
                BigDataStatMeth::hdf5Dataset* operationTarget = nullptr;

                if (target == "new") {
                    operationTarget = dsResult;
                    // operationTarget->createDataset(1, sizeInput, "real");  // Vector 1×N
                    operationTarget->inheritCompressionLevel(operationTarget->getCompressionLevel());
                    operationTarget->createDataset( sizeInput, 1, "real");  // Vector 1×N
                } else if (target == "input") {
                    if (isVectorInput) {
                        operationTarget = dsInput;  
                    } else {
                        // Para matriz, crear temp y después escribir diagonal
                        std::string tempNameResult = dsInput->getDatasetName() + "_temp_scalar_result";
                        // tempResult = new BigDataStatMeth::hdf5Dataset(dsInput->getFileptr(), dsInput->getGroup(), tempNameResult, true);
                        tempResult.reset(new BigDataStatMeth::hdf5Dataset(dsInput->getFileptr(), dsInput->getGroup(), tempNameResult, true));
                        // tempResult->createDataset(1, sizeInput, "real");
                        tempResult->inheritCompressionLevel(tempResult->getCompressionLevel());
                        tempResult->createDataset(sizeInput, 1, "real"); 
                        operationTarget = tempResult.get();
                    }
                }
                
                // Read input vector data
                std::vector<hsize_t> stride = {1, 1}, block = {1, 1};
                std::vector<double> input_data(sizeInput);
                
                if (finalInput->nrows() == 1) {
                    finalInput->readDatasetBlock({0, 0}, {1, sizeInput}, stride, block, input_data.data());
                } else {
                    finalInput->readDatasetBlock({0, 0}, {sizeInput, 1}, stride, block, input_data.data());
                }
                
                // Apply scalar operation
                if (bparal && sizeInput > 10000) {
                    #pragma omp parallel num_threads(get_threads(bparal, threads))
                    {
                    #pragma omp for schedule(static)
                        for (hsize_t i = 0; i < sizeInput; ++i) {
                            switch (operation) {
                            case 0: input_data[i] += scalar; break;  // add
                            case 1: input_data[i] -= scalar; break;  // subtract
                            case 2: input_data[i] *= scalar; break;  // multiply
                            case 3: input_data[i] /= scalar; break;  // divide
                            case 4: input_data[i] = std::pow(input_data[i], scalar); break;  // power
                            throw std::runtime_error("Unknown scalar operation");
                            }
                        }
                    }
                } else {
                    for (hsize_t i = 0; i < sizeInput; ++i) {
                        switch (operation) {
                        case 0: input_data[i] += scalar; break;
                        case 1: input_data[i] -= scalar; break;
                        case 2: input_data[i] *= scalar; break;
                        case 3: input_data[i] /= scalar; break;
                        case 4: input_data[i] = std::pow(input_data[i], scalar); break;
                        throw std::runtime_error("Unknown scalar operation");
                        }
                    }
                }
                
                // Write result
                operationTarget->writeDatasetBlock(input_data, {0, 0}, {1, sizeInput}, stride, block);
                
                // Write result back to matrix diagonal if needed
                if (target == "input" && !isVectorInput) {
                    writeDiagonalFromVector(tempResult.get(), dsInput);
                }
                
                // Cleanup
                // cleanup_temp_datasets(tempInput, nullptr);
                // if (tempResult) { delete tempResult; tempResult = nullptr; }
            } catch(std::exception& ex) {
                // cleanup_temp_datasets(tempInput, nullptr);
                // if (tempResult) { delete tempResult; tempResult = nullptr; }
                throw std::runtime_error(std::string("Error in scalarOperation: ") + ex.what());
            }
        }
    
    
    } // namespace DiagonalOps
} // namespace BigDataStatMeth

#endif // BIGDATASTATMETH_DIAGONAL_OPERATIONS_HPP
