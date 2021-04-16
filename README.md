# Big Data Statistical Methods (BigDataStatMeth)
 

This package implements basic Algebra methods using parallel algorithms to be used in big data problems such as omic data analyses. The functions will consider as input an object of class DelayedArray which is an array-like object (typically an on-disk object) specifically designed for Bioconductor. The methods will also accept other type of objects as those designed to deal with big matrices in the Matrix package. The implemented methods will include:

- Algebra

     - Matrix multiplication (by blocks, parallel, ...)
     - Cross-product and transpose cross-product
     - Matrix vector multiplication
     - Weighted cross-product and weighted transposed cross-product
     - Inverse Cholesky
     - Singular Value Decomposition
     - Jacobi decomposition


- Statistical methods
     - Data reduction (Principal Component Analysis, Canonical Correlation, Partial Least Squares, ...)
     - Linear regression models
     - Penalized methods for linear regression (Lasso, Elastic Net, ...)
     - Omic data integration (Penalized Generalized Canonical Correlation)
     - Supervised and non-supervised algorithms


## Implemented Functions
 
|                                                  |                        R Function                       | Delayed arrays | By blocks | Parallel | HDF5 |
|--------------------------------------------------|:-------------------------------------------------------:|----------------|:---------:|----------|------|
| **Basic functions with vectors and matrices**    |                                                         |                |           |          |      |
| Matrix Product                                   |             blockmult(…) / blockmult_hdf5(…)            |        ✔︎       |     ✔︎     |     ✔︎    |   ✔︎  |
| Sparse Matrix Product                            |      blockmult_sparse(…) / blockmult_sparse_hdf5(…)     |                |           |          |   ✔︎  |
| Matrix product with its transpose                | bdcrossprod(…) / Crossprod_hdf5(…) / tCrossprod_hdf5(…) |        ✔︎       |     ✔︎     |     ✔︎    |   ✔︎  |
| Matrix - Matrix weighted product (XWXt, XtWX)    |      tCrossprod_Weighted(…) / Crossprod_Weighted(…)     |        ✔︎       |     ✔︎     |     ✔︎    |   ✔︎  |
| Matrix - vector weighted product (XwXt, XtwX)    |           bdwproduct(…) / bdScalarwproduct(…)           |        ✔︎       |           |          |      |
| Matrix vector product                            |                   blockmult_vector(…)                   |        ✔︎       |     ✔︎     |     ✔︎    |      |
| Data Normalization (center, scale and both)      |          Normalize_Data(…) / Normalize_hdf5(…)          |        ✔︎       |     ✔︎     |          |   ✔︎  |
| **Other functions**                              |                                                         |                |           |          |      |
| Vector sum                                       |                   parallelVectorSum(…)                  |                |           |     ✔︎    |      |
| Pow(2) vector elements                           |                     parallelpow2(…)                     |                |           |     ✔︎    |      |
| **Lineal Algebra Functions**                     |                                                         |                |           |          |      |
| SVD matrix decomposition                         |                 bdSVD(…) / bdSVD_hdf5(…)                |        ✔︎       |     ✔︎     |     ✔︎    |   ✔︎  |
| QR matrix decomposition                          |                         bdQR(…)                         |        ✔︎       |           |          |      |
| Cholesky decomposition                           |                     bdInvCholesky(…)                    |        ✔︎       |           |          |      |
| Matrix Pseudoinverse                             |                      bdpseudoinv(…)                     |        ✔︎       |           |          |      |
| Solve matrix equation (A * X = B )               |                        bdSolve(…)                       |        ✔︎       |           |          |      |
| **Data Analysis**                                |                                                         |                |           |          |      |
| Principal Components Analysis (PCA)              |                      bdPCA_hdf5(…)                      |                |     ✔︎     |          |   ✔︎  |
| MLR-MR (Linear Regression Big Data)              |                       bdMLR_MR(…)                       |        ✔︎       |     ✔︎     |     ✔︎    |      |
| **HDF5 data files Utils**                        |                                                         |                |           |          |      |
| Remove rows or columns with hight missing values |                    bdRemovelowdata(…)                   |                |     ✔︎     |          |   ✔︎  |
| Impute missing data                              |                    bdImputeSNPHDF5(…)                   |                |     ✔︎     |          |   ✔︎  |
| Create hdf5 data file with one dataset inside    |                Create_HDF5_matrix_file(…)               |                |           |          |   ✔︎  |
| Add one dataset in hdf5 data file                |                  Create_HDF5_matrix(…)                  |                |           |          |   ✔︎  |
| Remove one dataset from hdf5 data file           |                  Remove_HDF5_element(…)                 |                |           |          |   ✔︎  |
| Import data from text file to HDF5               |                 Convert_text_to_HDF5(…)                 |                |           |          |   ✔︎  |