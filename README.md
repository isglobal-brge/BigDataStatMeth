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
     - Omics data integration (Penalized Generalized Canonical Correlation)
     - Supervised and non-supervised algorithms


## Install Package

To install BigDataStatMeth package we only have to run the following commands from the R command shell 

```{r, install, eval=FALSE}

# Install devtools and load library (if not previously installed)
install.packages("devtools") 
library(devtools)

# Install BigDataStatMeth and load package
install_github("isglobal-brge/BigDataStatMeth@HEAD")
library(BigDataStatMeth)

```

## Implemented Functions
 
|                                                  |                        R Function                       | Delayed arrays | By blocks | Parallel | HDF5 |
|--------------------------------------------------|:-------------------------------------------------------:|----------------|:---------:|----------|------|
| **Basic functions with vectors and matrices**    |                                                         |                |           |          |      |
| Matrix Product                                   |             bdblockmult(…) / bdblockmult_hdf5(…)            |        ✔︎       |     ✔︎     |     ✔︎    |   ✔︎  |
| Sparse Matrix Product                            |      bdblockmult_sparse(…) / bdblockmult_sparse_hdf5(…)     |                |           |          |   ✔︎  |
| Matrix product with its transpose                | bdCrossprod(…) / bdtCrossprod(…) / bdCrossprod_hdf5(…) / bdtCrossprod_hdf5(…) |        ✔︎       |     ✔︎     |     ✔︎    |   ✔︎  |
| Matrix - Matrix weighted product (XWXt, XtWX)    |      bdtCrossprod_Weighted(…) / bdCrossprod_Weighted(…)     |        ✔︎       |     ✔︎     |     ✔︎    |   ✔︎  |
| Matrix - vector weighted product (XwXt, XtwX)    |           bdwproduct(…) / bdScalarwproduct(…)           |        ✔︎       |           |          |      |
| Matrix vector product                            |                   bdblockmult_vector(…)                   |        ✔︎       |     ✔︎     |     ✔︎    |      |
| Data Normalization (center, scale and both)      |          bdNormalize_Data(…) / bdNormalize_hdf5(…)          |        ✔︎       |     ✔︎     |          |   ✔︎  |
| **Other functions**                              |                                                         |                |           |          |      |
| Vector sum                                       |                   bdparallelVectorSum(…)                  |                |           |     ✔︎    |      |
| Pow(2) vector elements                           |                     bdparallelpow2(…)                     |                |           |     ✔︎    |      |
| **Lineal Algebra Functions**                     |                                                         |                |           |          |      |
| SVD matrix decomposition                         |                 bdSVD(…) / bdSVD_hdf5(…)                |        ✔︎       |     ✔︎     |     ✔︎    |   ✔︎  |
| QR matrix decomposition                          |                         bdQR(…)                         |        ✔︎       |           |          |      |
| Cholesky decomposition                           |                     bdInvCholesky(…)                    |        ✔︎       |           |          |      |
| Matrix Pseudoinverse                             |                      bdpseudoinv(…)                     |        ✔︎       |           |          |      |
| Solve matrix equation (A * X = B )               |                        bdSolve(…)                       |        ✔︎       |           |          |      |
| **Data Analysis**                                |                                                         |                |           |          |      |
| Principal Components Analysis (PCA)              |                      bdPCA_hdf5(…)                      |                |     ✔︎     |          |   ✔︎  |
| MLR-MR (Linear Regression Big Data)              |                       bdlm_paral(…)                       |        ✔︎       |     ✔︎     |     ✔︎    |      |
| **HDF5 data files Utils**                        |                                                         |                |           |          |      |
| Remove rows or columns with hight missing values |                    bdRemovelowdata(…)                   |                |     ✔︎     |          |   ✔︎  |
| Impute missing data                              |                    bdImpute_snps_hdf5(…)                   |                |     ✔︎     |          |   ✔︎  |
| Create hdf5 data file with one dataset inside    |                bdCreate_hdf5_matrix_file(…)               |                |           |          |   ✔︎  |
| Add one dataset in hdf5 data file                |                  bdAdd_hdf5_matrix(…)                  |                |           |          |   ✔︎  |
| Split an hdf5 dataset in small datasets          |                  bdSplit_matrix_hdf5(…)                  |                |           |          |   ✔︎  |
| Reduce multiple datasets applying a function     |                  bdReduce_matrix_hdf5(…)                  |                |           |          |   ✔︎  |
| Merge multiple datasets by rows or columns       |                  bdBind_hdf5(…)                             |                |           |          |   ✔︎  |
| Apply a function to multiple datasets            |                  bdapply_Function_hdf5(…)                  |                |           |          |   ✔︎  |
| Get a list with all datasets inside a group      |                  bdgetDatasetsList_hdf5(…)                  |                |           |          |   ✔︎  |
| Remove one dataset from hdf5 data file           |                  bdRemove_hdf5_element(…)                 |                |           |          |   ✔︎  |
| Import data from text file to HDF5               |                 bdImport_text_to_hdf5(…)                 |                |           |          |   ✔︎  |
| Remove MAF by rows or columns                    |                    bdremove_maf_hdf5(…)                   |                |     ✔︎     |          |   ✔︎  |
| Get MAF by rows or columns                       |                    bdget_maf_hdf5(…)                   |                |     ✔︎     |          |   ✔︎  |



