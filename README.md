# Big Data Statistical Methods (BigDataStatMeth)

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/BigDataStatMeth)](https://CRAN.R-project.org/package=BigDataStatMeth)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/BigDataStatMeth)](https://CRAN.R-project.org/package=BigDataStatMeth)
<!-- badges: end -->

This package implements basic Algebra methods using parallel algorithms to be used in big data problems such as omic data analyses. The methods will accept objects as those designed to deal with big matrices in the Matrix package. The implemented methods will include:

- Algebra

     - Matrix multiplication (by blocks, parallel, ...)
     - Cross-product and transpose cross-product
     - Matrix vector multiplication
     - Apply vector to each matrix column/row (+, -, *, /)
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

### From CRAN (recommended)
```{r, installCRAN, eval=FALSE}
install.packages("BigDataStatMeth")
```

### From GitHub (development version)
```{r, install, eval=FALSE}
# Install devtools and load library (if not previously installed)
install.packages("devtools") 
library(devtools)

# Install BigDataStatMeth and load package
install_github("isglobal-brge/BigDataStatMeth@HEAD")
library(BigDataStatMeth)

```

## Implemented Functions
 
|                                                  |                        R Function                       | By blocks | Parallel | HDF5 |
|--------------------------------------------------|:-------------------------------------------------------:|:---------:|----------|------|
| **Basic functions with vectors and matrices**    |                                                         |           |          |      |
| Matrix Product                                   |             bdblockmult(…) / bdblockmult_hdf5(…)            |     ✔︎     |     ✔︎    |   ✔︎  |
| Matrix/Matrix and Matrix/Array Summatory         |                    blockSum_hdf5(…)                     |           |          |   ✔︎  |
| Sparse Matrix Product                            |      bdblockmult_sparse(…) / bdblockmult_sparse_hdf5(…)     |           |          |   ✔︎  |
| Matrix product with its transpose                | bdCrossprod(…) / bdtCrossprod(…) / bdCrossprod_hdf5(…) / bdtCrossprod_hdf5(…)  |     ✔︎     |     ✔︎    |   ✔︎  |
| Matrix - Matrix weighted product (XWXt, XtWX)    |      bdtCrossprod_Weighted(…) / bdCrossprod_Weighted(…)  |     ✔︎     |     ✔︎    |   ✔︎  |
| Matrix - vector weighted product (XwXt, XtwX)    |           bdwproduct(…) / bdScalarwproduct(…) / bdWeightedProduct_hdf5(...)        |     ✔︎     |          |   ✔︎  |
| Matrix vector product                            |                   bdblockmult_vector(…)                  |     ✔︎     |     ✔︎    |      |
| Apply calculus to each Matrix col/row (+, -, *, /)|           bdcomputeMatrixVector_hdf5(...)               |      ✔︎     |     ✔︎      |   ✔︎  |
| Data Normalization (center, scale and both)      |          bdNormalize_Data(…) / bdNormalize_hdf5(…)       |     ✔︎     |          |   ✔︎  |
| **Other functions**                              |                                                         |           |          |      |
| Vector sum                                       |                   bdparallelVectorSum(…)                |           |     ✔︎    |      |
| Get the diagonal from a hdf5 matrix dataset      |                    bdgetDiagonal_hdf5(…)                |           |     ︎      |   ✔︎  |
| Write the diagonal to a hdf5 matrix dataset      |                    bdWriteDiagonal_hdf5(…)                |           |     ︎      |   ✔︎  |
| Duplicate the lower/upper triangular hdf5 dataset|            bdWriteOppsiteTriangularMatrix_hdf5(…)       |           |     ︎      |   ✔︎  |
| Get mean and sd from a hdf5 dataset by cols/rows)|                    bdgetSDandMean_hdf5(…)                |     ✔︎     |          |   ✔︎  |
| Pow(2) vector elements                           |                     bdparallelpow2(…)                   |           |     ✔︎    |      |
| Sum two vectors                                  |                     bdparallelVectorSum(…)              |           |     ✔︎    |      |
| **Lineal Algebra Functions**                     |                                                         |           |          |      |
| SVD matrix decomposition                         |                 bdSVD(…) / bdSVD_hdf5(…)                 |     ✔︎     |     ✔︎    |   ✔︎  |
| QR matrix decomposition                          |                         bdQR(…)                          |           |          |      |
| Cholesky decomposition                           |            bdInvCholesky(…) / bdInvCholesky_hdf5(...)                     |          |   ✔︎  |
| Matrix Pseudoinverse                             |                      bdpseudoinv(…)                      |           |          |      |
| Solve matrix equation (A * X = B )               |                        bdSolve(…)                        |           |          |      |
| **Data Analysis**                                |                                                         |           |          |      |
| Principal Components Analysis (PCA)              |                      bdPCA_hdf5(…)                      |     ✔︎     |          |   ✔︎  |
| MLR-MR (Linear Regression Big Data)              |                       bdlm_paral(…)                      |     ✔︎     |     ✔︎    |      |
| **HDF5 data files Utils**                        |                                                         |           |          |      |
| Remove rows or columns with high missing values  |                    bdRemovelowdata(…)                   |     ✔︎     |          |   ✔︎  |
| Impute missing data                              |                    bdImpute_snps_hdf5(…)                |     ✔︎     |          |   ✔︎  |
| Create hdf5 data file with one dataset inside    |                bdCreate_hdf5_matrix_file(…)             |           |          |   ✔︎  |
| Add one dataset in hdf5 data file                |                  bdAdd_hdf5_matrix(…)                   |           |          |   ✔︎  |
| Split an hdf5 dataset in small datasets          |                  bdSplit_matrix_hdf5(…)                 |           |          |   ✔︎  |
| Sort an hdf5 dataset                             |                  bdSort_hdf5_dataset(…)                 |           |          |   ✔︎  |
| Reduce multiple datasets applying a function     |                  bdReduce_matrix_hdf5(…)                |           |          |   ✔︎  |
| Merge multiple datasets by rows or columns       |                  bdBind_hdf5(…)                         |           |          |   ✔︎  |
| Apply a function to multiple datasets            |                  bdapply_Function_hdf5(…)               |           |          |   ✔︎  |
| Get a list with all datasets inside a group      |                  bdgetDatasetsList_hdf5(…)              |           |          |   ✔︎  |
| Remove one dataset from hdf5 data file           |                  bdRemove_hdf5_element(…)               |           |          |   ✔︎  |
| Create a link to other dataset inside a hdf5 file|                  bdCreateLink_hdf5(…)                   |           |          |   ✔︎  |
| Create a group inside a hdf5 datafile            |                  bdCreateGroup_hdf5(…)                  |           |          |   ✔︎  |
| Create empty dataset inside a hdf5 datafile      |                  bdCreateEmptyDataset_hdf5(…)                   |           |          |   ✔︎  |
| Check if dataset exists in a hdf5 datafile       |                  bdExists_hdf5_element(…)                   |           |          |   ✔︎  |
| Import data from text file or url to HDF5        |                  bdImportData_hdf5(…)                   |     ✔︎     |          |   ✔︎  |
| Write dimnames inside the hdf5 datafile          |                  bdWriteDimnames_hdf5(…)                   |           |          |   ✔︎  |
| **Develop methods with hdf5 - examples**         |                                                         |          |          |      |
| Perform QR decomposition by blocks in hdf5       |                  getQRbyBlocks(…)                       |                 |          |   ✔︎  |
| Perform CCA by blocks in hdf5                    |                  bdCCA_hdf5(…)                          |     ✔︎     |          |   ✔︎  |



