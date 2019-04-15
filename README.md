# BigDataStatMeth

This package implements basic Algebra methods using parallel algorithms to be used in big data problems such as omic data analyses. The functions will consider as input an object of class `DelayedArray` which is an array-like object (typically an on-disk object) specifically designed for Bioconductor. The methods will also accept other type of objects as those designed to deal with big matrices in the `Matrix` pacKage. The implemented methods will include:

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


