% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bdCholesky_hdf5}
\alias{bdCholesky_hdf5}
\title{Cholesky Decomposition for HDF5-Stored Matrices}
\usage{
bdCholesky_hdf5(
  filename,
  group,
  dataset,
  outdataset,
  outgroup = NULL,
  fullMatrix = NULL,
  overwrite = NULL,
  threads = NULL,
  elementsBlock = 1000000L
)
}
\arguments{
\item{filename}{Character string. Path to the HDF5 file containing the input matrix.}

\item{group}{Character string. Path to the group containing the input dataset.}

\item{dataset}{Character string. Name of the input dataset to decompose.}

\item{outdataset}{Character string. Name for the output dataset.}

\item{outgroup}{Character string. Optional output group path. If not provided,
results are stored in the input group.}

\item{fullMatrix}{Logical. If TRUE, stores the complete matrix. If FALSE (default),
stores only the lower triangular part to save space.}

\item{overwrite}{Logical. If TRUE, allows overwriting existing results.}

\item{threads}{Integer. Number of threads for parallel computation.}

\item{elementsBlock}{Integer. Maximum number of elements to process in each block
(default = 100,000). For matrices larger than 5000x5000, automatically adjusted
to number of rows or columns * 2.}
}
\value{
No direct return value. Results are written to the HDF5 file in the
specified location with the following structure:
\describe{
\item{L}{The lower triangular Cholesky factor}
}
}
\description{
Computes the Cholesky decomposition of a symmetric positive-definite matrix stored
in an HDF5 file. The Cholesky decomposition factors a matrix A into the product
A = LL' where L is a lower triangular matrix.
}
\details{
The Cholesky decomposition is a specialized factorization for symmetric
positive-definite matrices that provides several advantages:
\itemize{
\item More efficient than LU decomposition for symmetric positive-definite matrices
\item Numerically stable
\item Useful for solving linear systems and computing matrix inverses
\item Important in statistical computing (e.g., for sampling from multivariate normal distributions)
}

This implementation features:
\itemize{
\item Block-based computation for large matrices
\item Optional storage formats (full or triangular)
\item Parallel processing support
\item Memory-efficient block algorithm
}

Mathematical Details:
For a symmetric positive-definite matrix A, the decomposition A = LL' has the following properties:
\itemize{
\item L is lower triangular
\item L has positive diagonal elements
\item L is unique
}

The elements of L are computed using:
\deqn{l_{ii} = \sqrt{a_{ii} - \sum_{k=1}^{i-1} l_{ik}^2}}
\deqn{l_{ji} = \frac{1}{l_{ii}}(a_{ji} - \sum_{k=1}^{i-1} l_{ik}l_{jk})}
}
\examples{
\dontrun{
library(rhdf5)

# Create a symmetric positive-definite matrix
set.seed(1234)
X <- matrix(rnorm(100), 10, 10)
A <- crossprod(X)  # A = X'X is symmetric positive-definite
    
# Save to HDF5
h5createFile("matrix.h5")
h5write(A, "matrix.h5", "data/matrix")
        
# Compute Cholesky decomposition
bdCholesky_hdf5("matrix.h5", "data", "matrix",
                outdataset = "chol",
                outgroup = "decompositions",
                fullMatrix = FALSE)
       
# Verify the decomposition
L <- h5read("matrix.h5", "decompositions/chol")
max(abs(A - L \%*\% t(L)))  # Should be very small
}

}
\references{
\itemize{
\item Golub, G. H., & Van Loan, C. F. (2013). Matrix Computations, 4th Edition.
Johns Hopkins University Press.
\item Higham, N. J. (2009). Cholesky factorization.
Wiley Interdisciplinary Reviews: Computational Statistics, 1(2), 251-254.
}
}
\seealso{
\itemize{
\item \code{\link{bdInvCholesky_hdf5}} for computing inverse using Cholesky decomposition
\item \code{\link{bdSolve_hdf5}} for solving linear systems
}
}
