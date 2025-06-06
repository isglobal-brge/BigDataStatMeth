% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bdblockSubstract_hdf5}
\alias{bdblockSubstract_hdf5}
\title{HDF5 dataset subtraction}
\usage{
bdblockSubstract_hdf5(
  filename,
  group,
  A,
  B,
  groupB = NULL,
  block_size = NULL,
  paral = NULL,
  threads = NULL,
  outgroup = NULL,
  outdataset = NULL,
  overwrite = NULL
)
}
\arguments{
\item{filename}{String indicating the HDF5 file path}

\item{group}{String indicating the group containing matrix A}

\item{A}{String specifying the dataset name for matrix A}

\item{B}{String specifying the dataset name for matrix B}

\item{groupB}{Optional string indicating group containing matrix B.
If NULL, uses same group as A}

\item{block_size}{Optional integer specifying block size for processing.
If NULL, automatically determined based on matrix dimensions}

\item{paral}{Optional boolean indicating whether to use parallel processing.
Default is false}

\item{threads}{Optional integer specifying number of threads for parallel processing.
If NULL, uses maximum available threads}

\item{outgroup}{Optional string specifying output group.
Default is "OUTPUT"}

\item{outdataset}{Optional string specifying output dataset name.
Default is "A_-_B"}

\item{overwrite}{Optional boolean indicating whether to overwrite existing datasets.
Default is false}
}
\value{
Modifies the HDF5 file in place, adding the subtraction result
}
\description{
Performs optimized block-wise subtraction between two datasets stored in HDF5
format. Supports both matrix-matrix and matrix-vector operations with
memory-efficient block processing.
}
\details{
The function implements optimized subtraction through:

Operation modes:
\itemize{
\item Matrix-matrix subtraction (A - B)
\item Matrix-vector subtraction
\item Vector-matrix subtraction
}

Block processing:
\itemize{
\item Automatic block size selection
\item Memory-efficient operations
\item Parallel computation support
}

Block size optimization based on:
\itemize{
\item Matrix dimensions
\item Available memory
\item Operation type (matrix/vector)
}

Error handling:
\itemize{
\item Dimension validation
\item Resource management
\item Exception handling
}
}
\examples{
\dontrun{
library(BigDataStatMeth)

# Create test matrices
N <- 1500
M <- 1500
set.seed(555)
a <- matrix(rnorm(N*M), N, M)
b <- matrix(rnorm(N*M), N, M)

# Save to HDF5
bdCreate_hdf5_matrix("test.hdf5", a, "data", "A",
                     overwriteFile = TRUE)
bdCreate_hdf5_matrix("test.hdf5", b, "data", "B",
                     overwriteFile = FALSE)

# Perform subtraction
bdblockSubstract_hdf5("test.hdf5", "data", "A", "B",
                      outgroup = "results",
                      outdataset = "diff",
                      block_size = 1024,
                      paral = TRUE)
}

}
