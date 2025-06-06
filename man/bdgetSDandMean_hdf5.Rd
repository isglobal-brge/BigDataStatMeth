% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bdgetSDandMean_hdf5}
\alias{bdgetSDandMean_hdf5}
\title{Compute Matrix Standard Deviation and Mean in HDF5}
\usage{
bdgetSDandMean_hdf5(
  filename,
  group,
  dataset,
  sd = NULL,
  mean = NULL,
  byrows = NULL,
  wsize = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{filename}{Character string. Path to the HDF5 file.}

\item{group}{Character string. Path to the group containing the dataset.}

\item{dataset}{Character string. Name of the dataset to analyze.}

\item{sd}{Logical (optional). Whether to compute standard deviation.
Default is TRUE.}

\item{mean}{Logical (optional). Whether to compute mean. Default is TRUE.}

\item{byrows}{Logical (optional). Whether to compute by rows (TRUE) or
columns (FALSE). Default is FALSE.}

\item{wsize}{Integer (optional). Block size for processing. Default is 1000.}

\item{overwrite}{Logical (optional). Whether to overwrite existing results.
Default is FALSE.}
}
\value{
No return value, called for side effects (statistics computation).
Results are stored in the HDF5 file under the 'mean_sd' group with names:
\itemize{
\item 'sd.\code{dataset}' for standard deviation
\item 'mean.\code{dataset}' for mean
}
}
\description{
Computes standard deviation and/or mean statistics for a matrix stored in
HDF5 format, with support for row-wise or column-wise computations.
}
\details{
This function provides efficient statistical computation capabilities with:
\itemize{
\item Computation options:
\itemize{
\item Standard deviation computation
\item Mean computation
\item Row-wise or column-wise processing
}
\item Processing features:
\itemize{
\item Block-based computation
\item Memory-efficient processing
\item Configurable block size
}
\item Implementation features:
\itemize{
\item Safe HDF5 file operations
\item Memory-efficient implementation
\item Comprehensive error handling
}
}

Results are stored in a new group 'mean_sd' within the HDF5 file.
}
\examples{
\dontrun{
library(BigDataStatMeth)

# Create test matrices
set.seed(123)
Y <- matrix(rnorm(100), 10, 10)
X <- matrix(rnorm(10), 10, 1)

# Save to HDF5
bdCreate_hdf5_matrix("test.hdf5", Y, "data", "matrix1",
                     overwriteFile = TRUE)
bdCreate_hdf5_matrix("test.hdf5", X, "data", "vector1",
                     overwriteFile = FALSE)

# Compute statistics
bdgetSDandMean_hdf5(
  filename = "test.hdf5",
  group = "data",
  dataset = "matrix1",
  sd = TRUE,
  mean = TRUE,
  byrows = TRUE,
  wsize = 500
)

# Cleanup
if (file.exists("test.hdf5")) {
  file.remove("test.hdf5")
}
}

}
\references{
\itemize{
\item The HDF Group. (2000-2010). HDF5 User's Guide.
\item Welford, B. P. (1962). Note on a method for calculating corrected
sums of squares and products. Technometrics, 4(3), 419-420.
}
}
\seealso{
\itemize{
\item \code{\link{bdCreate_hdf5_matrix}} for creating HDF5 matrices
}
}
