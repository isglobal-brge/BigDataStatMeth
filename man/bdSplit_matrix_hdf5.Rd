% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bdSplit_matrix_hdf5}
\alias{bdSplit_matrix_hdf5}
\title{Split hdf5 dataset}
\usage{
bdSplit_matrix_hdf5(
  filename,
  group,
  dataset,
  outgroup = NULL,
  outdataset = NULL,
  nblocks = NULL,
  blocksize = NULL,
  bycols = TRUE,
  force = FALSE
)
}
\arguments{
\item{filename, }{character array indicating the name of the file to create}

\item{group, }{character array indicating the input group where the data set to be imputed is.}

\item{dataset, }{character array indicating the input dataset to be imputed}

\item{outgroup, }{optional character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group.}

\item{outdataset, }{optional character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten.}

\item{nblocks, }{integer number of blocks in which we want to split the data}

\item{blocksize, }{integer, number of elements in each block}

\item{bycols, }{boolean by default = true, true indicates that the imputation will be done by columns, otherwise, the imputation will be done by rows}

\item{force, }{boolean if true, previous results in same location inside hdf5 will be overwritten.}
}
\value{
Original hdf5 data file with imputed data
}
\description{
Split hdf5 dataset by rows or columns and store splitted submatrices inside hdf5 file.
}
