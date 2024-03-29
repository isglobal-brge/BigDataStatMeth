% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImportData_hdf5.R
\name{bdImportData_hdf5}
\alias{bdImportData_hdf5}
\title{Import data from url or a file}
\usage{
bdImportData_hdf5(
  inFile,
  destFile,
  destGroup,
  destDataset,
  header = TRUE,
  rownames = FALSE,
  overwrite = FALSE,
  sep = NULL
)
}
\arguments{
\item{inFile}{string file name or url with data to import}

\item{destFile}{file name and path to store imported data}

\item{destGroup}{group name to store the dataset}

\item{destDataset}{dataset name to store the input file in hdf5}

\item{header}{(optional) either a logical value indicating whether the column names of x are to be written along with x, or a character vector of column names to be written. See the section on ‘CSV files’ for the meaning of col.names = NA.}

\item{rownames}{(optional) either a logical value indicating whether the row names of x are to be written along with x, or a character vector of row names to be written.}

\item{overwrite}{(optional) either a logical value indicating whether the output file can be overwritten or not.}

\item{sep}{(optional), by default = "\\t". The field separator string. Values within each row of x are separated by this string.}
}
\value{
none value returned, data are stored in a dataset inside an hdf5 data file.
}
\description{
This function download data from an url and decompress data (if needed), then imports the file to hdf5 data file
}
\examples{
   print ("Example in vignette")
   
}
