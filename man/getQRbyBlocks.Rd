% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_bdCCA.R
\name{getQRbyBlocks}
\alias{getQRbyBlocks}
\title{QR by blocks}
\usage{
getQRbyBlocks(strdataset, file, mblocks, center, scale, bcols, overwrt)
}
\arguments{
\item{strdataset}{string, dataset path within the hdf5 data file from which we want to calculate the QR}

\item{file}{string file name where dataset to normalize is stored}

\item{mblocks}{number of blocks in which we want to partition the matrix to perform the calculations}

\item{center, }{boolean, if true, dataset is centered to perform calculus}

\item{scale, }{boolean, if true, dataset is centered to perform calculus}

\item{bcols}{boolean if bcols = TRUE matrix it´s splitted by columns if bcols = FALSE, then matrix or dataset is splitted by rows.}

\item{overwrt, }{boolean, if true, datasets existing inside a file must be overwritten if we are using the same names}
}
\value{
hdf5 data file with CCA results,
}
\description{
This function is an application of the BigDataStatMeth functions to generate new methods. This function perform a QR
from two matrices stored in hdf5 data file. This function applies matrix partitioning, merge bloks to create a full matrix, apply a function to different blocks...
}
\examples{

   print ("Example in vignette")

}
