% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bdcrossprod}
\alias{bdcrossprod}
\title{Crossproduct and transposed crossproduct of DelayedArray}
\usage{
bdcrossprod(a, transposed = NULL)
}
\arguments{
\item{a}{numerical or Delayed Array matrix}

\item{transposed}{(optional, default = false) boolean indicating if we have to perform a crossproduct (transposed=false) or transposed crossproduct (transposed = true)}
}
\value{
numerical matrix with crossproduct or transposed crossproduct
}
\description{
This function performs a crossproduct or transposed crossproduct of numerical or DelayedArray matrix.
}
\examples{

library(DelayedArray)

n <- 100
p <- 60

X <- matrix(rnorm(n*p), nrow=n, ncol=p)

# without DelayedArray
bdcrossprod(X)
bdcrossprod(X, transposed = TRUE)

# with DelayedArray
XD <- DelayedArray(X)
bdcrossprod(XD)
bdcrossprod(XD, transposed = TRUE)

}
