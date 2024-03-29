% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bdScalarwproduct}
\alias{bdScalarwproduct}
\title{Matrix - Weighted Scalar Multiplication with numerical data}
\usage{
bdScalarwproduct(A, w, op)
}
\arguments{
\item{A}{numerical matrix}

\item{w}{scalar, weight}

\item{op}{string indicating if operation  "Xw" or "wX"}
}
\value{
numerical matrix
}
\description{
This function performs a weighted product of a matrix(X) with a weighted diagonal matrix (w)
}
\examples{

n <- 100
p <- 60

X <- matrix(rnorm(n*p), nrow=n, ncol=p)
w <- 0.75

bdScalarwproduct(X, w,"Xw")
bdScalarwproduct(X, w,"wX")

}
