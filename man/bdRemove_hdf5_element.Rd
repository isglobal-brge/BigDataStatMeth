% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bdRemove_hdf5_element}
\alias{bdRemove_hdf5_element}
\title{Remove element group or dataset from  hdf5 file}
\usage{
bdRemove_hdf5_element(filename, element)
}
\arguments{
\item{filename, }{character array indicating the name of the file to create}

\item{element}{path to element, character array indicating the complete route to the element to be removed (folder or dataset).}
}
\value{
none
}
\description{
Remove group or dataset from  hdf5 file
}
\examples{

matA <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), nrow = 3, byrow = TRUE)
matB <- matrix(c(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,5,3,4,5,2,6,2,3,4,
                   42, 23, 23, 423,1,2), nrow = 3, byrow = TRUE)
                   
bdCreate_hdf5_matrix_file("BasicMatVect.hdf5", matA, "INPUT", "matA")
bdAdd_hdf5_matrix(matB, "BasicMatVect.hdf5", "INPUT", "matB")

bdRemove_hdf5_element("BasicMatVect.hdf5", "INPUT/matA")


# Remove file (used as example)
  if (file.exists("BasicMatVect.hdf5")) {
    # Delete file if it exist
    file.remove("BasicMatVect.hdf5")
  }

}
