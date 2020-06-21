#' Close files from hdf5 objects
#' 
#' This function closes the files associated with the hdf5 objects of the operations performed on disk.
#' 
#' @export
#' 
#' @param results object from on disk execution
#' @return void
#' @examples
#' # with numeric matrix
#' 
#' m <- 350
#' n <- 150
#' k <- 50
#' 
#' Abig <- matrix(rnorm(m*k), nrow=m, ncol=k)
#' Bbig <- matrix(rnorm(k*n), nrow=k, ncol=n)
#' 
#' DAbig <- DelayedArray(Abig)
#' DBbig <- DelayedArray(Bbig)
#' 
#' # We consider a big matrix if number of rows or columns are > 500
#' C <- blockmult(DAbig, DBbig, bigmatrix = 1, paral = FALSE)
#' 
#' bdclose(C)
bdclose <- function(results)
{ 

  try (H5Fclose(results$res) )
  
}