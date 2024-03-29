#' Crossproduct 
#' 
#' This function performs a crossproduct or transposed crossproduct of numerical matrix.
#' 
#' @export
#' 
#' @param A numerical matrix
#' @param B optional, numerical matrix
#' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
#' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
#' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
#' @return numerical matrix with crossproduct
#' @examples
#' 
#' n <- 100
#' p <- 60
#' 
#' X <- matrix(rnorm(n*p), nrow=n, ncol=p)
#' res <- bdCrossprod(X)
#' 
#' all.equal(crossprod(X), res)
#' 
#' n <- 100
#' p <- 100
#' 
#' Y <- matrix(rnorm(n*p), nrow=n)
#' 
#' # With two matrices
#' res <- bdCrossprod(X,Y)
#' 
bdCrossprod <- function( A, B = NULL, block_size = 256, paral = TRUE, threads = NULL)
{ 
  
  res <- .Call('_BigDataStatMeth_bdCrossprod_generic', PACKAGE = 'BigDataStatMeth', A, B, FALSE, block_size, paral, threads)
  
  return (res)

}
