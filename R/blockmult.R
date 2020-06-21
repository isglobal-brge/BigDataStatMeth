#' Block matrix multiplication with Delayed Array Object
#' 
#' This function performs a block matrix-matrix multiplication with numeric matrix or Delayed Arrays
#' 
#' @export
#' 
#' @param a a double matrix.
#' @param b a double matrix.
#' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
#' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
#' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
#' @param bigmatrix (optiona, default = 5000) maximum number of rows or columns to consider as big matrix and work with
#' hdf5 files, by default a matrix with more than 5000 rows or files is considered big matrix and computation is made in disk 
#' @param mixblock_size (optiona, default = 128), only if we are working with big matrix and parallel computation = true. 
#' Block size for mixed computation in big matrix parallel. Size of the block to be used to perform parallelized memory 
#' memory of the block read from the disk being processed.
#' @param outfile (optional) file name to work with hdf5 if we are working with big matrix in disk.
#' @return numerical matrix
#' @examples
#' # with numeric matrix
#' m <- 500
#' k <- 1500
#' n <- 400
#' A <- matrix(rnorm(n*p), nrow=n, ncol=k)
#' B <- matrix(rnorm(n*p), nrow=k, ncol=n)
#' 
#' blockmult(A, B, 128, TRUE)
#' 
#' # with Delaeyd Array
#' AD <- DelayedArray(A)
#' BD <- DelayedArray(B)
#' 
#' blockmult( AD, BD, 128, TRUE)
blockmult <- function( a, b, block_size = 128, paral = TRUE, threads = NULL, bigmatrix = 5000, mixblock_size = 128, outfile = "./tmp_blockmult.hdf5")
{ 
  
  res <- .Call('_BigDataStatMeth_Bblockmult', PACKAGE = 'BigDataStatMeth', a, b, block_size, paral, threads, bigmatrix, outfile, mixblock_size)

  if (res$filename == '')
    return (res$matrix)
  else
    return( list("res" = H5Fopen(res$filename), "file" = res$filename) )

}
