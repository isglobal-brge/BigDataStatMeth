#' Block matrix-vector multiplication with Delayed Array Object
#' 
#' This function performs a block matrix-vector multiplication with R-Objects or Delayed Arrays
#' 
#' @export
#' 
#' @param A a double matrix.
#' @param b a double vector or array.
#' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
#' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
#' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
#' @param outfile (optional) file name to work with hdf5 if we are working with big matrix in disk.
#' @param onmemory (optional) if onmemory = TRUE the multiplication is forced to execute in memory
#' @return numerical matrix
#' @examples
#' # with numeric matrix
#' 
#' library(DelayedArray)
#' 
#' k <- 100
#' n <- 400
#' A <- matrix(rnorm(n*k), nrow=n, ncol=k)
#' B <- sample(1:100,50, replace = TRUE);
#' 
#' res <- bdblockmult_vector(A, B, 128, TRUE)
#' 
#' 
bdblockmult_vector <- function( A, b, block_size = 128, paral = TRUE, threads = NULL, outfile = "tmp_blockmult_vector.hdf5", onmemory = FALSE)
{ 
   
   if(is.array(b) || is.vector(b)) {
      b <- matrix(b);
   } else if( is.matrix(b) && dim(b)[1]>1 && dim(b)[2]>1) {
      stop("b should be a vector not a matrix")
   } 
   
   res <- .Call('_BigDataStatMeth_blockmult', PACKAGE = 'BigDataStatMeth', A, b, block_size, paral, threads, 999999999, 999999999999, outfile, onmemory)
   
   if (res$filename == '' | is.null(res$filename))
      return (res$matrix)
   else
      return( list( "file" = res$filename, "dataset" = res$dataset) )
   
}
