#' Block matrix multiplication with hdf5 datasets
#'
#' This function performs a block matrix-matrix multiplication with mattrix stored in HDF5 file
#'
#' @export
#' @param filename string file name where dataset to normalize is stored
#' @param group string with the group name where matrix is stored inside HDF5 file
#' @param a a double matrix.
#' @param b a double matrix.
#' @param groupB, string, (optional) group name where dataset b is stored
#' @param block_size (optional, defalut = 128) block size to make matrix 
#' multiplication, if `block_size = 1` no block size is applied 
#' (size 1 = 1 element per block)
#' @param paral, (optional, default = FALSE) paral = true --> TO BE IMPLEMENTED
#' @param threads (optional) only if bparal = true, number of concurrent threads 
#' in parallelization if threads is null then threads =  maximum number of 
#' threads available
#' @param outgroup (optional) string with group name where we want to store the 
#' result matrix
#' @param outdataset (optional) string with dataset name where we want to store 
#' the results
#' @return A list with an HDF5 object with numerical matrix and  HDF5 file name 
#' with results
#' \itemize{
#'  \item{"file"}{string with hdf5 file name with result and original matrixs}
#'  \item{"dataset"}{string complete path inside hdf5 file where results are stored.}
#' }
#' @examples
#'
#' library(BigDataStatMeth)
#' library(rhdf5)
#' 
#' 
#' # Sum matrix twice
#' 
#' # Prepare data and functions
#' X <- matrix(rnorm(200*10), nrow = 10, ncol = 200)
#' diag(X) <- 0.5
#' 
#' # Create hdf5 data file with  two X matrices ( named X and Y)
#' bdCreate_hdf5_matrix_file("test_file3.hdf5", X, "data", "X", force = TRUE)
#' bdAdd_hdf5_matrix(X, "test_file3.hdf5", "data", "Y", force = TRUE)
#' 
#' # Update diagonal
#' diagonal <- bdblockSum_hdf5( filename = "test_file3.hdf5", 
#'                    group = "data", a = "X", b = "Y", groupB = "data", 
#'                    block_size = 128, outgroup = "tmp", outdataset = "Z")
#' 
#'
#' # Sum matrix with an Array
#' K <- 1
#' N <- 300
#' 
#' M <- 300
#' L <- 2500
#' 
#' 
#' A <- rnorm(K*N)
#' B <- matrix(rnorm(M*L), nrow = M, ncol = L)
#' 
#' # Create hdf5 data file with A and B matrices
#' bdCreate_hdf5_matrix_file("test_file4.hdf5", as.matrix(A), "data", "A", force = TRUE)
#' bdAdd_hdf5_matrix(B, "test_file4.hdf5", "data", "B", force = TRUE)
#' 
#' results <- bdblockSum_hdf5( filename = "test_file4.hdf5", 
#'             group = "data", a = "A", b = "B", groupB = "data", 
#'             block_size = 128, outgroup = "tmp", outdataset = "Z")
#'

bdblockSum_hdf5 <- function( filename, group, a, b, groupB = NULL, block_size = 2048, paral = FALSE, threads = NULL,
                       outgroup = "OUTPUT", outdataset = NULL)
{

  res <- .Call('_BigDataStatMeth_blockSum_hdf5', PACKAGE = 'BigDataStatMeth', filename, group, a, b, groupB, block_size, paral, threads, outgroup, outdataset)

  if (res$result == 0)
    return( list( "file" = res$filename, "dataset" = res$dataset) )
  else
    warning("Error performing matrix sum")

}
