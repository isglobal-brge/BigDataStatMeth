#' Linear regression using MLR-MR algorithm
#' 
#' Linear regression for Big Data using MLR-MR algorithm to perform lm regression with big data :  https://isglobal-brge.github.io/Master_Modelling/dealing-with-big-data-in-r.html#linear-regression-for-big-data
#' 
#' @export
#' 
#' @param Y, numerical matrix column with response variable
#' @param model, numerical matrix with paired observations of the predictor variable X
#' @param blocks, integer with number of blocks we want to split matrix if null matrix is splited in blocks as maximum of 1000 variables per block
#' @param threads, threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
#' @return Lineal regression coefficients 
#' @examples
#' # with numeric matrix
#' 
#' library(BigDataStatMeth)
#' data(mtcars)
#' 
#' Y <- mtcars$mpg
#' X <- model.matrix(~ wt + cyl, data=mtcars)
#' m <- 7
#' 
#' 
#' res <- bdlm_paral( X, Y, m, 1)
#' res
#' 
bdlm_paral <- function( Y, model, blocks, threads = NULL)
{ 
   
   res <- .Call('_BigDataStatMeth_bdMLR_MR', PACKAGE = 'BigDataStatMeth', model, Y, blocks, threads)
   rownames(res) <- colnames(model)
   colnames(res) <- c("coef.")
   return(res)
   
}
