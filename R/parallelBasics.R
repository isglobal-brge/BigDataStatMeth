a <- matrix(rnorm(1000*100), nrow=1000)
parallel_crossprod <- function(X, chunk=100){
lapply(ss)  
  
}


.delayed_crossprod <- function(X, BPPARAM=SerialParam()) 
  # DelayedMatrix crossprod, 1000 rows at a time.
{
  CHUNK <- 1000L
  last <- 0L
  output <- 0
  finish <- nrow(X)
  
  repeat {
    previous <- last + 1L
    last <- min(last + CHUNK, finish)
    block <- as.matrix(X[previous:last,,drop=FALSE])
    output <- output + crossprod(block)
    if (last==finish) { break }
  }
  
  return(output)
}
