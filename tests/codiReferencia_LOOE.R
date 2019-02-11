# Y = a + bX
# A'A = X
# Solve the linear system (A'A + lambda * I)x = Ab using
# eigen decomposition
#
# example: X data Y outcome
#    M <- tcrossprod(X)
#    LOOE(M, Y)


devtools::install_github("jaredhuling/rfunctions")
library(rfunctions)


inversecpp <- function(X, lambda=1, eigen=TRUE,
                       Lambda, Q){
  if (eigen){
    ee <- eigen(X, symmetric = TRUE)  # mirar la librería BiocSingular
    Lambda <- ee$values
    Lambda[Lambda<0] <- 0
    Q <- ee$vectors
  }
  else
    if(missing(Lambda) | missing(Q))
      stop('SVD results should be provided. \n')
  
  if (lambda == 1)
    W <-  1/Lambda
  else
    W <- 1/(Lambda + lambda)
  
  Ginv <- rfunctions::crossprodcpp(t(Q), W)  # implementar xwxt  (en la librería rfunctions está xtwx) 
  # Q%*%diag(W)%*%t(Q)
  Ginv
}

solveEigen <- function(X, Y, lambda){   # X DelayedArray (HDF5), Y un vector
  XX <- tcrossprod(X)
  Ginv <- inversecpp(XX, lambda=lambda)  # Con DelayedArray (HDF5)
  coef <- t(Ginv%*%X)%*%Y
  coef
}


LOOE.i <- function(lambda, Lambda, Q, Y){
  Ginv <- inversecpp(lambda=lambda, Lambda=Lambda,
                     Q=Q, eigen=FALSE)
  cte <- Ginv%*%Y
  ans <- sum((cte/diag(Ginv))^2)
  ans
}



#
# Compute LOOE
#

LOOE <- function(X, Y, nlambdas=100, max.lambda=1, lambdas){
  
  if (missing(lambdas)){
    lambdas <- seq(0.01, max.lambda, length=nlambdas)
  }
  ee <- eigen(X, symmetric = TRUE)
  Lambda <- ee$values
  Lambda[Lambda<0] <- 0
  Q <- ee$vectors
  looe <- sapply(lambdas, LOOE.i, Lambda=Lambda, Q=Q, Y=Y)
  
  lambda.min <- lambdas[which.min(looe)]
  
  Ginv <- inversecpp(X, lambda.min)
  
  ans <- list(looe=looe, Ginv=Ginv, lambdas=lambdas,
              lambda.min=lambda.min)
  ans
}




X <- matrix(rnorm(10000), ncol=50)
M <- crossprod(X)
Y <- matrix(rnorm(10000), ncol=200)
dim(X);dim(Y)

LOOE(M, Y)

