library(RcppParallel)
library(rbenchmark) 


n <- 500
p <- 100
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]

RcppParallel::setThreadOptions(numThreads = 40)


res <- benchmark(solprodparal <- Prodparal::LOOE(M,Y),
                 solfin <- LOOE.all(tcrossprod(M),Y),
                 order="relative", replications = c(3))
res[,1:4] 

solprodparal

rbind(solfin[1:10], solprodparal$coef[1:10])



n <- 500
p <- 100
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
Y <- 2.8*M[,1] + 1.1*M[,2] + 0.3*M[,5]


res <- benchmark(solprodparal2 <- Prodparal::LOOE(M,Y),
                 solfin2 <- LOOE.all(tcrossprod(M),Y),
                 order="relative", replications = c(3))
res[,1:4] 

solprodparal2

rbind(solfin2[1:10], solprodparal2$coef[1:10])




## FUNCIONS R (originals)

inversecpp_orig <- function(X, lambda=1, eigen=TRUE,
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

solveEigen_orig <- function(X, Y, lambda){   # X DelayedArray (HDF5), Y un vector
  XX <- tcrossprod(X)
  Ginv <- inversecpp_orig(XX, lambda=lambda)  # Con DelayedArray (HDF5)
  coef <- t(Ginv%*%X)%*%Y
  coef
}


LOOE.i_orig <- function(lambda, Lambda, Q, Y){
  Ginv <- inversecpp_orig(lambda=lambda, Lambda=Lambda,
                          Q=Q, eigen=FALSE)
  cte <- Ginv%*%Y
  ans <- sum((cte/diag(Ginv))^2)
  ans
}



#
# Compute LOOE
#

LOOE_orig <- function(X, Y, nlambdas=100, max.lambda=1, lambdas){
  
  if (missing(lambdas)){
    lambdas <- seq(0.01, max.lambda, length=nlambdas)
  }
  ee <- eigen(X, symmetric = TRUE)
  Lambda <- ee$values
  
  Lambda[Lambda<0] <- 0
  Q <- ee$vectors
  looe <- sapply(lambdas, LOOE.i_orig, Lambda=Lambda, Q=Q, Y=Y)
  
  lambda.min <- lambdas[which.min(looe)]
  
  Ginv <- inversecpp_orig(X, lambda.min)
  
  ans <- list(looe=looe, Ginv=Ginv, lambdas=lambdas,
              lambda.min=lambda.min)
  ans
}

LOOE.all <- function(X, Y, nlambdas=100, max.lambda=1, lambdas){
  sol <- LOOE_orig(tcrossprod(X),Y)
  return(solveEigen_orig(X, as.matrix(Y), lambda=sol$lambda.min))
}
