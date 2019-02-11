
n <- 500
p <- 100
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
DM <- DelayedArray(M)
Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]

ans1 <- glmnet::cv.glmnet(M, Y) # Estimamos los betas de las 100 variables con `glmnet`:

mm <- which(ans1$lambda==ans1$lambda.min) # el resultado es:

ans1$glmnet.fit$beta[1:10,mm]

# V1         V2         V3         V4         V5         V6         V7         V8         V9        V10 
# 2.3423592  1.5403282  0.0000000  0.0000000 -0.3444285  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000

# ahora veamos nuestra implementación:

sol2 <- LOOE(partCrossProd(M), Y) #@@ sol <- LOOE(tcrossprod(M),Y)
ans2 <- solveEigen(M, Y, lambda=sol2$lambda.min)
ans2[1:10]
sol2$lambda.min
sol2$Ginv[1:10,1:10]


sol3 <- LOOE_orig(partCrossProd(M), Y) #@@ sol <- LOOE(tcrossprod(M),Y)
ans3 <- solveEigen_orig(M, as.matrix(Y), lambda=sol3$lambda.min)
ans3[1:10]
sol3$lambda.min
sol3$Ginv[1:10,1:10]



# [1]  2.399937e+00  1.599961e+00 -1.200198e-06  6.537065e-08  -3.999911e-01 -3.781633e-06  5.707634e-07 -1.428479e-06
# [9]  1.132573e-06 -4.996911e-06

ans3 <- rfunctions::cgls(M,Y, lambda=sol$lambda.min) # amb la función de `rfunctions` obtenim el mateix
ans3$x[1:10]





### PROVES BENCHMAT : 
library(RcppParallel)
library(DelayedArray)
library(rbenchmark) 
library(BiocSingular)

RcppParallel::setThreadOptions(numThreads = 40)

res <- benchmark(sol <- LOOE(partCrossProd(M), Y),
                 sol_orig <- LOOE_orig(partCrossProd(M), Y),
                 order="relative", replications = c(3))
res[,1:4] 

res <- benchmark(solucio <- solveEigen(M, Y, lambda=sol$lambda.min),
                 solucio_orig <- solveEigen_orig(M, as.matrix(Y), lambda=sol_orig$lambda.min),
                 order="relative", replications = c(3))
res[,1:4] 

solucio[1:10]
solucio_orig[1:10]



## Revisar velocitat partCrossProd
res <- benchmark(partCrossProd(X),
                 tcrossprod(X),
                 order="relative", replications = c(3))
res[,1:4] 

#####







inversecpp <- function(X, lambda=1, eigen=TRUE, Lambda, Q){
  if (eigen){
    ee<- runSVD(X, k=length(Y))
    Lambda <- ee$d
    Lambda[Lambda<0] <- 0
    Q <- ee$v
  }
  else
    if(missing(Lambda) | missing(Q))
      stop('SVD results should be provided. \n')
  
  if (lambda == 1)
    W <-  1/Lambda
  else
    W <- 1/(Lambda + lambda)

  Ginv <- parxwxt(Q,W)  
  Ginv
}


solveEigen <- function(X, Y, lambda){   
  XX <- partCrossProd(X)  

  Ginv <- inversecpp(XX, lambda=lambda) 
  coef <- parXYProd(parXYProd(Ginv,X),as.matrix(Y),"xty") 
  coef
}


LOOE.i <- function(lambda, Lambda, Q, Y){
  Ginv <- inversecpp(lambda=lambda, Lambda=Lambda,
                     Q=Q, eigen=FALSE)
  cte <- parXYProd(Ginv,Y) 
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
  
  # Lambda <- ee$values
  ee <- runSVD(X, k=length(Y))
  Lambda <- ee$d
  
  Lambda[Lambda<0] <- 0
  Q <- ee$v
  
  
  looe <- sapply(lambdas, LOOE.i, Lambda=Lambda, Q=Q, Y=Y)

  lambda.min <- lambdas[which.min(looe)]
  
  Ginv <- inversecpp(X, lambda.min)
  
  ans <- list(looe=looe, Ginv=Ginv, lambdas=lambdas,
              lambda.min=lambda.min)
  ans
}





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



