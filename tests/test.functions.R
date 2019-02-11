library(RcppParallel)
library(rbenchmark) 


#
# xwxt y xtwx
#

n <- 700
X <- matrix(runif(n*700), ncol = n)
D <- DelayedArray(X)
w <- runif(n)

RcppParallel::setThreadOptions(numThreads = 80)

res <- benchmark(X%*%diag(w)%*%t(X),
                 parxwxt(D,w),
                 parxwxt(X,w),
                 t(X)%*%diag(w)%*%X,
                 parxtwx(D,w),
                 parxtwx(X,w),
                 order="relative", replications = c(10))
res[,1:4] 

stopifnot(identical(X%*%diag(w)%*%t(X),parxwxt(X,w) )) # Comprovació xwxt
stopifnot(identical(t(X)%*%diag(w)%*%X,parxtwx(X,w))) # Comprovació xtwx

#
# Producte Matriu x vector (Xy)
#

n <- 10000
X <- matrix(runif(n*10000), ncol = n)
w <- runif(n)
D <- DelayedArray(X)


res <- benchmark(t(X)%*%w,
                 parXy(t(X),w),
                 X%*%w,
                 parXy(X,w),
                 order="relative", replications = c(10))
res[,1:4] 
stopifnot(identical(t(M)%*%Y,parXy(t(M),Y)))


#
# Producte matrius
#

n <- 750
X <- matrix(runif(n*750), ncol = n)
Y <- matrix(runif(n*750), nrow = n)
DX <- DelayedArray(X)
DY <- DelayedArray(Y)

stopifnot(identical(parXYProd(X,Y), X%*%Y))
stopifnot(identical(parXYProd(X,Y, "xty"), t(X)%*%Y))
stopifnot(identical(parXYProd(X,Y, "xyt"), X%*%t(Y)))

res <- benchmark(X%*%Y,
                 t(X)%*%Y,
                 X%*%t(Y),
                 parXYProd(X,Y),
                 parXYProd(X,Y, "xty"),
                 parXYProd(X,Y, "xyt"),
                 parXYProd(DX,DY),
                 parXYProd(DX,DY, "xty"),
                 parXYProd(DX,DY, "xyt"),
                 parXYProdBlock(X,Y),
                 parXYProdBlock(X,Y, "xty"),
                 parXYProdBlock(X,Y, "xyt"),
                 order="relative", replications = c(5))

res[,1:4] 



#
# CrossProd i tCrossProd
#

n <- 750
X <- matrix(runif(n*750), ncol = n)
DX <- DelayedArray(X)

stopifnot(identical(X%*%t(X), partCrossProd(X)))
stopifnot(identical(t(X)%*%X,parCrossProd(X)))

res <- benchmark(X%*%t(X),
                 t(X)%*%X,
                 partCrossProd(X),
                 parCrossProd(X),
                 partCrossProd(DX),
                 parCrossProd(DX),
                 order="relative", replications = c(5))
res[,1:4] 



#
# Descomposició svd
#


n <- 500
p <- 750
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
DM <- DelayedArray(M)

res <- benchmark(descBDsvd <- BDsvd(M, 10,100,FALSE),
                 descBDsvdDelayed <- BDsvd(M, 10,100,FALSE),
                 descsvd <- svd(tcrossprod(M)),
                 desceigen <- eigen(tcrossprod(M)),
                 order="relative", replications = c(5))
res[,1:4] 

dades <- rbind(descBDsvd$`d$`[1:10], descBDsvdDelayed$`d$`[1:10], descsvd$d[1:10], desceigen$values[1:10])
rownames(dades) <- c("descBDsvd","descBDsvdDelayed","descsvd", "desceigen"  )
dades



