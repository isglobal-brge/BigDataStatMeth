library(microbenchmark)
library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

N <- 1000
M <- 1000
nc <-  4
ntimes <- 2

set.seed(555)
mat <- matrix( rnorm( N*M, mean=0, sd=10), N, M) 
set.seed(444)
mat2 <- matrix( rnorm( N*M, mean=0, sd=10), N, M) 


# Matrix - CrossProd


R <-  crossprod(mat)
R2 <- crossprod(mat, mat2)

bd1 <- bdCrossprod(mat)
bd2 <- bdCrossprod(mat, mat2)

all.equal(R, bd1)
all.equal(R2, bd2)


bd_p1 <- bdCrossprod(mat, block_size = 500, paral = TRUE, threads = 2)
bd_p2 <- bdCrossprod(mat, mat2, block_size = 500, paral = TRUE, threads = 2)

all.equal(R, bd_p1)
all.equal(R2, bd_p2)


## Si tot funciona OK ==>

times <- microbenchmark::microbenchmark(
    R =  crossprod(mat),
    R2 = crossprod(mat, mat2),
    bd1 = bdCrossprod(mat),
    bd2 = bdCrossprod(mat, mat2),
    bd_p1 = bdCrossprod(mat, block_size = 500, paral = TRUE, threads = 2),
    bd_p2 = bdCrossprod(mat, mat2, block_size = 500, paral = TRUE, threads = 2),
    times = ntimes, unit = "s")
times



# Matrix - tCrossProd