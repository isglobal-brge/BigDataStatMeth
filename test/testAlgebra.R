
a <- matrix(seq(1:25), nrow = 5)

f <- matrix(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), nrow = 1)
f <- matrix(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), nrow = 3)

b <-  as.vector(c(1,2,3,4,5))
c <-  c(1,2,3,4,5)


library("BigDataStatMeth")
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")
bdCreate_hdf5_matrix_file("test_temp.hdf5", f, 
                          "pepet", "datasetpepet",
                          FALSE, TRUE)





# Test blockmult

library("BigDataStatMeth")
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

N = 20
M = 10

set.seed(555)
a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
b <- matrix( rnorm( N*M, mean=0, sd=1), M, N) 
v <- rnorm( N, mean=0, sd=1)

Y <- matrix(rnorm(100), 10, 10)
X <- matrix(rnorm(10), 10, 1)

X <- matrix(rnorm(100), 10, 10)
X.1 <- X
X[lower.tri(X)] <- 0
X <- X.1
X[upper.tri(X)] <- 0



# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = a, group = "pepet", 
                     dataset = "datasetpepet",
                     transp = FALSE,
                     overwriteFile = TRUE, 
                     overwriteDataset = FALSE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(b), 
                     group = "pepet", 
                     dataset = "tdatasetpepet",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(v), 
                     group = "pepet", 
                     dataset = "vdatasetpepet",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)


bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = X, 
                     group = "pepet", 
                     dataset = "triangular",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)


# Matrix-Vector Sum
# ----------------------

# Sum two matrix
bdblockSum_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", force = TRUE ) 
bdblockSum_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 5, force = TRUE ) 

# Sum matrix + Vector
bdblockSum_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "vdatasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", force = TRUE,block_size = 4 ) 
bdblockSum_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "vdatasetpepet", outgroup = "results", outdataset = "res", force = TRUE ) 


# Matrix-Vector Substract
# ----------------------

# Sum two matrix
bdblockSubstract_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", force = TRUE ) 
bdblockSubstract_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 5, force = TRUE ) 

# Sum matrix + Vector
bdblockSubstract_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "vdatasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", force = TRUE,block_size = 4 ) 
bdblockSubstract_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "vdatasetpepet", outgroup = "results", outdataset = "res", force = TRUE ) 



# Matrix Multiplication
# ----------------------

bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "resSparse", force = TRUE )
# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "resSparse", force = TRUE, paral = TRUE, threads = 4 )


# Matrix Multiplication - Sparse
# -----------------------------

bdbdblockmult_sparse_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "resSparse", force = TRUE )
bdbdblockmult_sparse_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "resSparse", force = TRUE, block_size = 10  )

## Parallel not working (at the moment --> need to be reviewed)
# .. TO BE REVIEWED .. # BigDataStatMeth::bdbdblockmult_sparse_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outdataset = "prod_sparse_parall", force = TRUE, paral = TRUE, threads = 4 )

# Matrix Crossprodduct
# ---------------------

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", force = TRUE ) # Crossprod A %*% t(A)
bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", block_size = 1024, force = TRUE ) # Crossprod A %*% t(A)
bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", block_size = 4, force = TRUE ) # Crossprod A %*% t(A)
bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", force = TRUE ) # Crossprod A %*% t(B)


# Matrix tCrossprodduct
# ---------------------

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdtCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", force = TRUE ) # tCrossprod t(A) %*% A
bdtCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 1024, force = TRUE ) # tCrossprod t(A) %*% A
bdtCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 4, force = TRUE ) # tCrossprod t(A) %*% A
bdtCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "res", force = TRUE ) # tCrossprod t(A) %*% B


# Matrix Diagonals - Get and Set
# -------------------------------

diagonal <- bdgetDiagonal_hdf5("test_temp.hdf5", "pepet", "datasetpepet")
bdWriteDiagonal_hdf5(diagonal + 2, "test_temp.hdf5", "pepet", "datasetpepet")

# Matrix Triangular - Write opposite 
# -----------------------------------
# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdWriteOppsiteTriangularMatrix_hdf5(filename = "test_temp.hdf5", group = "pepet", dataset = "triangular", copytolower = TRUE, elementsBlock = 4)
bdWriteOppsiteTriangularMatrix_hdf5(filename = "test_temp.hdf5", group = "pepet", dataset = "triangular", copytolower = FALSE, elementsBlock = 4)


# Matrix Vector - Calculus
# -----------------------------------

bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "*", byrows = TRUE,  paral=FALSE, threads = 1, force = TRUE)






##########
##########      CHECKING IN PROGRESS..... 
##########

devtools::reload(pkgload::inst("BigDataStatMeth"))

diagonal2 <- diagonal + 2
bdWriteDiagonal_hdf5(diagonal2, "test_temp.hdf5", "pepet", "datasetpepet")

file <- "test_temp.hdf5"
dataset <- "pepet/datasetpepet"
resr <- diag(a)
res <-  h5read(file,dataset)
# res[1:5,1:5]
# resr[1:5,1:5]
all.equal( round( res, 5), round( resr, 5))


all.equal( round( diag(res), 5), ound( diagonal2, 5) )

##########      END CHECKING IN PROGRESS....







library(rbenchmark)

res <- rbenchmark::benchmark( "R" = (a%*%b),
                              "BDSM" =  bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", force = TRUE ),
                              "BDSM_par2" =  bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", force = TRUE, paral = TRUE, threads = 2 ),
                              "BDSM_par4" =  bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", force = TRUE, paral = TRUE, threads = 4 ),
                              replications = 5,
                              columns = c("test", "replications", "elapsed",
                                          "relative", "user.self", "sys.self"))
print(res)

bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet",outdataset = "TOTAMAT")
bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet",outdataset = "resultatC", block_size = 1000, paral = T, threads = 4 )

f