
# Test SVD - InvCholesky

library("BigDataStatMeth")
library("microbenchmark")
library(rhdf5)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("~/PhD/dummy")

N = 5000
M = 300

# N = 300
# M = 5000

set.seed(555)
Y <- matrix(rnorm(N*M), N, M)
X <- matrix(rnorm(10), 10, 1)
# Ycp <- crossprod(Y)
# X <- c(1,2,3,4,5,6,7,8,9,10)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Y, group = "data", dataset = "matrix",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)
# devtools::reload(pkgload::inst("BigDataStatMeth"))
microbenchmark( svdR = svd(scale(Y)),
                SVD_block_q1k2 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL), 
                times = 10 )
res <- svd(Y)
file <- "test_temp.hdf5"
dataset.d <- "/SVD/matrix/d"
dataset.u <- "/SVD/matrix/u"

res.d <-  h5read(file, dataset.d)
res.u <-  h5read(file, dataset.u)




microbenchmark( SVD_block_full = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "full", dataset = "matrix", bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL), times = 3 )
microbenchmark( SVD_block_q1k2 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL), times = 3 )
microbenchmark( SVD_block_q1k4 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 4, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL), times = 3 )
microbenchmark( SVD_block_q4k1 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 1, q = 4, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL), times = 3 )


bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = X,  group = "data",  dataset = "vector",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)



# Matrix - SVD 
# -----------------------------------

bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = FALSE, bscale = FALSE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL); resr <- svd(Y)
bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = TRUE, bscale = FALSE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);  resr <- svd(scale(Y, center = TRUE, scale = FALSE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = FALSE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);  resr <- svd(scale(Y, center = FALSE, scale = TRUE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))
# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", method = "blocks", k = 4, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", method = "blocks", k = 2, q = 2, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))


# devtools::reload(pkgload::inst("BigDataStatMeth"))
times <- microbenchmark::microbenchmark( SVD_R = svd(scale(Y, center = TRUE, scale = TRUE)),
                                         SVD_block_full = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "full", dataset = "matrix", k = 4, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         # SVD_full = bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         SVD_block_q1k2 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         SVD_block_q1k4 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 4, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         SVD_block_q4k1 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 1, q = 4, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         # SVD_block_q1k2 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 1, q = 4, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         # SVD_block_q2k2 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 2, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         # SVD_block_q2k4 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 4, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         # SVD_block_q3k2 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 3, q = 2, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         # SVD_block_q4k2 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 4, q = 2, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         times = 1 )
times



# Matrix - SVD - by Blocks
# -----------------------------------

# N bigger than M (Vertical Matrix)
N = 5220 ; M = 22 ; set.seed(555)
Y <- matrix(rnorm(N*M), N, M)
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", object = Y, group = "data", dataset = "matrix", transp = FALSE, overwriteFile = TRUE, overwriteDataset = TRUE, unlimited = FALSE)

bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = FALSE, bscale = FALSE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL); resr <- svd(scale(Y, center = FALSE, scale = FALSE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = TRUE, bscale = FALSE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);  resr <- svd(scale(Y, center = TRUE, scale = FALSE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = FALSE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);  resr <- svd(scale(Y, center = FALSE, scale = TRUE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 1, q = 4, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 2, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))


times <- microbenchmark::microbenchmark( SVD_R = svd(scale(Y, center = TRUE, scale = TRUE)),
                                         SVD_full = bdSVD_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", k = 1, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         SVD_block_q1k4 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 1, q = 4, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         SVD_block_q2k2 = bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 2, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL),
                                         times = 3 )

times





# M bigger than N (Horizontal Matrix)
N = 220 ; M = 5220 ; set.seed(555)
Y <- matrix(rnorm(N*M), N, M)
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", object = Y, group = "data", dataset = "matrix", transp = FALSE, overwriteFile = TRUE, overwriteDataset = TRUE, unlimited = FALSE)

bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 1, bcenter = FALSE, bscale = FALSE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);  resr <- svd(Y)
bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 2, bcenter = TRUE, bscale = FALSE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);  resr <- svd(scale(Y, center = TRUE, scale = FALSE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 1, q = 3, bcenter = FALSE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);  resr <- svd(scale(Y, center = FALSE, scale = TRUE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 1, q = 1, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 1, q = 4, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))
bdSVD_hdf5( "test_temp.hdf5", group = "data", method = "blocks", dataset = "matrix", k = 2, q = 2, bcenter = TRUE, bscale = TRUE, rankthreshold = 0.0, overwrite  = TRUE, threads = NULL);   resr <- svd(scale(Y, center = TRUE, scale = TRUE))

# Check
file <- "test_temp.hdf5"
dataset.d <- "/SVD/matrix/d"
dataset.u <- "/SVD/matrix/u"

res.d <-  h5read(file, dataset.d)
res.u <-  h5read(file, dataset.u)

res.d; resr$d

all.equal( round(as.vector(res.d), 5)[1:length(res.d)], round(resr$d, 5)[1:length(resr$d)])
all.equal( round(abs(res.u), 5)[1:nrow(res.u),1:ncol(res.u)], round(abs(resr$u), 5)[1:nrow(resr$u),1:ncol(resr$u)])



# Matrix - InvCholesky 
# -----------------------------------

library(BigDataStatMeth)
library(rhdf5)

# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")
set.seed(1234)
Y <- matrix(sample.int(10, 100, replace = TRUE), ncol = 10)
Ycp <- crossprod(Y)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Ycp, group = "data", dataset = "matrix",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

# Get Inverse Cholesky
res <- bdInvCholesky_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "invmatrix", outgroup = "InvCholesky", fullMatrix = FALSE, overwrite = TRUE)
res <- bdInvCholesky_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "invmatrix", outgroup = "InvCholesky", fullMatrix = TRUE, overwrite = TRUE)
res <- bdInvCholesky_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "invmatrix", outgroup = "InvCholesky", fullMatrix = TRUE, overwrite = TRUE, elementsBlock = 50)

# Check
file <- "test_temp.hdf5"
dataset <- "/InvCholesky/invmatrix"
res <-  h5read(file, dataset)
resr <- solve(Ycp)
all.equal(resr[upper.tri(resr)], res[upper.tri(res)])

microbenchmark::microbenchmark( T <- solve(Ycp),
                                res <- res <- bdInvCholesky_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "invmatrix", outgroup = "InvCholesky", fullMatrix = FALSE, overwrite = TRUE),
                                times = 5 )


# Matrix - Cholesky Decomposition
# -----------------------------------

setwd("/Users/mailos/PhD/dummy")
set.seed(1234)
Y <- matrix(sample.int(10, 10000, replace = TRUE), ncol = 100)
Ycp <- crossprod(Y)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Ycp, group = "data", dataset = "matrix",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCholesky_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrixDec", outgroup = "Cholesky_Dec", fullMatrix = FALSE, overwrite = TRUE)
bdCholesky_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrixDec", outgroup = "Cholesky_Dec", fullMatrix = FALSE, overwrite = TRUE, elementsBlock = 50)

res <-  h5read("test_temp.hdf5", "Cholesky_Dec/matrixDec")
all.equal(res, chol(Ycp))




# Matrix - QR Decomposition
# -----------------------------------


# On-memory execution

res <- bdQR(Y, thin = FALSE)

all.equal(res$Q, qr.Q(qr(Y)))
all.equal(res$R[1:10,1:10], qr.R(qr(Y))[1:10,1:10])


# HDF5 data file execution
# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdQR_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "QR_Dec",  overwrite = TRUE)
bdQR_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "QR_Dec",  overwrite = TRUE, block_size = 256)
bdQR_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "QR_Dec", kelements = 4 , overwrite = TRUE)


resR <-  h5read("test_temp.hdf5", "QR_Dec/R.matrix")
resQ <-  h5read("test_temp.hdf5", "QR_Dec/Q.matrix")

all.equal(resR, qr.R(qr(Y))[1:10,1:10])
all.equal(resQ, qr.Q(qr(Y)))


microbenchmark::microbenchmark( T <- qr(Y),
                                res = bdQR_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "QR_Dec", thin = FALSE , overwrite = TRUE),
                                resp = bdQR_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "QR_Dec", thin = FALSE , overwrite = TRUE, block_size = 1024),
                                rest = bdQR_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "QR_Dec", thin = TRUE , overwrite = TRUE),
                                respt = bdQR_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "QR_Dec", thin = TRUE , overwrite = TRUE, block_size = 1024),
                                times = 3 )



# Matrix - pseudoinverse
# -----------------------------------

# On-memory execution
resm <- bdpseudoinv(Y)
all.equal(resm, pseudoinverse(Y))


# HDF5 data file execution
# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdpseudoinv_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "PseudoInv",  overwrite = TRUE)

resf <-  h5read("test_temp.hdf5", "PseudoInv/matrix")
all.equal(resf, pseudoinverse(Y))


microbenchmark::microbenchmark( R <- pseudoinverse(Y),
                                resf = bdpseudoinv_hdf5(filename = "test_temp.hdf5", group = "data", dataset = "matrix", outdataset = "matrix", outgroup = "PseudoInv",  overwrite = TRUE),
                                resm = bdpseudoinv(Y),
                                times = 5 )




# Matrix - solve equation
# -----------------------------------

library("BigDataStatMeth")
library(rhdf5)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
setwd("/Users/mailos/PhD/dummy")

N = 1800
M = 1800

set.seed(555)
Y <- matrix(rnorm(N*M), N, M)
X <- matrix(rnorm(N), N, 1)
Ycp <- crossprod(Y)

# On-memory execution
resm <- bdSolve(Ycp, X)
resr <- solve(Ycp, X)

all.equal( resm, resr)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Ycp, group = "data", dataset = "A",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = X,  group = "data",  dataset = "B",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdSolve_hdf5( filename = "test_temp.hdf5", groupA = "data", datasetA = "A", groupB = "data", datasetB = "B", outgroup = "Solved", outdataset = "A_B", overwrite = TRUE )
resf <-  h5read("test_temp.hdf5", "Solved/A_B")

all.equal( as.numeric(resf), as.numeric(resr))


microbenchmark::microbenchmark( R <- solve(Ycp, X),
                                resf = bdSolve_hdf5( filename = "test_temp.hdf5", groupA = "data", datasetA = "A", groupB = "data", datasetB = "B", outgroup = "Solved", outdataset = "A_B", overwrite = TRUE ),
                                resm = bdSolve(Ycp, X),
                                times = 5 )




# Principal Component Analysis - PCA
# -----------------------------------

library("BigDataStatMeth")
library(rhdf5)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
setwd("/Users/mailos/PhD/dummy")

N = 200
M = 80

set.seed(555)
Y <- matrix(rnorm(N*M), N, M)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Y, group = "data", dataset = "Y",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdPCA_hdf5("test_temp.hdf5", "data", "Y", ncomponents = 0L, bcenter = TRUE, bscale = TRUE, k = 2L, q = 1L, rankthreshold = 0.0, SVDgroup = NULL, overwrite = TRUE, threads = NULL) 

res_lambda <-  h5read("test_temp.hdf5", "/PCA/Y/lambda")
res_var <-  h5read("test_temp.hdf5", "PCA/Y/variance")
res_cumvar <-  h5read("test_temp.hdf5", "PCA/Y/cumvar")
res_var.coord <-  h5read("test_temp.hdf5", "PCA/Y/var.coord")
res_var.cos2 <-  h5read("test_temp.hdf5", "PCA/Y/var.cos2")
res_components <-  h5read("test_temp.hdf5", "PCA/Y/components")




library(stats)
resr <- prcomp(Y, center = TRUE, scale. = TRUE)

library(psych)
resr <- principal(scale(Y, center = TRUE, scale = TRUE))
resr$values
psych::print.psych(resr, cut = 1, sort=TRUE)
#















##########
##########      CHECKING IN PROGRESS..... 
##########


## FER UNA FUNCIÓ QUE PERMETI LLEGIR LES DADES DELS FITXERS HDF5 TAL QUAL
##      SENSE FER LA MODIFICACIÓ DE FILES - COLUMNES PER ADAPTAR-HO A LES DADES DE R .... 
##      REVISAR EL PROCÉS PERQUÈ PRIMER HAURIA DE COMPROVAR QUE REALMENT ESTIC FENT EL CAPGIRAMENT...
##      HAURIA D'IMPLEMENTAR LA FUNCIÓ QUE TENIA IMPLEMENTADA COM: 
##      GetCurrentBlock_hdf5_Original() A L'ANTIGA VERSIÓ DE BIGDATASTATMETH
##
# devtools::reload(pkgload::inst("BigDataStatMeth"))



file <- "test_temp.hdf5"
dataset.d <- "/SVD/matrix/d"
dataset.u <- "/SVD/matrix/u"

res.d <-  h5read(file, dataset.d)
res.u <-  h5read(file, dataset.u)

res.d; resr$d
res.u[1:10,1:10]; resr$u[1:10,1:10]

all.equal( round(as.vector(res.d), 5)[1:13], round(resr$d, 5)[1:13])
all.equal( round(abs(res.u), 5)[1:10,1:10], round(abs(resr$u), 5)[1:10,1:10])





##########      END CHECKING IN PROGRESS....

