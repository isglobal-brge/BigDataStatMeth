## ----cleanup, echo=FALSE, include=FALSE---------------------------------------
if(file.exists('delayed.hdf5'))
    file.remove('delayed.hdf5')
if(file.exists('robjects.hdf5'))
    file.remove('robjects.hdf5')

## ----setup, include = FALSE---------------------------------------------------
library(BiocStyle)
knitr::opts_chunk$set(collapse = TRUE, comment = "", cache=TRUE, message = FALSE, width = 180)

## ----r_setup, include = FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 180)

## ----load-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(BigDataStatMeth)
library(DelayedArray)
library(rhdf5)
library(DT)

## ----load2------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(microbenchmark)

## ----hdf5Img, fig.align = 'center', fig.cap = "HDF5 hierarchical structure", echo=FALSE-------------------------------------------------------------------------------------------
knitr::include_graphics("imgs/hdf5_squema.jpg")

## ----hdf5Create-------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(rhdf5)
library(DelayedArray)

set.seed(5234)
n <- 100
m <- 10000
A <- matrix(rnorm(n*m,mean=0,sd=1), n,m)

Ad <- DelayedArray(A)

# Create a file with a dataset from DelayedMatrix Ad in INPUT group
Create_HDF5_matrix_file("delayed.hdf5", Ad, "INPUT", "A")

# We also can create a dataset from R matrix object
Create_HDF5_matrix_file("robjects.hdf5", A, "INPUT", "A")

## ----hdf5AddDataset---------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(5234)
n <- 50
m <- 12000
B <- matrix(rnorm(n*m,mean=3,sd=0.5), n,m)

Bd <- DelayedArray(B)

# Create dataframe B with Bd matrix in delayed.hdf5 file at INPUT group
Create_HDF5_matrix(Bd, "delayed.hdf5", "INPUT", "B");

# Create dataframe data with Ad matrix in delayed.hdf5 file at OMIC group
set.seed(5234)
n <- 150000
m <- 50
odata <- matrix(rnorm(n*m,mean=0,sd=1), n,m)

Create_HDF5_matrix(odata, "delayed.hdf5", "OMICS", "data")


## ----hdf5Open---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Examine hierarchy before open file
h5ls("delayed.hdf5")

# Open file
h5fdelay = H5Fopen("delayed.hdf5")
# Show hdf5 hierarchy (groups)
h5fdelay


## ----hdf5Dataset------------------------------------------------------------------------------------------------------------------------------------------------------------------
Adata = h5fdelay$OMICS$data
Adata[1:3,1:6]

## ----hdf5Close--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Close delayed.hdf5 file
H5Fclose(h5fdelay)

# Open 2 files and close all
h5fdelay = H5Fopen("delayed.hdf5")
h5fr = H5Fopen("robjects.hdf5")

h5closeAll()

## ----mat_sim----------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123456)
n <- 500
p <- 750
A <- matrix(rnorm(n*n), nrow=n, ncol=n)
B <- matrix(rnorm(n*p), nrow=n, ncol=p)

n <- 1000
p <- 10000
Abig <- matrix(rnorm(n*n), nrow=n, ncol=n)
Bbig <- matrix(rnorm(n*p), nrow=n, ncol=p)

## ----getDelayed-------------------------------------------------------------------------------------------------------------------------------------------------------------------
DA <- DelayedArray(A)
DB <- DelayedArray(B)

## ----mat_mult---------------------------------------------------------------------------------------------------------------------------------------------------------------------
AxB <- blockmult(A, B)
AxBDelay <- blockmult(DA, DB)
AxBDelay [1:5,1:5]

## ----check------------------------------------------------------------------------------------------------------------------------------------------------------------------------
all.equal(AxB, AxBDelay)
all.equal(A%*%B, AxBDelay)

## ----noblockmultparal-------------------------------------------------------------------------------------------------------------------------------------------------------------
AxB <- blockmult(A, B, paral = TRUE)
AxBDelay <- blockmult(DA, DB, paral = TRUE) 

all.equal(AxB, AxBDelay)

## ----benchmark1-------------------------------------------------------------------------------------------------------------------------------------------------------------------
microbenchmark(
  R = A%*%B,
  noparal = blockmult(A, B),
  paral = blockmult(A, B, paral=TRUE),
  times = 10)

## ----noblockmultparalthreads------------------------------------------------------------------------------------------------------------------------------------------------------
AxB <- blockmult(A, B, paral = TRUE, threads = 2)
AxBDelay <- blockmult(DA, DB, paral = TRUE, threads = 3) 

all.equal(AxB, AxBDelay)

## ----benchmark1b------------------------------------------------------------------------------------------------------------------------------------------------------------------
microbenchmark(
  AxBnThread2 = blockmult(DA, DB, paral = TRUE, threads = 2) ,
  AxBnThread3 = blockmult(DA, DB, paral = TRUE, threads = 3) ,
  AxBnThread4 = blockmult(DA, DB, paral = TRUE, threads = 4) ,
  times = 10)

## ----blockmult--------------------------------------------------------------------------------------------------------------------------------------------------------------------
AxB <- blockmult(A, B, block_size = 10)
AxBDelay <- blockmult(DA, DB, block_size = 10 )

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
all.equal(AxBDelay,A%*%B)
all.equal(AxB, AxBDelay)

## ----blockmultparal---------------------------------------------------------------------------------------------------------------------------------------------------------------
AxB <- blockmult(A, B, block_size = 10, paral = TRUE)
AxBDelay <- blockmult(DA, DB, block_size = 10, paral = TRUE )
AxBDelayT <- blockmult(DA, DB, block_size = 10, paral = TRUE, threads =  3 )

all.equal(AxBDelay,A%*%B)
all.equal(AxB, AxBDelay)
all.equal(AxBDelay, AxBDelayT)

## ----blockmultbm1-----------------------------------------------------------------------------------------------------------------------------------------------------------------
DAbig <- DelayedArray(Abig)
DBbig <- DelayedArray(Abig)

# We consider a big matrix if number of rows or columns are > 500
AxBBig3000 <- blockmult(DAbig, DBbig, bigmatrix = 500)

# We want to force it to run in memory
AxBNOBig <- blockmult(DAbig, DBbig, bigmatrix = 100000) 


## ----blockmultresmat--------------------------------------------------------------------------------------------------------------------------------------------------------------
class(AxBNOBig)
AxBNOBig[1:5,1:5]

## ----blockmultresfile-------------------------------------------------------------------------------------------------------------------------------------------------------------
class(AxBBig3000)

# Assign hdf5 result content to R variable reshdf5
reshdf5 <- AxBBig3000$res$OUTPUT$C
reshdf5[1:5,1:5]

all.equal(reshdf5, AxBNOBig)

## ----blockmultresfileclose--------------------------------------------------------------------------------------------------------------------------------------------------------
bdclose(AxBBig3000)

## ----benchmark2, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------
bench1 <- microbenchmark(
  noblockParal = blockmult(Abig, Bbig, paral = TRUE),
  blockParalfile = blockmult(Abig, Bbig, block_size = 256, paral=TRUE),
  blockParalmem = blockmult(Abig, Bbig, block_size = 256, paral=TRUE, bigmatrix = 100000),
  times = 3 ) # Tornar-ho a posar a 10 quan tot OK per estalviar temps d'execució !!!

bench1

## ---- bench1, fig.height=3, fig.width=3-------------------------------------------------------------------------------------------------------------------------------------------
ggplot2::autoplot(bench1)

## ----crossprod--------------------------------------------------------------------------------------------------------------------------------------------------------------------
n <- 500
m <- 250
A <- matrix(rnorm(n*m), nrow=n, ncol=m)
DA <- DelayedArray(A)

# Cross Product of a standard R matrix
cpA <- bdcrossprod(A)  
# Result with DelayedArray data type
cpDA <- bdcrossprod(DA)  

## ----check_cp---------------------------------------------------------------------------------------------------------------------------------------------------------------------
all.equal(cpDA, crossprod(A))

## ----nocrossprod------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Transposed Cross Product R matrices
tcpA <- bdcrossprod(A, transposed = TRUE)
# With DelayeArray data types
tcpDA <- bdcrossprod(DA, transposed = TRUE) 

## ----check_tcp--------------------------------------------------------------------------------------------------------------------------------------------------------------------
all.equal(tcpDA, tcrossprod(A))

## ----benchmark_bdcrossprod--------------------------------------------------------------------------------------------------------------------------------------------------------
bench <- microbenchmark(
  bdcrossp = bdcrossprod(DA, transposed = TRUE),
  rcrossp = tcrossprod(A),
  times = 1)  # Only DEBUG MODE 
  # times = 10) # PRODUCTION MODE
bench

## ----wcrossprod-------------------------------------------------------------------------------------------------------------------------------------------------------------------
n <- 250
X <- matrix(rnorm(n*n), nrow=n, ncol=n)
u <- runif(n)
w <- u * (1 - u)
DX <- DelayedArray(X)
Dw <- DelayedArray(as.matrix(w))
  
wcpX <- bdwproduct(X, w,"xwxt")
wcpDX <- bdwproduct(DX, Dw,"xwxt") # with DelayedArray

wcpDX[1:5,1:5]

## ----check_wcp--------------------------------------------------------------------------------------------------------------------------------------------------------------------
all.equal( wcpDX, X%*%diag(w)%*%t(X) )

## ----wtcrossprod------------------------------------------------------------------------------------------------------------------------------------------------------------------
wtcpX <- bdwproduct(X, w,"xtwx")
wtcpDX <- bdwproduct(DX, Dw,"xtwx") # with DelayedArray

wtcpDX[1:5,1:5]

## ----check_wtcp-------------------------------------------------------------------------------------------------------------------------------------------------------------------
all.equal(wtcpDX, t(X)%*%diag(w)%*%X)

## ----invChols---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate a positive definite matrix
Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

A <- Posdef(n = 500, ev = 1:500)
DA <- DelayedArray(A)

invchol <- bdInvCholesky(A)
Dinvchol <- bdInvCholesky(DA)

round(invchol[1:5,1:5],8)

## ----invCholsequal----------------------------------------------------------------------------------------------------------------------------------------------------------------
all.equal(Dinvchol, solve(A))

## ----svd_default------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Matrix simulation
set.seed(413)
n <- 500
A <- matrix(rnorm(n*n), nrow=n, ncol=n)
# Get a Delayed Array object
DA <- DelayedArray(A)

# SVD
bsvd <- bdSVD(A)
Dbsvd <- bdSVD(DA)

# Singular values, and right and left singular vectors
bsvd$d[1:5]
bsvd$u[1:5,1:5]
bsvd$v[1:5,1:5]

## ----check_svd--------------------------------------------------------------------------------------------------------------------------------------------------------------------
all.equal( sqrt( svd( tcrossprod( scale(A) ) )$d[1:10] ), bsvd$d[1:10] ) 
all.equal( sqrt( svd( tcrossprod( scale(A) ) )$d[1:10] ), Dbsvd$d[1:10] )

## ----svd_nonorm-------------------------------------------------------------------------------------------------------------------------------------------------------------------
bsvd <- bdSVD(A, bcenter = FALSE, bscale = FALSE)
Dbsvd <- bdSVD(DA, bcenter = FALSE, bscale = FALSE)


bsvd$d[1:5]
bsvd$u[1:5,1:5]
bsvd$v[1:5,1:5]

all.equal( sqrt(svd(tcrossprod(A))$d[1:10]), bsvd$d[1:10] ) 
all.equal( sqrt(svd(tcrossprod(A))$d[1:10]), Dbsvd$d[1:10] )
  

## ----BSVDImg, fig.align = 'center', fig.cap = "Flowchart for a two-level hierarchical Block SVD algorithm", echo=FALSE------------------------------------------------------------
knitr::include_graphics("imgs/blocksvd.png")

## ----BlockSVDNorm-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Direct from hdf5 data file
svdh5 <- bdSVD_hdf5("delayed.hdf5", "OMICS", "data")
svdh5$d[1:7]

# with R implementation from data in memory
fprova <- H5Fopen("delayed.hdf5")
omdata <- fprova$OMICS$data
h5closeAll()

svd <- svd(scale(omdata))
svd$d[1:7]

all.equal(svd$d,svdh5$d )

## ----BlockSVDNotNorm--------------------------------------------------------------------------------------------------------------------------------------------------------------
# Direct from hdf5 data file
svdh5 <- bdSVD_hdf5("delayed.hdf5", "OMICS", "data", bcenter = FALSE, bscale = FALSE)
svdh5$d[1:7]

# with R implementation from data in memory
svd <- svd(omdata)
svd$d[1:7]

all.equal(svd$d,svdh5$d )

## ----BlockSVDk4-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Block decomposition with 1 level and 4 local SVDs at each level
svdh5 <- bdSVD_hdf5("delayed.hdf5", "OMICS", "data", q=1, k=4 )
svdh5$d[1:7]

# with R implementation from data in memory
svd <- svd(scale(omdata))
svd$d[1:7]

all.equal(svd$d,svdh5$d )

## ----svdperform, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------
bdSVD_hdf5_Normalized_k2q1 <- "bdSVD_hdf5('delayed.hdf5',group='OMICS',dataset='data')"
bdSVD_hdf5_Normalized_k4q1 <- "bdSVD_hdf5('delayed.hdf5',group='OMICS',dataset='data',k=4)"
bdSVD_hdf5_Normalized_k4q2 <- "bdSVD_hdf5('delayed.hdf5',group='OMICS',dataset='data',q=2,k=4)"
bdSVD_memory_Normalized <- 'bdSVD_lapack(omdata)'
svd_R_memory_Normalized <- 'svd(scale(omdata))'

# With normalization process
res <- microbenchmark( eval(parse(text=bdSVD_hdf5_Normalized_k2q1)),
                       eval(parse(text=bdSVD_hdf5_Normalized_k4q1)),
                       eval(parse(text=bdSVD_hdf5_Normalized_k4q2)),
                       eval(parse(text=bdSVD_memory_Normalized)),
                       eval(parse(text=svd_R_memory_Normalized)),
                       times = 3, unit = "s")

print(summary(res)[, c(1:7)],digits=3)

## ----svdperform2, cache=TRUE------------------------------------------------------------------------------------------------------------------------------------------------------
bdSVD_hdf5_Not_Normalized_k2q1 <- "bdSVD_hdf5('delayed.hdf5',group='OMICS',dataset='data',bcenter=FALSE,bscale=FALSE)"
bdSVD_hdf5_Not_Normalized_k4q1 <- "bdSVD_hdf5('delayed.hdf5',group='OMICS',dataset='data',k=4,bcenter=FALSE,bscale=FALSE)"
bdSVD_hdf5_Not_Normalized_k4q2 <- "bdSVD_hdf5('delayed.hdf5',group='OMICS',dataset='data',q=2,k=4,bcenter=FALSE,bscale=FALSE)"
bdSVD_memory_Not_Normalized <- 'bdSVD_lapack(omdata, bcenter = TRUE, bscale = TRUE)'
svd_R_memory_Not_Normalized <- 'svd(omdata)'

# Without normalization process
res <- microbenchmark( eval(parse(text=bdSVD_hdf5_Not_Normalized_k2q1)),
                       eval(parse(text=bdSVD_hdf5_Not_Normalized_k4q1)),
                       eval(parse(text=bdSVD_hdf5_Not_Normalized_k4q2)),
                       eval(parse(text=bdSVD_memory_Not_Normalized)),
                       eval(parse(text=svd_R_memory_Not_Normalized)),
                       times = 3, unit = "s")

print(summary(res)[, c(1:7)],digits=3)

## ----mirNA------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data(miRNA)
dim(miRNA)

## ----tab--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
table(cancer)

## ----pca--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pc <- prcomp(miRNA)

## ----plot_pca---------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(pc$x[, 1], pc$x[, 2],  
     main = "miRNA data on tumor samples", 
     xlab = "PC1", ylab = "PC2", type="n")
abline(h=0, v=0, lty=2)
points(pc$x[, 1], pc$x[, 2], col = cancer, 
       pch=16, cex=1.2)
legend("topright", levels(cancer), pch=16, col=1:3)

## ----cia_da-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
miRNAD <- DelayedArray(miRNA)
miRNAD.c <- DelayedArray::sweep(miRNAD, 2,
      DelayedArray::colMeans(miRNAD), "-")
svd.da <- bdSVD(miRNAD.c, bcenter = FALSE, bscale = FALSE)

## ----plot_svd_da------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(svd.da$u[, 1], svd.da$u[, 2],  
     main = "miRNA data on tumor samples", 
     xlab = "PC1", ylab = "PC2", type="n")
abline(h=0, v=0, lty=2)
points(svd.da$u[, 1], svd.da$u[, 2], col = cancer, 
       pch=16, cex=1.2)
legend("topright", levels(cancer), pch=16, col=1:3)

## ----simul_lasso------------------------------------------------------------------------------------------------------------------------------------------------------------------
# number of samples
n <- 500
# number of variables
p <- 200
# covariates
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
# outcome (only variables 1, 2 and 5 are different from 0)
Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5] 

## ----looe-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get DealyedArray matrices
MD <- DelayedArray(M)
YD <- DelayedArray(as.matrix(Y))

# Model
mod <- LOOE(MD, YD, paral=FALSE)
mod$coef[abs(mod$coef)>mean(mod$coef)]

## ----glmnet-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(glmnet)
mod.cv <- cv.glmnet(M, Y)
mod.glmnet <- glmnet(M, Y, lambda=mod.cv$lambda.min)
mod.glmnet$beta[1:7,]

## ----impute-----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----sesinfo----------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

