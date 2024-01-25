
# Test Apply functions

library("BigDataStatMeth")
library(rhdf5)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

# setwd("~/PhD/TREBALLANT/BigDataStatMethAPIsV3")
setwd("/Users/mailos/PhD/dummy")

# Prepare data and functions
set.seed(123)
Y <- matrix(rnorm(100), 10, 10)
X <- matrix(rnorm(400), 20, 20)
Z <- matrix(rnorm(300), 20, 15)
filename <- "test_temp.hdf5"

# devtools::reload(pkgload::inst("BigDataStatMeth"))

# Create hdf5 datasets
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Y, group = "data", dataset = "Y",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = X,  group = "data",  dataset = "X",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5",
                     object = Z,  group = "data",  dataset = "Z",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE,
                     unlimited = FALSE)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
dsets <- bdgetDatasetsList_hdf5("test_temp.hdf5", group = "data")
dsets


# -----------------------------------
# Apply function :  QR Decomposition
# -----------------------------------

bdapply_Function_hdf5(filename = filename,
                      group = "data",datasets = dsets,
                      outgroup = "QR",func = "QR",
                      force = TRUE)

QX <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/QR/X.Q")
RX <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/QR/X.R")

QY <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/QR/Y.Q")
RY <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/QR/Y.R")

QZ <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/QR/Z.Q")
RZ <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/QR/Z.R")

all.equal(QY, qr.Q(qr(Y)))
all.equal(RY, qr.R(qr(Y)))

all.equal(QX, qr.Q(qr(X)))
all.equal(RX, qr.R(qr(X)))

all.equal(QZ, qr.Q(qr(Z)))
all.equal(RZ, qr.R(qr(Z)))


# ---------------------------
# Apply function :  CrossProd
# ---------------------------

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdapply_Function_hdf5(filename = filename,
                      group = "data",datasets = dsets,
                      outgroup = "Crossprod",func = "CrossProd",
                      force = TRUE)

CpX <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/Crossprod/X")
CpY <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/Crossprod/Y")
CpZ <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/Crossprod/Z")

all.equal(CpX, crossprod(X))
all.equal(CpY, crossprod(Y))
all.equal(CpZ, crossprod(Z))


# ---------------------------
# Apply function : tCrossProd
# ---------------------------

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdapply_Function_hdf5(filename = filename,
                      group = "data",datasets = dsets,
                      outgroup = "tCrossprod",func = "tCrossProd",
                      force = TRUE)

tCpX <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossprod/X")
tCpY <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossprod/Y")
tCpZ <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossprod/Z")

all.equal(tCpX, tcrossprod(X))
all.equal(tCpY, tcrossprod(Y))
all.equal(tCpZ, tcrossprod(Z))


# ----------------------------------
# Apply function : Inverse Cholesky
# ----------------------------------

# En aquest exemple no es poden carregar les matrius: 
#       - X i Y: no definides positives
#       - Z: no es quadrada
# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdapply_Function_hdf5(filename = filename,
                      group = "data",datasets = dsets,
                      outgroup = "invChol",func = "invChol",
                      force = TRUE)


# Redefinim matrius per a que X i Y siguin quadrades positives
#       - X i Y: no definides positives
#       - Z: no es quadrada
# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = tcrossprod(X),  group = "PositiveDefinite",  dataset = "XP",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = crossprod(Y), group = "PositiveDefinite", dataset = "YP",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)


dsetsPositive <- bdgetDatasetsList_hdf5("test_temp.hdf5", group = "PositiveDefinite")
dsetsPositive

bdapply_Function_hdf5(filename = filename,
                      group = "PositiveDefinite",datasets = dsetsPositive,
                      outgroup = "invChol", func = "invChol",
                      force = TRUE)

invXP <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/invChol/XP")
invYP <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/invChol/YP")

resr <- solve(tcrossprod(X))
all.equal(resr[upper.tri(resr)], invXP[upper.tri(invXP)])
resr <- solve(crossprod(Y))
all.equal(upper.tri(resr), upper.tri(invYP))


# ---------------------------------------
# Apply function : Cholesky Decomposition
# ---------------------------------------

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdapply_Function_hdf5(filename = filename,
                      group = "PositiveDefinite",datasets = dsetsPositive,
                      outgroup = "descChol", func = "descChol",
                      force = TRUE)

descXP <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/descChol/XP")
descYP <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/descChol/YP")

resr <- chol(tcrossprod(X))
all.equal(resr, descXP)
resr <- chol(crossprod(Y))
all.equal(resr, descYP)


# ---------------------------------------------------
# Apply function : Matrix multiplication (2 matrices)
# ---------------------------------------------------

# devtools::reload(pkgload::inst("BigDataStatMeth"))
N = 50
M = 30

set.seed(555)
a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
b <- matrix( rnorm( N*M, mean=0, sd=1), M, N) 

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = a,  group = "groupA",  dataset = "aa",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = b, group = "groupA", dataset = "ab",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = a,  group = "groupB",  dataset = "ba",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = b, group = "groupB", dataset = "bb",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)


dsets <- bdgetDatasetsList_hdf5("test_temp.hdf5", group = "groupA")
dsetsb <- bdgetDatasetsList_hdf5("test_temp.hdf5", group = "groupB")



# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = dsets, outgroup = "blockmult", func = "blockmult",  b_group = "groupB",  b_datasets = dsetsb, force = TRUE) 

#       No results - no es poden computar

bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = dsets, outgroup = "blockmult", func = "blockmult",  b_group = "groupB",  transp_dataset = FALSE, transp_bdataset = TRUE, b_datasets = dsetsb, force = TRUE) 

res1 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/blockmult/aa_ba")
res2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/blockmult/ab_bb")
all.equal(a%*%t(a), res1)
all.equal(b%*%t(b), res2)

bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = c("aa"), outgroup = "blockmult", func = "blockmult",  b_group = "groupB",  transp_dataset = FALSE, transp_bdataset = FALSE, b_datasets = c("bb"), force = TRUE) 

res <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/blockmult/aa_bb")
all.equal(a%*%b, res)

bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = c("aa"), outgroup = "blockmult", func = "blockmult",  b_group = "groupB",  transp_dataset = TRUE, transp_bdataset = FALSE, b_datasets = c("ba"), force = TRUE) 

res <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/blockmult/aa_ba")
all.equal( t(a)%*%a, res)


# ---------------------------------------------------
# Apply function : CrossProd (2 matrices)
# ---------------------------------------------------


# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = dsets, outgroup = "CrossProdD", func = "CrossProd",  b_group = "groupB",  b_datasets = dsetsb, force = TRUE) 

res1 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/CrossProdD/Cross_aa_ba")
res2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/CrossProdD/Cross_ab_bb")
all.equal(crossprod(a,a), res1)
all.equal(crossprod(b,b), res2)



#
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(a),  group = "groupA",  dataset = "aa",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(b), group = "groupA", dataset = "ab",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)
bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = dsets, outgroup = "CrossProdD", func = "CrossProd",  b_group = "groupB",  b_datasets = dsetsb, transp_dataset = TRUE, force = TRUE) 

res1 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/CrossProdD/Cross_aa_ba")
res2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/CrossProdD/Cross_ab_bb")
all.equal(crossprod(a,a), res1)
all.equal(crossprod(b,b), res2)


#
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(a),  group = "groupB",  dataset = "ba",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(b), group = "groupB", dataset = "bb",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)
bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = dsets, outgroup = "CrossProdD", func = "CrossProd",  b_group = "groupB",  b_datasets = dsetsb, transp_dataset = TRUE, transp_bdataset = TRUE, force = TRUE) 

res1 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/CrossProdD/Cross_aa_ba")
res2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/CrossProdD/Cross_ab_bb")
all.equal(crossprod(a,a), res1)
all.equal(crossprod(b,b), res2)



# ---------------------------------------------------
# Apply function : tCrossProd (2 matrices)
# ---------------------------------------------------

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = dsets, outgroup = "tCrossProdD", func = "tCrossProd",  b_group = "groupB",  b_datasets = dsetsb, force = TRUE) 

res1 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossProdD/tCross_aa_ba")
res2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossProdD/tCross_ab_bb")
all.equal(tcrossprod(a,a), res1)
all.equal(tcrossprod(b,b), res2)


#
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(a),  group = "groupA",  dataset = "aa",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(b), group = "groupA", dataset = "ab",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)
bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = dsets, outgroup = "tCrossProdD", func = "tCrossProd",  b_group = "groupB",  b_datasets = dsetsb, transp_dataset = TRUE, force = TRUE) 

res1 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossProdD/tCross_aa_ba")
res2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossProdD/tCross_ab_bb")
all.equal(tcrossprod(a,a), res1)
all.equal(tcrossprod(b,b), res2)

#
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(a),  group = "groupB",  dataset = "ba",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(b), group = "groupB", dataset = "bb",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)
bdapply_Function_hdf5(filename = filename, group = "groupA", datasets = dsets, outgroup = "tCrossProdD", func = "tCrossProd",  b_group = "groupB",  b_datasets = dsetsb, transp_dataset = TRUE, transp_bdataset = TRUE, force = TRUE) 

res1 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossProdD/tCross_aa_ba")
res2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/tCrossProdD/tCross_ab_bb")
all.equal(tcrossprod(a,a), res1)
all.equal(tcrossprod(b,b), res2)



# ---------------------------------------------------
# Apply function : Solve Matrix equation
# ---------------------------------------------------


# devtools::reload(pkgload::inst("BigDataStatMeth"))

N = 1800
M = 1800

set.seed(555)
Y <- matrix(rnorm(N*M), N, M)
X <- matrix(rnorm(N), N, 1)
Ycp <- crossprod(Y)


# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Ycp, group = "matrices", dataset = "A",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = X,  group = "vectors",  dataset = "vA",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

resr1 <- solve(Ycp, X)

# Add another group Matrx - Vector

set.seed(42314)
Y <- matrix(rnorm(N*M), N, M)
X <- matrix(rnorm(N), N, 1)
Ycp <- crossprod(Y)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Ycp, group = "matrices", dataset = "B",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = X,  group = "vectors",  dataset = "vB",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)
resr2 <- solve(Ycp, X)

MatrixDatasets <- bdgetDatasetsList_hdf5("test_temp.hdf5", group = "matrices")
VecrorsDatasets <- bdgetDatasetsList_hdf5("test_temp.hdf5", group = "vectors")



# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdapply_Function_hdf5(filename = filename, group = "matrices", datasets = MatrixDatasets, outgroup = "Solved", func = "solve",  b_group = "vectors",  b_datasets = VecrorsDatasets, force = TRUE) 

res1 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/Solved/A_eq_vA")
res2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/Solved/B_eq_vB")

all.equal( resr1, res1)
all.equal( resr2, res2)



# ---------------------------------------------------
# Apply function : get sd and mean
# ---------------------------------------------------


# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

# Prepare data and functions
set.seed(123)
Y <- matrix(rnorm(100), 10, 10)
X <- matrix(rnorm(400), 20, 20)
Z <- matrix(rnorm(300), 20, 15)
filename <- "test_temp.hdf5"

# devtools::reload(pkgload::inst("BigDataStatMeth"))

# Create hdf5 datasets
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Y, group = "data", dataset = "Y",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = X,  group = "data",  dataset = "X",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5",
                     object = Z,  group = "data",  dataset = "Z",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE,
                     unlimited = FALSE)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
dsets <- bdgetDatasetsList_hdf5("test_temp.hdf5", group = "data")
dsets


bdapply_Function_hdf5(filename = filename, group = "data", datasets = dsets, outgroup = "sdmean", func = "sdmean", byrows = FALSE,  force = TRUE) 

meanX <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/mean.X")
sdX <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/sd.X")
meanY <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/mean.Y")
sdY <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/sd.Y")
meanZ <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/mean.Z")
sdZ <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/sd.Z")

all.equal(as.numeric(apply(X,2,mean)), as.numeric(meanX)); all.equal(as.numeric(apply(X,2,sd)), as.numeric(sdX))
all.equal(as.numeric(apply(Y,2,mean)), as.numeric(meanY)); all.equal(as.numeric(apply(Y,2,sd)), as.numeric(sdY))
all.equal(as.numeric(apply(Z,2,mean)), as.numeric(meanZ)); all.equal(as.numeric(apply(Z,2,sd)), as.numeric(sdZ))


bdapply_Function_hdf5(filename = filename, group = "data", datasets = dsets, outgroup = "sdmean", func = "sdmean", byrows = TRUE,  force = TRUE) 

rmeanX <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/mean.X")
rsdX <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/sd.X")
rmeanY <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/mean.Y")
rsdY <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/sd.Y")
rmeanZ <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/mean.Z")
rsdZ <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/sdmean/sd.Z")

all.equal(as.numeric(apply(X,1,mean)), as.numeric(rmeanX)); all.equal(as.numeric(apply(X,1,sd)), as.numeric(rsdX))
all.equal(as.numeric(apply(Y,1,mean)), as.numeric(rmeanY)); all.equal(as.numeric(apply(Y,1,sd)), as.numeric(rsdY))
all.equal(as.numeric(apply(Z,1,mean)), as.numeric(rmeanZ)); all.equal(as.numeric(apply(Z,1,sd)), as.numeric(rsdZ))



