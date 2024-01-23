
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

# Apply function :  QR Decomposition

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


# Apply function :  CrossProd

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


# Apply function : tCrossProd

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




