# Matrix Crossprodduct
# ---------------------



# Test blockmult

library("BigDataStatMeth")
library("rhdf5")
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

N <- 2171
K <- 300
M <- 2171
L <- 739

# N <- 21
# K <- 9
# M <- 21
# L <- 6


set.seed(555)
a <- matrix( rnorm( N*K, mean=0, sd=1), N, K) 
b <- matrix( rnorm( L*M, mean=0, sd=1), M, L) 


# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(a), group = "pepet", 
                     dataset = "ta",
                     transp = FALSE,
                     overwriteFile = TRUE, 
                     overwriteDataset = FALSE, 
                     unlimited = FALSE)


bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = t(b), 
                     group = "pepet", 
                     dataset = "tb",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)





# .---------------------

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdtCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "ta", B = "tb", overwrite = TRUE , paral = TRUE) # tCrossprod A %*% t(B) 

file <- "test_temp.hdf5"
dataset <- "OUTPUT/tCrossProd_ta_x_tb"
resr <- tcrossprod(t(a),t(b))
res <-  h5read(file,dataset)
res[1:5,1:5]
resr[1:5,1:5]
all.equal( res, resr)

all.equal( res[1:5,1:5], resr[1:5,1:5])


# .---------------------

filename <- "/Users/mailos/PhD/dummy/Analyses/TCGA_CCA/cca_tcga_small.hdf5"

res <- bdCrossprod_hdf5(filename = filename,
                        group = "Step6", A = "XQ",
                        groupB = "Step6", B = "YQ",
                        outgroup = "Step7")

# -----     TEST
#
h5f <-  H5Fopen(filename)
XQ <- h5f$Step6$XQ
YQ <- h5f$Step6$YQ

res <- h5f$Step7$CrossProd_XQ_x_YQ
h5closeAll()

resr <- crossprod(XQ,YQ)
# -----



