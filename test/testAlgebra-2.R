
# a <- matrix(seq(1:25), nrow = 5)
# 
# f <- matrix(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), nrow = 1)
# f <- matrix(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), nrow = 3)
# 
# b <-  as.vector(c(1,2,3,4,5))
# c <-  c(1,2,3,4,5)
# 
# 
# library("BigDataStatMeth")
# # devtools::reload(pkgload::inst("BigDataStatMeth"))
# 
# setwd("/Users/mailos/PhD/dummy")
# bdCreate_hdf5_matrix_file("test_temp.hdf5", f, 
#                           "pepet", "datasetpepet",
#                           FALSE, TRUE)



# Test blockmult

library("BigDataStatMeth")
library(rhdf5)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

N = 20
M = 10

set.seed(555)
Y <- matrix(rnorm(110), 11, 10)
X <- matrix(rnorm(10), 10, 1)
# X <- c(1,2,3,4,5,6,7,8,9,10)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = Y, group = "data", dataset = "matrix",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = FALSE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = X,  group = "data",  dataset = "vector",
                     transp = FALSE,
                     overwriteFile = FALSE, overwriteDataset = TRUE, 
                     unlimited = FALSE)


# Matrix Vector - Calculus
# -----------------------------------

# - by Columns
bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "*", byrows = FALSE,  paral=FALSE, threads = 1, force = TRUE)
bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "-", byrows = FALSE,  paral=FALSE, threads = 1, force = TRUE)
bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "+", byrows = FALSE,  paral=FALSE, threads = 1, force = TRUE)
bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "/", byrows = FALSE,  paral=FALSE, threads = 1, force = TRUE)

# - by Rows
bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "*", byrows = TRUE,  paral=FALSE, threads = 1, force = TRUE)
bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "-", byrows = TRUE,  paral=FALSE, threads = 1, force = TRUE)
bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "+", byrows = TRUE,  paral=FALSE, threads = 1, force = TRUE)
bdcomputeMatrixVector_hdf5("test_temp.hdf5", group = "data", dataset = "matrix", vectorgroup = "data", vectordataset = "vector",  
                           outgroup = "results", outdataset = "res",  func = "/", byrows = TRUE,  paral=FALSE, threads = 1, force = TRUE)

# Matrix - Sd and mean by rows or Cols
# -----------------------------------

# - by Rows
bdgetSDandMean_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", sd = TRUE, mean = TRUE, byrows = TRUE, force = TRUE)
bdgetSDandMean_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", sd = TRUE, mean = TRUE, byrows = TRUE, wsize = 3, force = TRUE)

# - by Columns
bdgetSDandMean_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", sd = TRUE, mean = TRUE, byrows = FALSE, force = TRUE)
bdgetSDandMean_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", sd = TRUE, mean = TRUE, byrows = FALSE, wsize = 3, force = TRUE)

# # -- Test especial Sd and mean
#                 file <- "test_temp.hdf5"
#                 dataset.mean <- "mean_sd/mean.matrix"
#                 dataset.sd <- "mean_sd/sd.matrix"
# 
#                 resr.mean <- apply(Y, 2, mean)
#                 resr.sd <- apply(Y, 2, sd)
#                 res.mean <-  h5read(file,dataset.mean)
#                 res.sd <-  h5read(file,dataset.sd)
#                 res.mean; resr.mean
#                 res.sd; resr.sd
#                 all.equal( round(as.vector(res.mean), 5), round(resr.mean, 5))
#                 all.equal( round(as.vector(res.sd), 5), round(resr.sd, 5))
                          

# Matrix - Normalization
# -----------------------------------

# - by Colums
# - - Only center
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = TRUE, bscale = FALSE, byrows = FALSE, force  = TRUE) ; resr <- scale(Y, center = TRUE, scale = FALSE ) ; all.equal(resr, h5read("test_temp.hdf5", "NORMALIZED/data/matrix"), check.attributes = FALSE ) 
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = TRUE, bscale = FALSE, byrows = FALSE, wsize = 3, force  = TRUE); resr <- scale(Y, center = TRUE, scale = FALSE ) ; all.equal(resr, h5read("test_temp.hdf5", "NORMALIZED/data/matrix"), check.attributes = FALSE ) 
# - - Only scale
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = FALSE, bscale = TRUE, byrows = FALSE, force  = TRUE) 
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = FALSE, bscale = TRUE, byrows = FALSE, wsize = 3, force  = TRUE)
# - - Center and Scale
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = TRUE, bscale = TRUE, byrows = FALSE, force  = TRUE); resr <- scale(Y, center = TRUE, scale = TRUE ) ; all.equal(resr, h5read("test_temp.hdf5", "NORMALIZED/data/matrix"), check.attributes = FALSE ) 
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = TRUE, bscale = TRUE, byrows = FALSE, wsize = 3, force  = TRUE); resr <- scale(Y, center = TRUE, scale = TRUE ) ; all.equal(resr, h5read("test_temp.hdf5", "NORMALIZED/data/matrix"), check.attributes = FALSE ) 

# - by Rows
# - - Only center
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = TRUE, bscale = FALSE, byrows = TRUE, force  = TRUE); resr <- scale(t(Y), center = TRUE, scale = FALSE ) ; all.equal(t(resr), h5read("test_temp.hdf5", "NORMALIZED/data/matrix"), check.attributes = FALSE ) 
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = TRUE, bscale = FALSE, byrows = TRUE, wsize = 3, force  = TRUE); resr <- scale(t(Y), center = TRUE, scale = FALSE ) ; all.equal(t(resr), h5read("test_temp.hdf5", "NORMALIZED/data/matrix"), check.attributes = FALSE ) 
# - - Only scale
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = FALSE, bscale = TRUE, byrows = TRUE, force  = TRUE)
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = FALSE, bscale = TRUE, byrows = TRUE, wsize = 3, force  = TRUE)
# - - Center and Scale
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = TRUE, bscale = TRUE, byrows = TRUE, force  = TRUE); resr <- scale(t(Y), center = TRUE, scale = TRUE ) ; all.equal(t(resr), h5read("test_temp.hdf5", "NORMALIZED/data/matrix"), check.attributes = FALSE ) 
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = TRUE, bscale = TRUE, byrows = TRUE, wsize = 3, force  = TRUE); resr <- scale(t(Y), center = TRUE, scale = TRUE ) ; all.equal(t(resr), h5read("test_temp.hdf5", "NORMALIZED/data/matrix"), check.attributes = FALSE ) 

# # -- Test especial Normalization - (only scale)
        # Only Scale:
        #..# resr <- scale(Y, center = FALSE, scale = apply(Y, 2, sd, na.rm = TRUE))





##########
##########      CHECKING IN PROGRESS..... 
##########

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = FALSE, bscale = TRUE, byrows = TRUE, force  = TRUE)
bdNormalize_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", bcenter = FALSE, bscale = TRUE, byrows = TRUE, wsize = 3, force  = TRUE)


file <- "test_temp.hdf5"
dataset <- "NORMALIZED/data/matrix"

# Only Scale:
#..# resr <- scale(t(Y), center = FALSE, scale = apply(t(Y), 2, sd, na.rm = TRUE))
# Center or Center + Scale
resr <- scale(t(Y), TRUE, FALSE)

res <-  h5read(file, dataset)
res ; resr
all.equal( round(res, 5), t(round(resr, 5)), check.attributes = FALSE)




## -----------------------------------------------------------------------
##  !!! --ATENCIÓ, A PARTIR D'AQUÍ COMENCEN A SER PARAULES MAJORS !!!!
## -----------------------------------------------------------------------



# Matrix - InvCholesky 
# -----------------------------------


# devtools::reload(pkgload::inst("BigDataStatMeth"))

# -->
# --> Estic aquí, programant i testejant .... a decidir quina funció programo 
# -->










##########
##########      CHECKING IN PROGRESS..... 
##########

devtools::reload(pkgload::inst("BigDataStatMeth"))

bdgetSDandMean_hdf5( "test_temp.hdf5", group = "data", dataset = "matrix", sd = TRUE, mean = TRUE, byrows = FALSE, wsize = 3, force = TRUE)


file <- "test_temp.hdf5"
dataset.mean <- "mean_sd/mean.matrix"
dataset.sd <- "mean_sd/sd.matrix"

resr.mean <- apply(Y, 2, mean)
resr.sd <- apply(Y, 2, sd)
res.mean <-  h5read(file,dataset.mean)
res.sd <-  h5read(file,dataset.sd)
res.mean; resr.mean
res.sd; resr.sd
all.equal( round(as.vector(res.mean), 5), round(resr.mean, 5))
all.equal( round(as.vector(res.sd), 5), round(resr.sd, 5))




##########      END CHECKING IN PROGRESS....



