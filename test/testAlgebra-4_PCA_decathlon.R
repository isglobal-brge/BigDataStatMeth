library(BigDataStatMeth)
library(rhdf5)

library("FactoMineR")
library("factoextra")

data(decathlon2)

setwd("/Users/mailos/PhD/dummy")


decathlon2.active <- decathlon2[1:23, 1:10]

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = as.matrix(decathlon2.active), group = "data", dataset = "decathlon",
                     transp = FALSE,
                     overwriteFile = TRUE, overwriteDataset = TRUE, 
                     unlimited = FALSE)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdPCA_hdf5("test_temp.hdf5", "data", "decathlon", ncomponents = 0L, bcenter = TRUE, bscale = TRUE, k = 2L, q = 1L, rankthreshold = 0.0, SVDgroup = NULL, overwrite = TRUE, threads = NULL) 

res_ind.dist <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "PCA/decathlon/ind.dist"); res_ind.dist
res_components <-  h5read("test_temp.hdf5", "PCA/decathlon/components"); res_components[1:5,1:5]
res_ind.coord <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "PCA/decathlon/ind.coord"); res_ind.coord[1:5,1:5]
res_ind.cos2 <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "PCA/decathlon/ind.cos2"); res_ind.cos2[1:5,1:5]
res_ind.contrib <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "PCA/decathlon/ind.contrib"); res_ind.contrib[1:5,1:5]


res_lambda <-  h5read("test_temp.hdf5", "/PCA/decathlon/lambda")
res_var <-  h5read("test_temp.hdf5", "PCA/decathlon/variance")
res_cumvar <-  h5read("test_temp.hdf5", "PCA/decathlon/cumvar")
res_var.coord <-  h5read("test_temp.hdf5", "PCA/decathlon/var.coord")
res_var.cos2 <-  h5read("test_temp.hdf5", "PCA/decathlon/var.cos2")







matrix <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/NORMALIZED_T/data/decathlon")
u <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/SVD/decathlon/u")
v <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/SVD/decathlon/v")
d <-  h5read("/Users/mailos/PhD/dummy/test_temp.hdf5", "/SVD/decathlon/d")

sweep( crossprod(as.matrix(matrix), u ),2,d,FUN="/")


resr <- PCA(decathlon2.active, scale.unit = TRUE, ncp = 5, graph = FALSE)


all.equal()




decathlon2.active
scale(decathlon2.active, scale = TRUE, center = FALSE)
scale(decathlon2.active, scale = FALSE, center = TRUE)

scale(scale(decathlon2.active, scale = FALSE, center = TRUE), scale = TRUE, center = FALSE)

svd(scale(decathlon2.active, scale = TRUE, center = TRUE))


# setwd(/Users/mailos/PhD/TREBALLANT/BigDataStatMethAPIsV3)

setwd("/Users/mailos/PhD/dummy")
bdSVD_hdf5("test_temp.hdf5", "data", "decathlon", bcenter = TRUE, bscale = TRUE, k = 2L, q = 1L, force = TRUE) 
