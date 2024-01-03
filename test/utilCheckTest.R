library(rhdf5)
library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

N = 1000
M = 800

set.seed(555)
a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
b <- matrix( rnorm( N*M, mean=0, sd=1), M, N) 
# b <- t(b)
v <- rnorm( N, mean=0, sd=1)



# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = a, 
                     group = "pepet", 
                     dataset = "datasetpepet",
                     transp = FALSE,
                     overwriteFile = TRUE, 
                     overwriteDataset = FALSE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = b, 
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


blockSum_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 5, force = TRUE ) 

file <- "test_temp.hdf5"
dataset <- "results/res"
# dataset <- "pepet/datasetpepet"

resr <- a + a
res <-  h5read(file,dataset)

all.equal( round( res, 5), round( resr, 5))





blockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "resSparse", force = TRUE,block_size = 5)


# devtools::reload(pkgload::inst("BigDataStatMeth"))
blockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "resSparse", force = TRUE, paral = TRUE, threads = 4, block_size = 5 )

# bdblockmult_sparse_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "resSparse2", force = TRUE, block_size = 5  )
# blockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", B = "tdatasetpepet", outgroup = "results", outdataset = "resSparse", force = TRUE,block_size = 5)


# bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet",  outgroup = "results", outdataset = "res", force = TRUE ) # Crossprod A %*% t(A)
# bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 1024, force = TRUE ) # Crossprod A %*% t(A)
# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 140, force = TRUE ) # Crossprod A %*% t(A)
# bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", B = "tdatasetpepet", force = TRUE, block_size = 140 ) # Crossprod A %*% t(B)

file <- "test_temp.hdf5"
dataset <- "results/res"
dataset <- "pepet/datasetpepet"

resr <- crossprod(a)
res <-  h5read(file,dataset)

all.equal( round( res, 5), round( resr, 5))





res <- rbenchmark::benchmark( "R" = crossprod(a, b),
                              "BDSM" =  bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", B = "tdatasetpepet", force = TRUE, block_size = 140 ), # Crossprod A %*% t(B)
                              "R_a" = crossprod(a),
                              "BDSM_a" =  bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 140, force = TRUE ) ,
                              replications = 5
                              # columns = c("test", "replications", "elapsed","relative", "user.self", "sys.self")
                              )
print(res)



res <- microbenchmark::microbenchmark( "R" = crossprod(a, b),
                              "BDSM" =  bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", B = "tdatasetpepet", force = TRUE, block_size = 140 ), # Crossprod A %*% t(B)
                              "R_a" = crossprod(a),
                              "BDSM_a" =  bdCrossprod_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "datasetpepet", outgroup = "results", outdataset = "res", block_size = 140, force = TRUE ) ,
                              times = 5
                              # columns = c("test", "replications", "elapsed","relative", "user.self", "sys.self")
                            )
print(res)

