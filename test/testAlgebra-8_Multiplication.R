library("BigDataStatMeth")
library("rhdf5")
# devtools::reload(pkgload::inst("BigDataStatMeth"))

# setwd("/Users/mailos/PhD/dummy")
setwd("/scratch/dpelegri/BigDataStatMeth/Benchmark/hdf5_files")

N = 1500
M = 1500

set.seed(555)
a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
b <- matrix( rnorm( N*M, mean=0, sd=1), M, N) 



# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = a, group = "pepet", 
                     dataset = "A",
                     transp = FALSE,
                     overwriteFile = TRUE, 
                     overwriteDataset = FALSE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "test_temp.hdf5", 
                     object = b, 
                     group = "pepet", 
                     dataset = "B",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)


bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "A", B = "B", 
                 outgroup = "results", outdataset = "res", 
                 force = TRUE, paral = TRUE, threads = 4 ,block_size = 1024)


times <- microbenchmark::microbenchmark(
    R = a%*%b,
    hdf5 =  bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "A", B = "B", 
                             outgroup = "results", outdataset = "res", 
                             force = TRUE, paral = TRUE, threads = 4),
    hdf5B =  bdblockmult_hdf5(filename = "test_temp.hdf5",group = "pepet", A = "A", B = "B", 
                             outgroup = "results", outdataset = "res", 
                             force = TRUE, paral = TRUE, threads = 4 ,block_size = 1024),
    mem = bdblockMult(a, b, paral = FALSE), # Runs with blocks nthreads: 4
    memB = bdblockMult(a, b, paral = TRUE, threads = 4, block_size = 1024), # Runs with blocks nthreads: 4
    times = 3, unit = "s")
times


# Checks

file <- "test_temp.hdf5"
dataset <- "results/res"
resr <- a%*%b
res <-  h5read(file,dataset)

all.equal( res, resr)
