library(rhdf5)
library(microbenchmark)
library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))


N <- 16000
M <- 16000
nc <-  4
ntimes <- 2

set.seed(555)
mat <- matrix( rnorm( N*M, mean=0, sd=10), N, M) 
set.seed(555)
vect <- rnorm( N, mean=0, sd=10)


fname <- "/Users/mailos/PhD/dummy/benchmark/BasicFunctions_tmp.hdf5"
gname <- "InputDatasets"
dname <- "dataset" # Matrix

vname <- "vdataset" # Vector

bdCreate_hdf5_matrix(filename = fname, 
                     object = mat, group = gname, 
                     dataset = dname,
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)   


bdCreate_hdf5_matrix(filename = fname, 
                     object = vect, group = gname, 
                     dataset = vname,
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)   




# Matrix - Vector Sum
times <- microbenchmark::microbenchmark(
    SumR = mat + mat, 
    Sum_mem_blocks = BigDataStatMeth::bdblockSum(mat, mat, paral = TRUE, threads = nc), # Runs with blocks
    Sum_mem = BigDataStatMeth::bdblockSum(mat, mat, paral = TRUE, threads = nc, byBlocks = FALSE), # Forced to run without blocks
    SumR_mv = mat + vect, 
    Sum_mv = BigDataStatMeth::bdblockSum(mat, vect, paral = FALSE, threads = nc),
    Sum_vm = BigDataStatMeth::bdblockSum(vect, mat, paral = FALSE, threads = nc),
    Sum_mv_blocks = BigDataStatMeth::bdblockSum(mat, vect, paral = TRUE, threads = nc),
    Sum_vm_blocks = BigDataStatMeth::bdblockSum(vect, mat, paral = TRUE, threads = nc),
    SumR_v = vect + vect, 
    Sum_v = BigDataStatMeth::bdblockSum(vect, vect, paral = TRUE, threads = nc),
    SumHDF5_2 =   bdblockSum_hdf5(filename = fname, group = gname,  A = dname, B = dname, 
                                outgroup = "results", outdataset = "res",  threads = 2, force = TRUE ),
    SumHDF5_4 =   bdblockSum_hdf5(filename = fname, group = gname,  A = dname, B = dname, 
                                outgroup = "results", outdataset = "res",  threads = 4, force = TRUE ),
    SumHDF5_6 =   bdblockSum_hdf5(filename = fname, group = gname,  A = dname, B = dname, 
                                outgroup = "results", outdataset = "res",  threads = 6, force = TRUE ),
    times = ntimes, unit = "s")
times



# Matrix - Vector Substract

times <- microbenchmark::microbenchmark(
    SubR = mat - mat, 
    Sub_mem_blocks = BigDataStatMeth::bdblockSubstract(mat, mat, paral = TRUE, threads = nc), # Runs with blocks
    Sub_mem = BigDataStatMeth::bdblockSubstract(mat, mat, paral = TRUE, threads = nc, byBlocks = FALSE), # Forced to run without blocks
    SubR_mv = mat + vect, 
    Sub_mv = BigDataStatMeth::bdblockSubstract(mat, vect, paral = FALSE, threads = nc),
    Sub_vm = BigDataStatMeth::bdblockSubstract(mat, vect, paral = FALSE, threads = nc),
    Sub_mv_blocks = BigDataStatMeth::bdblockSubstract(mat, vect, paral = TRUE, threads = nc),
    Sub_vm_blocks = BigDataStatMeth::bdblockSubstract(vect, mat, paral = TRUE, threads = nc),
    SubR_v = vect + vect, 
    Sub_v = BigDataStatMeth::bdblockSubstract(vect, vect, paral = TRUE, threads = nc),
    SumHDF5_2 =   bdblockSubstract_hdf5(filename = fname, group = gname,  A = dname, B = dname, 
                                  outgroup = "results", outdataset = "res",  threads = 2, force = TRUE ),
    SumHDF5_4 =   bdblockSubstract_hdf5(filename = fname, group = gname,  A = dname, B = dname, 
                                  outgroup = "results", outdataset = "res",  threads = 4, force = TRUE ),
    SumHDF5_6 =   bdblockSubstract_hdf5(filename = fname, group = gname,  A = dname, B = dname, 
                                  outgroup = "results", outdataset = "res",  threads = 6, force = TRUE ),
    times = ntimes, unit = "s")
times



# --- Checks


# Memory
library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

nc = 2
N <- 1000
M <- 800

M <- 1000
N <- 800

set.seed(555)
mat <- matrix( rnorm( N*M, mean=0, sd=10), N, M) 
set.seed(555)
vect <- rnorm( N, mean=0, sd=10)

SumR_mv = mat + vect
Sum_mv =  BigDataStatMeth::bdblockSum(mat, vect, paral = FALSE, threads = nc)
Sum_mv_B = BigDataStatMeth::bdblockSum(mat, vect, paral = TRUE, threads = nc)

all.equal(SumR_mv, Sum_mv)
all.equal(SumR_mv, Sum_mv_B)


SumR <-  mat + mat
Sum_mem_blocks <-  BigDataStatMeth::bdblockSum(mat, mat, paral = TRUE, threads = nc) # Runs with blocks
Sum_mem = BigDataStatMeth::bdblockSum(mat, mat, paral = TRUE, threads = nc, byBlocks = FALSE) # Forced to run without blocks

all.equal(SumR[1:5000,], Sum_mem_blocks[1:5000,])
all.equal(SumR[1:5000,], Sum_mem[1:5000,])
all.equal(SumR[15000:16000,], Sum_mem[15000:16000,])

Sum_v <-  BigDataStatMeth::bdblockSum(vect, vect, paral = TRUE, threads = nc)
SumR_v <-  vect + vect

all.equal(SumR_v[1:5000], Sum_v[1:5000])



# HDF5 data files

library(rhdf5)
library(microbenchmark)
library(BigDataStatMeth)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
nc = 2
N <- 1000
M <- 800

M <- 1000
N <- 800

set.seed(555); mat <- matrix( rnorm( N*M, mean=0, sd=10), N, M) 
set.seed(555)
vect <- rnorm( M, mean=0, sd=10)


fname <- "/Users/mailos/PhD/dummy/benchmark/BasicFunctions_tmp.hdf5"
gname <- "InputDatasets"
dname <- "dataset" # Matrix

vname <- "vdataset" # Vector

bdCreate_hdf5_matrix(filename = fname, 
                     object = mat, group = gname, 
                     dataset = dname,
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)   


bdCreate_hdf5_matrix(filename = fname, 
                     object = vect, group = gname, 
                     dataset = vname,
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)   



SumR_mv <- mat + vect
SumHDF5_mv = bdblockSum_hdf5(filename = fname, group = gname,  A = dname, B = vname, 
                             outgroup = "results", outdataset = "res",  threads = 6, force = TRUE )

res <-  h5read(fname, "/results/res")
all.equal(SumR, res)

SubsR <- mat - vect
SubsHDF5_mv = bdblockSubstract_hdf5(filename = fname, group = gname,  A = dname, B = vname, 
                             outgroup = "results", outdataset = "res",  threads = 6, force = TRUE )

res <-  h5read(fname, "/results/res")
all.equal(SubsR, res)



SubsR_mm <- mat - mat
SubsHDF5_mv = bdblockSubstract_hdf5(filename = fname, group = gname,  A = dname, B = dname, 
                                    outgroup = "results", outdataset = "res",  threads = 6, force = TRUE )

res_mm <-  h5read(fname, "/results/res")
all.equal(SubsR_mm, res_mm)


SumR_mm <- mat + mat
SubsHDF5_m = bdblockSum_hdf5(filename = fname, group = gname,  A = dname, B = dname,
                              outgroup = "results", outdataset = "res",  threads = 6, force = TRUE )

res_mm <-  h5read(fname, "/results/res")
all.equal(SumR_mm, res_mm)
