library(rhdf5)
library(microbenchmark)
library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

N <- 10000
M <- 10000
nc <-  4
ntimes <- 2

set.seed(555)
mat <- matrix( rnorm( N*M, mean=0, sd=10), N, M) 
set.seed(555)
vect <- rnorm( N, mean=0, sd=10)


# Matrix-Matrix and Matrix-Vector Multiplication




# Matrix - Matrix Multiplication

# devtools::reload(pkgload::inst("BigDataStatMeth"))
times <- microbenchmark::microbenchmark(
    Mult_R = mat %*% mat, # Runs with blocks nthreads: 2
    Mult_mem_blocks_NP = BigDataStatMeth::bdblockMult(mat, mat, paral = FALSE, threads = 1, block_size = 1000), # Runs with blocks nthreads: 2
    Mult_mem_blocks_1 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 1, block_size = 1000), # Runs with blocks nthreads: 2
    Mult_mem_blocks_2_1k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 2, block_size = 1000), # Runs with blocks nthreads: 4
    Mult_mem_blocks_2_2k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 2, block_size = 2000), # Runs with blocks nthreads: 4
    Mult_mem_blocks_2_3k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 2, block_size = 3000), # Runs with blocks nthreads: 2
    Mult_mem_blocks_2_5k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 2, block_size = 5000), # Runs with blocks nthreads: 2
    Mult_mem_blocks_3 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 3, block_size = 1000), # Runs with blocks nthreads: 4
    Mult_mem_blocks_4_1k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 4, block_size = 1000), # Runs with blocks nthreads: 4
    Mult_mem_blocks_4_2k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 4, block_size = 2000), # Runs with blocks nthreads: 4
    Mult_mem_blocks_4_25k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 4, block_size = 2500), # Runs with blocks nthreads: 4
    Mult_mem_blocks_4_3k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 4, block_size = 3000), # Runs with blocks nthreads: 4
    Mult_mem_blocks_4_5k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 4, block_size = 5000), # Runs with blocks nthreads: 4
    Mult_mem_blocks_6_1k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 6, block_size = 1000), # Runs with blocks nthreads: 4
    Mult_mem_blocks_6_2k = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 6, block_size = 2000), # Runs with blocks nthreads: 4
    times = ntimes, unit = "s")
times



# devtools::reload(pkgload::inst("BigDataStatMeth"))
res2 <- BigDataStatMeth::bdblockMult(mat, mat, paral = FALSE, threads = 1) # Runs with blocks nthreads: 4
res1 <- BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 4, block_size = 2000) # Runs with blocks nthreads: 4
res3 <- BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 1, block_size = 2000) # Runs with blocks nthreads: 4

resR <- mat %*% mat
all.equal(resR, res1)
all.equal(resR, res2)
all.equal(resR, res3)


# Matrix - Vector Multiplication

# devtools::reload(pkgload::inst("BigDataStatMeth"))
times <- microbenchmark::microbenchmark(
    MultR = mat * vect, 
    Mult_mem_NP = BigDataStatMeth::bdblockMult(mat, vect, paral = FALSE), # Runs with blocks nthreads: 2
    Mult_mem_blocks_1 = BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = 1), # Runs with blocks nthreads: 2
    Mult_mem_blocks_2 = BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = 2), # Runs with blocks
    Mult_mem_blocks_3 = BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = 3), # Runs with blocks
    Mult_mem_blocks_4 = BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = 4), # Runs with blocks
    Mult_mem_blocks_5 = BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = 5), # Runs with blocks
    Mult_mem_blocks_6 = BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = 6), # Runs with blocks
    Mult_mem = BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = nc, byBlocks = FALSE), # Forced to run without blocks
    times = ntimes, unit = "s")
times



# devtools::reload(pkgload::inst("BigDataStatMeth"))
res2 <- BigDataStatMeth::bdblockMult(mat, vect, paral = FALSE, threads = 1) # Runs with blocks nthreads: 4
res1 <- BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = 4, block_size = 2000) # Runs with blocks nthreads: 4
res3 <- BigDataStatMeth::bdblockMult(mat, vect, paral = TRUE, threads = 1, block_size = 2000) # Runs with blocks nthreads: 4
res4 <- BigDataStatMeth::bdblockMult(mat, vect, paral = FALSE, threads = 1, block_size = 2000) # Runs with blocks nthreads: 4

resR <- mat * vect
all.equal(resR, res1)
all.equal(resR, res2)
all.equal(resR, res3)
all.equal(resR, res3)






# Matrix - Vector Sum
times <- microbenchmark::microbenchmark(
    SumR = mat %*% mat, 
    Sum_mem_blocks_1 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 1), # Runs with blocks nthreads: 2
    Sum_mem_blocks_2 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 2), # Runs with blocks
    Sum_mem_blocks_3 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 3), # Runs with blocks
    Sum_mem_blocks_4 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 4), # Runs with blocks
    Sum_mem_blocks_5 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 5), # Runs with blocks
    Sum_mem_blocks_6 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 6), # Runs with blocks
    Sum_mem = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = nc, byBlocks = FALSE), # Forced to run without blocks
    times = ntimes, unit = "s")
times

# devtools::reload(pkgload::inst("BigDataStatMeth"))
times <- microbenchmark::microbenchmark(
    Sum_mem_blocks_1 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 1), # Runs with blocks nthreads: 2
    Sum_mem_blocks_2 = BigDataStatMeth::bdblockMult(mat, mat, paral = TRUE, threads = 2), # Runs with blocks nthreads: 2
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
