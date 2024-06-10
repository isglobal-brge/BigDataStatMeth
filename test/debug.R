devtools::clean_dll()
devtools::load_all()


# devtools::reload(pkgload::inst("BigDataStatMeth"))
library(BigDataStatMeth)
library(rhdf5)
N <- 40; M <- 20

set.seed(555); mat <- matrix( rnorm( N*M, mean=0, sd=10), N, M) 
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



SumR <- mat + vect
SumHDF5_mv = bdblockSum_hdf5(filename = fname, group = gname,  A = dname, B = vname, 
                             outgroup = "results", outdataset = "res",  threads = 6, force = TRUE )