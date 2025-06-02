
## ########################################################################## ##
##   Check Reduce datasets inside hdf5 data file
## ########################################################################## ##

# Creem un grup amb diferents matrius on aplicarem el reduce '+' i '-'



library("BigDataStatMeth")
library(rhdf5)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

N = 20
M = 10

set.seed(555)
a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
set.seed(111)
b <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 

N = 30
M = 10
set.seed(222)
c <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
set.seed(333)
d <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 



# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "testReduce.hdf5", 
                     object = a, group = "pepet", 
                     dataset = "dataseta",
                     transp = FALSE,
                     overwriteFile = TRUE, 
                     overwriteDataset = FALSE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "testReduce.hdf5", 
                     object = b, group = "pepet", 
                     dataset = "datasetb",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "testReduce.hdf5", 
                     object = c, group = "pepet", 
                     dataset = "datasetc",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "testReduce.hdf5", 
                     object = d, group = "pepet", 
                     dataset = "datasetd",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)



# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdReduce_hdf5_dataset(filename = "testReduce.hdf5", group = "pepet", 
                      reducefunction = "+", 
                      outgroup = "pepetR", outdataset = "reduit", 
                      overwrite = TRUE, remove = TRUE)



# Check
##########

file <- "testReduce.hdf5"
dataset <- "pepetR/reduit"
resr <- a+b+c+d
res <-  h5read(file,dataset)
res

all.equal( res, resr)
all.equal( t(res), resr)


a+b+c+d
