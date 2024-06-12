
## ########################################################################## ##
##   Check Bind datasets inside hdf5 data file
## ########################################################################## ##

# Creem un grup amb diferents matrius on aplicarem el bind per files i per columnes



library("BigDataStatMeth")
library(rhdf5)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

N = 30
M = 20

set.seed(555)
a <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
set.seed(111)
b <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 

N = 30
M = 20
set.seed(222)
c <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 
set.seed(333)
d <- matrix( rnorm( N*M, mean=0, sd=1), N, M) 



# devtools::reload(pkgload::inst("BigDataStatMeth"))
bdCreate_hdf5_matrix(filename = "testBind.hdf5", 
                     object = a, group = "pepet", 
                     dataset = "dataseta",
                     transp = FALSE,
                     overwriteFile = TRUE, 
                     overwriteDataset = FALSE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "testBind.hdf5", 
                     object = b, group = "pepet", 
                     dataset = "datasetb",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "testBind.hdf5", 
                     object = c, group = "pepet", 
                     dataset = "datasetc",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)

bdCreate_hdf5_matrix(filename = "testBind.hdf5", 
                     object = d, group = "pepet", 
                     dataset = "datasetd",
                     transp = FALSE,
                     overwriteFile = FALSE, 
                     overwriteDataset = TRUE, 
                     unlimited = FALSE)



# ==>   byRows
#################

# devtools::reload(pkgload::inst("BigDataStatMeth"))
datasets <- c("dataseta", "datasetb", "datasetc", "datasetd")
bdBind_hdf5_datasets(filename = "testBind.hdf5", group = "pepet", 
                     datasets = datasets, 
                     outgroup = "pepetbind", outdataset = "FullpepetR", 
                     func = "bindRows", overwrite = TRUE)


# Check - byRows

file <- "testBind.hdf5"
dataset <- "pepetbind/FullpepetR"
resr <- rbind( a, b, c, d)
res <-  h5read(file,dataset)


all.equal( res, resr)



# ==>   byCols
#################

# devtools::reload(pkgload::inst("BigDataStatMeth"))
datasets <- c("dataseta", "datasetb", "datasetc", "datasetd")
bdBind_hdf5_datasets(filename = "testBind.hdf5", group = "pepet", 
                     datasets = datasets, 
                     outgroup = "pepetbind", outdataset = "FullpepetC", 
                     func = "bindCols", overwrite = TRUE)


# Check - byRows

file <- "testBind.hdf5"
dataset <- "pepetbind/FullpepetC"
resr <- cbind( a, b, c, d)
res <-  h5read(file,dataset)


all.equal( res, resr)

