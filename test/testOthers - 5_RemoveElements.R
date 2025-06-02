
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
datasets <- c("pepet/dataseta", "pepet/datasetc")

bdgetDatasetsList_hdf5(filename = "testBind.hdf5", group = "pepet")

bdRemove_hdf5_element(filename = "testBind.hdf5", elements = datasets)

bdgetDatasetsList_hdf5(filename = "testBind.hdf5", group = "pepet")



