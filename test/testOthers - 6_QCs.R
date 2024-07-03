library(BigDataStatMeth)
library(rhdf5)

# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

# CREATE DATA
set.seed(108432)
geno.sim <- matrix(sample(0:3, 10000, replace = TRUE), byrow = TRUE, ncol = 20)
bdCreate_hdf5_matrix(filename = "delayed.hdf5", object = geno.sim, group = "OMICS", dataset = "geno", overwriteFile = TRUE)

file <- "delayed.hdf5"
dataset <- "OMICS/geno"
geno <-  h5read(file,dataset)

geno[1:5,1:10]


# QC - IMPUTE DATA

bdImpute_snps_hdf5(filename = "delayed.hdf5", group="OMICS", dataset="geno",
                   outgroup="OMICS", outdataset="imputedgeno", overwrite = TRUE)

# Get imputed data and show the first 5 rows
file <- "delayed.hdf5"
dataset <- "OMICS/imputedgeno"
imputedgeno <- h5read(file,dataset)

imputedgeno[1:5,1:10]


## QC - REMOVE LOW DATA


### by Cols
bdRemovelowdata_hdf5("delayed.hdf5", group="OMICS", dataset="geno",
                     outgroup="OMICS", outdataset="removedgeno",
                     bycols = TRUE, pcent = 0.25, overwrite = TRUE)


# Get data and show the first 5 rows
file <- "delayed.hdf5"
dataset <- "OMICS/removedgeno"
geno <-  h5read(file,dataset)

geno[1:5,]

dim(geno.sim)
dim(geno)


### by Rows
bdRemovelowdata_hdf5("delayed.hdf5", group="OMICS", dataset="geno",
                     outgroup="OMICS", outdataset="removedgeno",
                     bycols = FALSE, pcent = 0.25, overwrite = TRUE)


# Get data and show the first 5 rows
file <- "delayed.hdf5"
dataset <- "OMICS/removedgeno"
geno <-  h5read(file,dataset)

geno[1:5,]

dim(geno.sim)
dim(geno)




########### CHECKEJAT FINS AQUÍ !!!!!!!!!!











## NORMALIZE DATA

bdNormalize_hdf5(filename = "delayed.hdf5", group="OMICS", dataset="imputedgeno", bcenter = TRUE, bscale = TRUE, force = TRUE)

bdNormalize_hdf5()

# Geno after normalization
h5fsvd = H5Fopen("delayed.hdf5")
genonormalized <- h5fsvd$NORMALIZED$OMICS$geno
h5closeAll()

genonormalized[1:5,1:10]




# Get imputed data and show the first 5 rows
h5fsvd = H5Fopen("delayed.hdf5")
removedgeno <- h5fsvd$OMICS$removedgeno
h5closeAll()

removedgeno[1:5,1:10]