## ########################################################################## ##
##   Check Sort datasets inside hdf5 data file
## ########################################################################## ##

# Utilitzem el fitxer creat amb mgcca per no haver de crear-ne un de nou, així també
# testejem i comprovem que ho ordena tal i com s'esperaria per a funcionar amb MGCCA.
# Hi ha una còpia del fitxer a: 
#       /Users/mailos/PhD/dummy/data/mgccsChecks/

library(rhdf5)
library(BigDataStatMeth)

# devtools::reload(pkgload::inst("BigDataStatMeth"))

load("~/PhD/dummy/data/mgccsChecks/dadesFiltrades.RData")

setwd("/Users/mailos/PhD/dummy")

# devtools::reload(pkgload::inst("BigDataStatMeth"))
BigDataStatMeth:::bdSort_hdf5_dataset( filename = "gettables.hdf5", 
                                       group = "NORMALIZED/MGCCA_IN",
                                       dataset = "ACC_RNASeq2GeneNorm-20160128",
                                       outdataset = "ACC_RNASeq2GeneNorm-20160128", 
                                       outgroup = "Sorted",
                                       blockedSortlist = blocks,
                                       func = "sortRows",
                                       overwrite = TRUE)




