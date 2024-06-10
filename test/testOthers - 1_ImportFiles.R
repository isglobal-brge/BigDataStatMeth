## ########################################################################## ##
##   Check text file import
## ########################################################################## ##

library(rhdf5)
library(BigDataStatMeth)
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("/Users/mailos/PhD/dummy")

bdImportTextFile_hdf5(filename = "colesterol.csv",
                      outputfile = "colesterol_file.hdf5", 
                      outGroup = "COLESTEROL", 
                      outDataset = "COLESTEROLDATA", 
                      sep = ",", 
                      header = TRUE, rownames = FALSE, 
                      overwrite = TRUE)



## ----read_data_col_HDF5-------------------------------------------------------

file <- "/Users/mailos/PhD/dummy/colesterol_file.hdf5"
dataset <- "COLESTEROL/COLESTEROLDATA"
res <-  h5read(file,dataset)

colesteroldata <- read.table("/Users/mailos/PhD/dummy/colesterol.csv", header = T, sep = ",")

h5ls("colesterol_file.hdf5")

# We can open the file and have access to the data
h5Colesterol <- H5Fopen("colesterol_file.hdf5")

# Show hdf5 content dataset
h5Colesterol$COLESTEROL$COLESTEROLDATA[1:5, 1:6]

# Show colnames
head(h5Colesterol$COLESTEROL$.COLESTEROLDATA_dimnames$`2`)

H5Fclose(h5Colesterol)



# ####################################################
# ### --- TETEJAR IMPORT AMB ROWNAMES I COLNAMES
# ####################################################

colesteroldata <- read.table("/Users/mailos/PhD/dummy/colesterol.csv", header = T, sep = ",")
rownames(colesteroldata) <- paste0("sampleId_", seq_len(nrow(colesteroldata)) )
write.table(colesteroldata, file = "/Users/mailos/PhD/dummy/colesterol_rownames.csv", sep = ",", row.names = TRUE)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
setwd("/Users/mailos/PhD/dummy")

bdImportTextFile_hdf5(filename = "colesterol_rownames.csv",
                      outputfile = "colesterol_file.hdf5", 
                      outGroup = "COLESTEROL_WR", 
                      outDataset = "COLESTEROLDATA", 
                      sep = ",", 
                      header = TRUE, rownames = TRUE, 
                      overwrite = TRUE)


# ##############################################################################
# ### --- TETEJAR IMPORT AMB ROWNAMES I COLNAMES 
#           + fitxers mes grans on faci particions x blocks de 1000 registres
# ##############################################################################

colesteroldata <- read.table("/Users/mailos/PhD/dummy/colesterol.csv", header = T, sep = ",")

colesterolgran <- rbind(colesteroldata, colesteroldata)

rownames(colesterolgran) <- paste0("sampleId_", seq_len(nrow(colesterolgran)) )
write.table(colesterolgran, file = "/Users/mailos/PhD/dummy/colesterol_rc_gran.csv", sep = ",", row.names = TRUE)

# devtools::reload(pkgload::inst("BigDataStatMeth"))
setwd("/Users/mailos/PhD/dummy")

bdImportTextFile_hdf5(filename = "colesterol_rc_gran.csv",
                      outputfile = "colesterol_file.hdf5", 
                      outGroup = "COLESTEROL_RC_GRAN", 
                      outDataset = "COLESTEROLDATA", 
                      sep = ",", 
                      header = TRUE, rownames = TRUE, 
                      overwrite = TRUE)




# ##############################################################################
# ### --- TETEJAR ESPAIS EN BLANC + ROWNAMES + COLNAMES
#           + fitxers mes grans on faci particions x blocks de 1000 registres
#           
#  Link de descárrega del fitxer: 
#   Enllaç: https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q4&filename=CRISPRInferredModelGrowthRate.csv
#   Web: https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q4
# ##############################################################################

# devtools::reload(pkgload::inst("BigDataStatMeth"))
setwd("/Users/mailos/PhD/dummy")

bdImportTextFile_hdf5(filename = "CRISPR.csv",
                      outputfile = "colesterol_file.hdf5", 
                      outGroup = "CRISPR", 
                      outDataset = "crisprData", 
                      sep = ",", 
                      header = TRUE, rownames = TRUE, 
                      overwrite = TRUE)

