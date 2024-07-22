library(gdsfmt)
library(DelayedArray)
library(SeqArray)
library(BigDataStatMeth)
library(rhdf5)
library(grDevices)


library(ggplot2)
library(ggthemes)
library(rasterpdf)


## ##########################################
## 0.- Install BigDataStatMeth package
## ##########################################

#..# library(devtools)
#..# install_github("isglobal-brge/BigDataStatMeth")

## ##########################################
## 0.- Create paths
## ##########################################

dir.create("/scratch/dpelegri/working/InversionCalling/results/data/",  recursive = TRUE, showWarnings = FALSE)
dir.create("/home/isglobal.lan/dpelegri/working/PhD/BigDataStatMeth/Analysis/InversionCalling/plots",  recursive = TRUE, showWarnings = FALSE)


## ##########################################
## 1.- Load gds file as DelayedMatrix
## ##########################################

#.. SERVER 05 ..# setwd("/home/isglobal.lan/dpelegri/homews/results") #.. SERVER 05 ..# 
setwd("/home/isglobal.lan/dpelegri/working/PhD/BigDataStatMeth/Analysis/InversionCalling") #.. Yamabuki Server ..# 
filename <- "/home/isglobal.lan/dpelegri/working/TFM/UkbBiobank/chr8inv.gds"

# Open gds file
gds<-gdsfmt::openfn.gds(filename)

# Get genotype data
geno <- index.gdsn(gds, "genotype")
feno <- index.gdsn(gds, "sample.annot/phenotype")
# genotype <- DelayedArray( read.gdsn( geno ))
genotype <- read.gdsn( geno )
fenotype <- read.gdsn(feno)

# Set colnames and rownames
rownames(genotype) <- read.gdsn(index.gdsn(gds, "sample.id"))
colnames(genotype) <- read.gdsn(index.gdsn(gds, "snp.id"))

# Close gds file
gdsfmt::closefn.gds(gds)


## ##########################################
## 2.- Write DelayedMatrix as hdf5 Data File
## ##########################################

print(paste0(format(Sys.time(), "%a %b %d %X %Y - "), "Reading data from GDS File and writing to hdf5 file...."))
# bdCreate_hdf5_matrix_file(filename = "/scratch/dpelegri/working/InversionCalling/results/data/chr8.hdf5", object = genotype, group = "invs", dataset = "geno")

bdCreate_hdf5_matrix( filename = "/scratch/dpelegri/working/InversionCalling/results/data/chr8.hdf5", 
          object = genotype, group = "invs", dataset = "geno", overwriteFile = TRUE, overwriteDataset = TRUE)


## ##########################################
## 3.- Remove NA with percent NA>5%
## ##########################################

# Remove Samples with more or equal to 5 % of missing data

print(paste0(format(Sys.time(), "%a %b %d %X %Y - "), "Removing Samples with 5% missing data ... "))
start_time <- Sys.time()

remsamples <- bdRemovelowdata_hdf5("/scratch/dpelegri/working/InversionCalling/results/data/chr8.hdf5", 
                              group = "invs", 
                              dataset = "geno", 
                              outgroup = "invs", 
                              outdataset = "genofilter", 
                              bycols = FALSE, pcent = 0.05)

end_time <- Sys.time()
paste0("Time spent removing Samples ... ", end_time - start_time)
print(paste("Removed Samples : ", remsamples))


# chr8 :  32069 Samples Removed (SNPs) with pcent = 0.05
# chr8   122962 Samples Removed (SNPs) with pcent = 0.03
# chr8 :  140 Samples Removed (Individuals)


## ##########################################
## 3.2 - Remove MAF 5%
## ##########################################

# Remove Samples with more or equal to 10 % of missing data

print(paste0(format(Sys.time(), "%a %b %d %X %Y - "), "Removing SNPs with maf 5% ... "))
start_time <- Sys.time()

ressnps <- bdRemoveMAF_hdf5 ("/scratch/dpelegri/working/InversionCalling/results/data/chr8.hdf5", 
                                 group = "invs", 
                                 dataset = "genofilter", 
                                 outgroup = "invs", 
                                 outdataset = "genofilter_maf", 
                                 blocksize = 1000,
                                 bycols = TRUE, maf = 0.05, 
                                 overwrite = TRUE)
end_time <- Sys.time()
print(paste("Time spent removing SNPs : ", end_time - start_time))


print(paste("Removed SNPs : ", ressnps))


## ##########################################
## 4.- Impute Missing Data
## ##########################################

print(paste0(format(Sys.time(), "%a %b %d %X %Y - "), "Imputing missings..... "))
start_time <- Sys.time()

bdImputeSNPs_hdf5( "/scratch/dpelegri/working/InversionCalling/results/data/chr8.hdf5", 
                    group = "invs", 
                    dataset = "genofilter_maf", 
                    bycols = TRUE, 
                    outgroup = "invs", 
                    outdataset = "genofilter_maf",
                    overwrite = TRUE)


end_time <- Sys.time()
paste0("Time spent imputing Samples ... ", end_time - start_time)

## ##########################################
## 5.- Perform SVD decomposition from QC Data
## ##########################################

print(paste0(format(Sys.time(), "%a %b %d %X %Y - "), "Performing SVD..... "))
start_time <- Sys.time()

res <- bdSVD_hdf5(file = "/scratch/dpelegri/working/InversionCalling/results/data/chr8.hdf5", 
                  group = "invs", 
                  dataset = "genofilter_maf", 
                  k = 40, 
                  threads = 8, 
                  bcenter = TRUE, 
                  bscale = TRUE)

end_time <- Sys.time()
paste0("Time spent performing SVD ... ", end_time - start_time)

print(paste0(format(Sys.time(), "%a %b %d %X %Y - "), ".....SVD Do it "))


# Times : 
#     Start : dc mai 12 08:31:01 2021 - Performing SVD.....
#     End :   dc mai 12 10:38:20 2021 - .....SVD Do it

# A PARTIR D'AQUÍ, TOT COMENTAT.....


## ##########################################
## 6.- Get data from hdf5 data file
## ##########################################

# Examine hierarchy before open file
h5ls(res$file)




# Open file and get data, all data is stored under SVD group
#..# h5fsvd = H5Fopen("/scratch/dpelegri/working/InversionCalling/results/data/chr8.hdf5")
h5fsvd = H5Fopen(res$file)

v <- h5fsvd$SVD$genofilter_maf$v
d <- h5fsvd$SVD$genofilter_maf$d
u <- h5fsvd$SVD$genofilter_maf$u

colnames(u) <- sprintf("PC%s",seq(1:dim(u)[2]))
colnames(v) <- sprintf("PC%s",seq(1:dim(v)[2]))

# dim(v)
dim(u)

h5closeAll()



## #######################################################
## 7.- Analyze PCA (colzes --> percent Variança acumulada)
## #######################################################

# Get variance explained by principal components, the chosen factors should explain at least 70 to 80% of variance.
rho <- d^2/sum(d^2)
cumsum(rho[1:10])

print(paste("Cumulated variance  10 first PCs : ", cumsum(rho[1:10])))

# First 2 components explains...%


fname <- paste0("plots/chr8_cumvar.png")
png(fname, width = 8, height = 8, units = 'in', res = 600, pointsize = 4)
plot(cumsum(rho))
abline(h=0.7, col="red", lty = 2)
dev.off()


## ##########################################
## 8.- Plot 2 first components
## ##########################################



# R default plot
fname <- paste0("plots/chr8_PCA_1and2_U.png")
png(fname, width = 8, height = 8, units = 'in', res = 600, pointsize = 4)
plot(u[,"PC2"], u[,"PC1"])
dev.off()


# ggplot plot
fname <- paste0("plots/chr8_PCA_1and2_ggplot_U.png")
png(fname, width = 8, height = 8, units = 'in', res = 600, pointsize = 4)
p <- ggplot(as.data.frame(u), aes(x=PC2, y=PC1))+
  geom_point()
print(p)
dev.off()


# ggplot plot
fname <- paste0("plots/chr8_PCA_1and2_ggplot_V.png")
png(fname, width = 8, height = 8, units = 'in', res = 600, pointsize = 4)
p <- ggplot(as.data.frame(v), aes(x=PC2, y=PC1))+
  geom_point()
print(p)
dev.off()



# ggplot plot
fname <- paste0("plots/chr8_PCA_1and2_ggplot_Ui.png")
png(fname, width = 8, height = 8, units = 'in', res = 600, pointsize = 4)
p <- ggplot(as.data.frame(u), aes(x=PC1, y=PC2))+
  geom_point()
print(p)
dev.off()





# Read data again
library(rhdf5)
#..# setwd("/home/isglobal.lan/dpelegri/working/TFM/results") #.. SERVER 06 ..#
#..# filehdf <-  "data/chr8.h5"
#..# h5fsvd = H5Fopen(filehdf)
#..# u <- h5fsvd$SVD$genofilter_maf$u
#..# u <- as.data.frame(u)
#..# colnames(u) <- sprintf("PC%s",seq(1:dim(u)[2]))
#..# h5closeAll()


# NEW !! -- > ggplot plot2



# CODI PER A EXECUCIÓ GRÀFIC -> SERVER06
fname <- paste0("plots/chr8_PCA.pdf")
#fname <- paste0("plots/PCA/chr8_PCA.png")


#png(fname, width = 10, height = 8, units = 'in', res = 900, pointsize = 4)
#..# pdf(fname, res = 600)
rasterpdf::raster_pdf(fname, res = 600, width = 12, height = 7)
p <- ggplot(as.data.frame(u), aes(PC2, PC1, color = PC2) )+
  geom_point(  size = 0.5 ) +
  scale_color_gradientn(colours = c('yellow','orange','darkred'),
                        n.breaks = 3,
                        breaks = c(-0.002,-0.001,0.001),
                        labels = c("NI/NI", "NI/II", "II/II") ) +
  labs(x = "PC2",
       y = "PC1") +
  theme(
    text = element_text( size = 15),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    plot.background = element_blank(),
    # axis.line = element_blank(),
    panel.grid = element_blank(),
    # plot.title = element_text(size = 12, hjust = 0.5),
    legend.title=element_blank()) +
  geom_rangeframe()#  +
# ylim(-0.0018,-0.008)

print(p)
dev.off()





# # NOTA :: CODI PER A EXECUCIÓ -> SERVER06
# #       Al servidor 5 aquesta gràfica no funciona bé.
# fname <- paste0("plots/chr8_PCA.pdf")
# # png(fname, width = 10, height = 8, units = 'in', res = 900, pointsize = 4)
# pdf.options(reset = TRUE, onefile = FALSE)
# pdf(fname)
# p <- ggplot(as.data.frame(u), aes(PC2, PC1, color = PC2) )+
#   geom_point(  size = 1 ) +
#   geom_rangeframe() +
#   scale_color_gradientn(colours = c('yellow','orange','darkred'),
#                         n.breaks = 3,
#                         breaks = c(-0.01,0.03,0.08),
#                         labels = c("NI/NI", "NI/II", "II/II") ) +
#   labs(title = "Chr8 - Inverse SNPs",
#        x = "PC2",
#        y = "PC1") +
#   theme(legend.background = element_blank(),
#         legend.key = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_blank(),
#         strip.background = element_blank(),
#         plot.background = element_blank(),
#         axis.line = element_blank(),
#         panel.grid = element_blank(),
#         plot.title = element_text(size = 18, hjust = 0.5),
#         plot.caption = element_text(size = 20),
#         legend.title=element_blank()) +
#   ylim(-0.037,-0.08)
#
# print(p)
# dev.off()






#
#
# # # ####### COMPROVACIONS BÀSIQUES #####
# # # geno <- h5fsvd$invs$geno
# # #
# # #
# setwd("/home/isglobal.lan/dpelegri/working/TFM/results") #.. SERVER 06 ..#
# filename <- "/home/isglobal.lan/dpelegri/working/TFM/data/UkbBiobank/chr8inv.gds"
# filehdf <-  "data/chr8.h5"
# #
# # Open gds file
#
# # gds<-gdsfmt::openfn.gds(filename)
# # igenn <- index.gdsn(gds, "genotype")
# # genn<- read.gdsn( igenn)
# # genn[1:14,1:13]
#
# #
# # h5ls(filehdf)
# # h5fsvd = H5Fopen(filehdf)
# # geno <- h5fsvd$invs$genofilter
# # geno[488370:488377,1:13]
# #
# # h5closeAll()
