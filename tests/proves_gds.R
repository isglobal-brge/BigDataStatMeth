library(SNPRelate)
snpgdsBED2GDS("colorectal.bed", "colorectal.fam", 
              "colorectal.bim", "test.gds")

genofile <- snpgdsOpen("test.gds")

snpgdsSummary(genofile)

read.gdsn(index.gdsn(genofile, "genotype"), start=c(1000,1000), 
          count=c(15,5))


closefn.gds(genofile)


# de la funció 'read.gdsn'
.Call("gdsObjReadData", node=index.gdsn(genofile, "genotype"), 
      start=c(1,1), count=c(5,5), simplify="none", .useraw=FALSE, 
      list(.value=NULL, .substitute=NULL))