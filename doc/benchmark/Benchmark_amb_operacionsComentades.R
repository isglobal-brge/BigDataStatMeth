library(microbenchmark)
library(parallel)
library(ggplot2)
library(DelayedArray)
library(BigDataStatMeth)



setwd("~/Library/Mobile Documents/com~apple~CloudDocs/UAB/Curs 2018-2019/PFG/Proves/BigDataStatMeth_1_20190227.DD")

# --    PRODUCTE DE MATRIUS x BLOCS   --
# --------------------------------------

# Creem data.frame per emmagatzemar els resultats
result.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                         max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),"A(n)"=numeric(),
                         "A(m)"=numeric(), "B(m)"=numeric(), ncores = numeric(), nexec=numeric())
repet <- 1 # Repeticions microbenchmark
tipusop <- c('par','seq') # tipus execució paral·lel, sequencial
tipusdades <- c('R type','Delayed') # tipus execució paral·lel, sequencial
ncores <- detectCores() # Cores de la màquina

# Fem una execució per diferents tamanys de matrius, i per diferents tamanys de blocs
# per saber com actuen els algoritmes en les diferents situacions
matrixsize <- c(1000,2500)
blocksize <- c(2, 4, 8, 16, 32, 64, 128, 256, 512)


for (k in 1:length(matrixsize))
{
  n <- matrixsize[k]
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  B <- matrix(rnorm(n*n), nrow=n, ncol=n)
  
  for ( i in seq(2, 6, by=4)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
  {
     bparal <- TRUE
     for (j in 1:length(tipusop))
     {
       if(j==1) 
         bparal <- TRUE
       else
         bparal <- FALSE
       for(l in 1:length(tipusdades))
       {
         if(tipusdades=="Delayed")
           res <- microbenchmark(MPBD <- blockmult(AD,BD,i, bparal), times = repet, unit = "s")
          else
            res <- microbenchmark(MPB <- blockmult(A,B,i, bparal), times = repet, unit = "s")
             
         #res <- microbenchmark(MPB <- blockmult(A,B,i, bparal),
         #                       MPBD <- blockmult(AD,BD,i, bparal), 
         #                       times = repet, unit = "s")
         #res <- microbenchmark(MPBP <- blockmult(A,B,i, TRUE),
         #                      MPB <- blockmult(A,B,i, FALSE),
         #                      MPBP <- blockmult(DA,DB,i, TRUE),
         #                      MPB <- blockmult(DA,DB,i, FALSE),
         #                      times = repet, unit = "s")
         resdata <- as.data.frame(summary(res)[, c(1:7)])
         resdata <- cbind(resdata,tipusop[j],tipusdades[l],i,dim(A)[1],dim(A)[2],dim(B)[2],ncores, repet)
         result.df <- rbind(result.df,resdata)
       }
      
      }
    
  }
}

# Guardem les dades per reutilitzar-les si convé.
   write.csv(result.df,"./doc/benchmark/multi5blocksize.csv")

# Llegim els resultats de l'i7 del fitxer
#   setwd("~/Library/Mobile Documents/com~apple~CloudDocs/UAB/Curs 2018-2019/PFG/Proves/ProvesEigen")
#   result.df <- read.csv("./benchmark/multi7blocksize.csv")
# Readaptem el nom de les columnes del data.frame per tenir-los controlats -- OPCIÓ QUAN LLEGIM DEL FITXER, SINÓ VA SENSE EL 'id'
#   colnames(result.df) <-  c('id','expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
#                              'tipusop', 'bloc_size','M', 'K', 'N', 'ncores', 'nexec')


   ######
# ANOTACIONS IMPORTANTS microbenchmark : 
#   ## Change default unit to relative runtime
#   options(microbenchmark.unit="relative")
#   print(res)
#   ## Change default unit to evaluations per second
#   options(microbenchmark.unit="eps")
#   print(res)
   ######

colnames(result.df) <-  c('expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
                          'tipusop','tipusdades', 'bloc_size','M', 'K', 'N', 'ncores', 'nexec')

# Grafiquem els resultats
p <- ggplot(result.df, aes( x = result.df$bloc_size,
                            y = result.df$mean, 
                            group= interaction(tipusop,tipusdades,M))) + 
  geom_line(aes(color=interaction(tipusop,tipusdades,M))) +
  xlab('NomBlocs') +
  ylab('temps')

print(p)


# Busquem el nombre de blocs on temps execució és mínim per poder fer la propera prova amb 
# tamany matrius

## mins <- aggregate( result.df$mean ~ tipusop, data = result.df, FUN = function(x) min(x))
## minims <- result.df[which(round(result.df$mean,8)==round(mins[,2],8)),]



resultsize.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                             max=numeric(),tipusop=character(),bloc_size=numeric(),"A(n)"=numeric(), "A(m)"=numeric(), "B(m)"=numeric(),
                             ncores = numeric(), nexec=numeric())

repet <- 2


for ( i in seq(1000, 8000, by=1000)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
{
  n <- i
  bloc_size <- 128
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  B <- matrix(rnorm(n*n), nrow=n, ncol=n)
  
#  for (j in 1:length(tipusop))
#  {
#    if(j==1) 
#    {
#      bparal <- TRUE
#      #bloc_size <- minims[which(tipusop=='par'),2]
#      bloc_size <- 128
#    } else
#    {
#      bparal <- FALSE
#      #bloc_size <- minims[which(tipusop=='seq'),2]
#      bloc_size <- 128
#    }
    
    #res <- microbenchmark(CPP <- blockmult(A,B,blocs[i], bparal), times = repet)
    res <- microbenchmark(MPBP <- blockmult( A, B, bloc_size, TRUE), 
                          MPBS <- blockmult( A, B, bloc_size, FALSE), 
                          times = repet, unit = "s")
    resdata <- as.data.frame(summary(res)[, c(1:7)])
    resdata <- cbind(resdata, tipusop[j], bloc_size, dim(A)[1], dim(A)[2], dim(B)[2], ncores, repet)
    resultsize.df <- rbind(resultsize.df,resdata)
  # }
  
}


#   # Guardem les dades per reutilitzar-les si convé.
     write.csv(resultsize.df,"./doc/benchmark/multBLOCKi5matrixsize.csv")

#   Llegim les dades de l'i7 referent als temps per calcular diferents tamanys de matrius
#   amb una unitat fixada de blocs, per a totes el mateix.
#     resultsize.df <- read.csv("./benchmark/multi7matrixsize.csv")
#     colnames(resultsize.df) <-  c('id','expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
#                                   'tipusop', 'bloc_size','M', 'K', 'N', 'ncores', 'nexec')


# Readaptem el nom de les columnes del data.frame per tenir-los controlats
colnames(resultsize.df) <-  c('expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
                              'tipusop', 'bloc_size','M', 'K', 'N', 'ncores', 'nexec')


# Grafiquem els resultats
p <- ggplot(resultsize.df, aes( x = resultsize.df$M, 
                                y = resultsize.df$mean, 
                                group= interaction(tipusop))) + 
  geom_line(aes(color=interaction(tipusop))) +
  xlab('NomBlocs') +
  ylab('temps')

print(p)


# --    PRODUCTE DE MATRIUS - NO BLOCS   --
# ------------------------------------------


resultsize.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                             max=numeric(),tipusop=character(),bloc_size=numeric(),"A(n)"=numeric(), "A(m)"=numeric(), "B(m)"=numeric(),
                             ncores = numeric(), nexec=numeric())

repet <- 2

for ( i in seq(1000, 3000, by=1000)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
{
  n <- i
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  B <- matrix(rnorm(n*n), nrow=n, ncol=n)
  DA <- DelayedArray(A)
  DB <- DelayedArray(B)
  
  res <- microbenchmark(parXYProd(A,B),
                        parXYProd(DA,DB),
                        parXYProdBlock(A,B),
                        parXYProdBlock(DA,DB), 
                        A%*%B,
                        times = repet, unit = "s")
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop[j], bloc_size, dim(A)[1], dim(A)[2], dim(B)[2], ncores, repet)
  resultsize.df <- rbind(resultsize.df,resdata)
  # }
  
}

write.csv(resultsize.df,"./doc/benchmark/multi5matrixsize.csv")

# Readaptem el nom de les columnes del data.frame per tenir-los controlats
colnames(resultsize.df) <-  c('expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
                              'tipusop', 'bloc_size','M', 'K', 'N', 'ncores', 'nexec')


# Grafiquem els resultats
p <- ggplot(resultsize.df, aes( x = resultsize.df$M, 
                                y = resultsize.df$mean, 
                                group= interaction(tipusop))) + 
  geom_line(aes(color=interaction(tipusop))) +
  xlab('NomBlocs') +
  ylab('temps')

print(p)

