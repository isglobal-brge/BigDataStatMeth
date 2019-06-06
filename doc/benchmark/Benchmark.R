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
                         max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(), M=numeric(),
                         K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())
repet <- 1 # Repeticions microbenchmark
tipusop <- c('par','seq') # tipus execució paral·lel, sequencial
tipusdades <- c('R type','Delayed') # tipus execució paral·lel, sequencial
ncores <- detectCores() # Cores de la màquina

# Fem una execució per diferents tamanys de matrius, i per diferents tamanys de blocs
# per saber com actuen els algoritmes en les diferents situacions
matrixsize <- c(1000,2500)
#..# blocksize <- c(2, 4, 8, 16, 32, 64, 128, 256, 512)


for (k in 1:length(matrixsize))
{
  n <- matrixsize[k]
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  B <- matrix(rnorm(n*n), nrow=n, ncol=n)
  AD <- DelayedArray(A)
  BD <- DelayedArray(B)
  
  for ( i in seq(2, 160, by=8)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
  {
     
     for (j in 1:length(tipusop))
     {
       if(j==1) 
         bparal <- TRUE
       else
         bparal <- FALSE
       for(l in 1:length(tipusdades))
       {
         if(tipusdades[l]=="Delayed")
           res <- microbenchmark(MPBD <- blockmult(AD,BD,i, bparal), times = repet, unit = "s")
          else
            res <- microbenchmark(MPB <- blockmult(A,B,i, bparal), times = repet, unit = "s")

         resdata <- as.data.frame(summary(res)[, c(1:7)])
         resdata <- cbind(resdata, tipusop[j], 
                          tipusdades[l],i,M=dim(A)[1], K=dim(A)[2], N=dim(B)[2], ncores, repet)
         result.df <- rbind(result.df,resdata)
       }
      
      }
    
  }
}

# Guardem les dades per reutilitzar-les si convé.
   write.csv(result.df,"./doc/benchmark/multi5blocksize.csv")

# Llegim els resultats de l'i7 del fitxer
#   setwd("~/Library/Mobile Documents/com~apple~CloudDocs/UAB/Curs 2018-2019/PFG/Proves/ProvesEigen")
   result.df <- read.csv("./doc/benchmark/multi5blocksize.csv")
# Readaptem el nom de les columnes del data.frame per tenir-los controlats -- OPCIÓ QUAN LLEGIM DEL FITXER, SINÓ VA SENSE EL 'id'
   colnames(result.df) <-  c('id','expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
                              'tipusop','tipusdades', 'bloc_size','M', 'K', 'N', 'ncores', 'nexec')

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
  xlab('mida bloc') +
  ylab('temps(s)')

print(p)

result.par <- result.df[which(result.df[,11]>2),]
p <- ggplot(result.par, aes( x = result.par$bloc_size,
                            y = result.par$mean, 
                            group= interaction(tipusop,tipusdades,M))) + 
  geom_line(aes(color=interaction(tipusop,tipusdades,M))) +
  xlab('mida bloc') +
  ylab('temps(s)')

print(p)



# Busquem el nombre de blocs on temps execució és mínim per poder fer la propera prova amb 
# tamany matrius

## mins <- aggregate( result.df$mean ~ tipusop, data = result.df, FUN = function(x) min(x))
## minims <- result.df[which(round(result.df$mean,8)==round(mins[,2],8)),]




result.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                         max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(), M=numeric(),
                         K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

repet <- 2

for ( i in seq(1000, 5000, by=1000)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
{
  n <- i
  bloc_size <- 128
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  B <- matrix(rnorm(n*n), nrow=n, ncol=n)
  AD <- DelayedArray(A)
  BD <- DelayedArray(B)

  for (j in 1:length(tipusop))
  {
    if(j==1) 
      bparal <- TRUE
    else
      bparal <- FALSE
    for(l in 1:length(tipusdades))
    {
      if(tipusdades[l]=="Delayed")
        res <- microbenchmark(MPBD <- blockmult(AD,BD,bloc_size, bparal), times = repet, unit = "s")
      else
        res <- microbenchmark(MPB <- blockmult(A,B,bloc_size, bparal), times = repet, unit = "s")
      
      resdata <- as.data.frame(summary(res)[, c(1:7)])
      resdata <- cbind(resdata, tipusop[j], tipusdades[l], i, M=dim(A)[1], K=dim(A)[2], N=dim(B)[2], ncores, repet)
      result.df <- rbind(result.df,resdata)
    }
    
  }

}

write.csv(result.df,"./doc/benchmark/multBLOCKi5matrixsize.csv")
colnames(resultsize.df) <-  c('expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
                              'tipusop', 'tipusdades','bloc_size','M', 'K', 'N', 'ncores', 'nexec')

# Grafiquem els resultats
p <- ggplot(resultsize.df, aes( x = resultsize.df$M, 
                                y = resultsize.df$mean, 
                                group= interaction(tipusop,tipusdades))) + 
  geom_line(aes(color=interaction(tipusop,tipusdades))) +
  xlab('mida matriu (nxn)') +
  ylab('temps (s)')

print(p)



# --    PRODUCTE DE MATRIUS - NO BLOCS   --
# ------------------------------------------


resultsize.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                         max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                         K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

repet <- 2

for ( i in seq(1000, 2500, by=500)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
{
  n <- i
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  B <- matrix(rnorm(n*n), nrow=n, ncol=n)
  DA <- DelayedArray(A)
  DB <- DelayedArray(B)
  
  for(l in 1:length(tipusdades))
  {
    if(tipusdades[l]=="Delayed")
      res <- microbenchmark(parXYProd(DA,DB),
                            blockmult(DA,DB,block_size = n, bparal), # Forcem un únic bloc per tota la matriu equiv a no blocs
                            times = repet, unit = "s")
    else
      res <- microbenchmark(parXYProd(A,B), 
                            blockmult(A,B,block_size = n, bparal), # Forcem un únic bloc per tota la matriu equiv a no blocs
                            times = repet, unit = "s")
    
    resdata <- as.data.frame(summary(res)[, c(1:7)])
    resdata <- cbind(resdata, tipusop = 'par',tipusdades = tipusdades[l], 0, M=dim(A)[1], K=dim(A)[2], N=dim(B)[2], ncores, repet)
    resultsize.df <- rbind(resultsize.df,resdata)
  }
  
  res <- microbenchmark(A%*%B, times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = 'seq', tipusdades = 'R type', 0, M = dim(A)[1], K = dim(A)[2], N = dim(B)[2], ncores, repet)
  resultsize.df <- rbind(resultsize.df,resdata)

}

write.csv(resultsize.df,"./doc/benchmark/multi5matrixsize.csv")

# Readaptem el nom de les columnes del data.frame per tenir-los controlats
colnames(resultsize.df) <-  c('expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
                              'tipusop','tipusdades', 'bloc_size','M', 'K', 'N', 'ncores', 'nexec')


# Grafiquem els resultats
p <- ggplot(resultsize.df, aes( x = resultsize.df$M, 
                                y = resultsize.df$mean, 
                                group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)



# --    PRODUCTE DE MATRIUS - BLOCS vs NO BLOCS (PARAL·LEL)   --
# --------------------------------------------------------------


resultbvsnb.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                             max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                             K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

repet <- 2
bparal <- TRUE

for ( i in seq(1000, 3000, by=1000)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
{
  n <- i
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  B <- matrix(rnorm(n*n), nrow=n, ncol=n)
  
  for(j in seq(1,2,by=1) )
  {
    if (j==1) 
    {
      blocksize = 128
      res <- microbenchmark(blockmult(A,B,blocksize, bparal), times = repet, unit = "s")
    }
    else {
      blocksize = 0
      res <- microbenchmark(parXYProd(A,B),
                            blockmult(A,B,block_size = n, bparal),
                            times = repet, unit = "s")
    }
  
    resdata <- as.data.frame(summary(res)[, c(1:7)])
    resdata <- cbind(resdata, tipusop = 'par',tipusdades = 'R type', bloc_size = blocksize, 
                     M=dim(A)[1], K=dim(A)[2], N=dim(B)[2], ncores, repet)
    resultbvsnb.df <- rbind(resultbvsnb.df,resdata)
  }
}


write.csv(resultbvsnb.df,"./doc/benchmark/multi5blocvsnobloc.csv")

#..# # Readaptem el nom de les columnes del data.frame per tenir-los controlats
#..# colnames(resultbvsnb.df) <-  c('expr', 'min', 'lq', 'mean', 'median', 'uq', 'max', 
#..#                               'tipusop','tipusdades', 'bloc_size','M', 'K', 'N', 'ncores', 'nexec')


# Grafiquem els resultats
p <- ggplot(resultbvsnb.df, aes( x = resultbvsnb.df$M, 
                                y = resultbvsnb.df$mean, 
                                group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)



# --    CROSSPROD i TCROSSPROD    --
# ----------------------------------

# crossprod

results.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                          max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                          K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

repet <- 4

for ( i in seq(500, 3000, by=500)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
{
  n <- i
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  AD <- DelayedArray(A)
  
  
  res <- microbenchmark(bdcrossprod(A, transposed = FALSE), # crossprod
                        crossprod(A),
                        t(A)%*%A,
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'R type', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=0, ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  res <- microbenchmark(bdcrossprod(AD, transposed = FALSE), # crossprod
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'Delayed', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=0, ncores, repet)
  results.df <- rbind(results.df,resdata)
  
}

write.csv(results.df,"./doc/benchmark/crossprodi5.csv")

# Grafiquem els resultats
p <- ggplot(results.df, aes( x = results.df$M, 
                                 y = results.df$mean, 
                                 group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)



#   tcrossprod

results.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                          max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                          K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())


for ( i in seq(500, 3000, by=500)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
{
  n <- i
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  AD <- DelayedArray(A)
  
  
  res <- microbenchmark(bdcrossprod(A, transposed = TRUE), # crossprod
                        tcrossprod(A),
                        A%*%t(A),
                        times = repet, unit = "ms")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'R type', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=0, ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  res <- microbenchmark(bdcrossprod(AD, transposed = TRUE), # crossprod
                        times = repet, unit = "ms")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'Delayed', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=0, ncores, repet)
  results.df <- rbind(results.df,resdata)
  
}

write.csv(results.df,"./doc/benchmark/tcrossprodi5.csv")

# Grafiquem els resultats
p <- ggplot(results.df, aes( x = results.df$M, 
                             y = results.df$mean, 
                             group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)




# --    CROSSPROD i TCROSSPROD AMB PESOS    --
# --------------------------------------------

# crossprod amb pesos
results.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                          max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                          K=numeric(), w=numeric(), ncores = numeric(), nexec=numeric())

repet <- 10

for ( i in seq(500, 2500, by=500))  
{
  n <- i
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  u <- runif(n)
  w <- u * (1 - u)
  AD <- DelayedArray(A)
  wD <- DelayedArray(as.matrix(w))
  
  res <- microbenchmark(bdwproduct(A, w,"xtwx"),
                        crossprod(A, w*A),
                        t(A)%*%diag(w)%*%A,
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'R type', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], w=length(w), ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  res <- microbenchmark(bdwproduct(AD, w,"xtwx"),
                        bdwproduct(AD, wD,"xtwx"),
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'Delayed', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], w=length(w), ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  
}


write.csv(results.df,"./doc/benchmark/crossweightprodi5.csv")

# Grafiquem els resultats
p <- ggplot(results.df, aes( x = results.df$M, 
                             y = results.df$mean, 
                             group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)



# tcrossprod amb pesos
results.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                          max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                          K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

repet <- 5

for ( i in seq(500, 3000, by=500)) ## Apliquem blocs quadrats 2x2,4x4,6x6.... fins  
{
  n <- i
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  u <- runif(n)
  w <- u * (1 - u)
  AD <- DelayedArray(A)
  wD <- DelayedArray(as.matrix(w))
  
  res <- microbenchmark(bdwproduct(A, w,"xwxt"),
                        A%*%diag(w)%*%t(A),
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'R type', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=length(w), ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  res <- microbenchmark(bdwproduct(AD, w,"xwxt"),
                        bdwproduct(AD, wD,"xwxt"),
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'Delayed', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=length(w), ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  
}


write.csv(results.df,"./doc/benchmark/tcrossweightprodi5.csv")

# Grafiquem els resultats
p <- ggplot(results.df, aes( x = results.df$M, 
                             y = results.df$mean, 
                             group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)


# --    INVERSA CHOLESKY    --
# --------------------------------------------

# Genera matrius "definides positives"
Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}



results.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                          max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                          K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

repet <- 4

for ( i in seq(500, 2000, by=500)) 
{
  
  A <- Posdef(n=i, ev=1:i)
  AD <- DelayedArray(A)

  res <- microbenchmark(inversechol_par(A),
                        bdInvCholesky_LDL_eigen(A),
                        solve(A),
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'R type', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=0, ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  res <- microbenchmark(inversechol_par(AD),
                        bdInvCholesky_LDL_eigen(AD),
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'Delayed', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=0, ncores, repet)
  results.df <- rbind(results.df,resdata)
  
}


write.csv(results.df,"./doc/benchmark/invcholi5.csv")

# Grafiquem els resultats
p <- ggplot(results.df, aes( x = results.df$M, 
                             y = results.df$mean, 
                             group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)



# --    DESCOMPOSICIÓN SVD    --
# ------------------------------



results.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                          max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                          K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

repet <- 2

for ( i in seq(300, 1500, by=300)) 
{
  n <- i
  A <- matrix(rnorm(n*n), nrow=n, ncol=n)
  AD <- DelayedArray(A)
  
  res <- microbenchmark( bdSVD( A, n-1, n, FALSE), # No normalitza la matriu
                        svd(tcrossprod(A)),
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'R type', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=0, ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  res <- microbenchmark(bdSVD(AD, n-1,n,FALSE),
                        times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'Delayed', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=0, ncores, repet)
  results.df <- rbind(results.df,resdata)
  
}


write.csv(results.df,"./doc/benchmark/svdi5.csv")

# Grafiquem els resultats
p <- ggplot(results.df, aes( x = results.df$M, 
                             y = results.df$mean, 
                             group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)




# --    LOOE    --
# ------------------------------




# FUNCIONS R ORIGINALS : (Afegida una funció que crida les 2 funcions necessaries per obtenir els coefficients per tal de 
#                         poder obtenir el temps total.)

inversecpp_orig <- function(X, lambda=1, eigen=TRUE,
                            Lambda, Q){
  if (eigen){
    ee <- eigen(X, symmetric = TRUE)  # mirar la librería BiocSingular
    Lambda <- ee$values
    Lambda[Lambda<0] <- 0
    Q <- ee$vectors
  }
  else
    if(missing(Lambda) | missing(Q))
      stop('SVD results should be provided. \n')
  
  if (lambda == 1)
    W <-  1/Lambda
  else
    W <- 1/(Lambda + lambda)
  
  Ginv <- rfunctions::crossprodcpp(t(Q), W)  # implementar xwxt  (en la librería rfunctions está xtwx) 
  # Q%*%diag(W)%*%t(Q)
  Ginv
}

solveEigen_orig <- function(X, Y, lambda){   # X DelayedArray (HDF5), Y un vector
  XX <- tcrossprod(X)
  Ginv <- inversecpp_orig(XX, lambda=lambda)  # Con DelayedArray (HDF5)
  coef <- t(Ginv%*%X)%*%Y
  coef
}


LOOE.i_orig <- function(lambda, Lambda, Q, Y){
  Ginv <- inversecpp_orig(lambda=lambda, Lambda=Lambda,
                          Q=Q, eigen=FALSE)
  print(dim(Ginv));
  print(dim(Y));
  cte <- Ginv%*%Y
  ans <- sum((cte/diag(Ginv))^2)
  ans
}



#
# Compute LOOE
#

LOOE_orig <- function(X, Y, nlambdas=100, max.lambda=1, lambdas){
  
  if (missing(lambdas)){
    lambdas <- seq(0.01, max.lambda, length=nlambdas)
  }
  ee <- eigen(X, symmetric = TRUE)
  Lambda <- ee$values
  
  Lambda[Lambda<0] <- 0
  Q <- ee$vectors
  looe <- sapply(lambdas, LOOE.i_orig, Lambda=Lambda, Q=Q, Y=Y)
  
  lambda.min <- lambdas[which.min(looe)]
  
  Ginv <- inversecpp_orig(X, lambda.min)
  
  ans <- list(looe=looe, Ginv=Ginv, lambdas=lambdas,
              lambda.min=lambda.min)
  ans
}

LOOE.all <- function(X, Y){
  sol <- LOOE_orig(tcrossprod(X),Y)
  return(solveEigen_orig(X, as.matrix(Y), lambda=sol$lambda.min))
}

### FI FUNCIONS ORIGINALS R ###




results.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                          max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                          K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

repet <- 2

for ( i in seq(300, 1500, by=300)) 
{
  if(i>=900) p <- 500
  else p <- 100
  n <- i
  A <- matrix(rnorm(n*p), nrow=n, ncol=p)
  
  Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]
  AD <- DelayedArray(M)
  YD <- DelayedArray(as.matrix(Y))
  
  res <- microbenchmark( LOOE_BLAST(A,Y,paral=TRUE),
                         LOOE_BLAST(A,Y,paral=FALSE),
                         LOOE.all(A,Y),
                         times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'R type', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=length(Y), ncores, repet)
  results.df <- rbind(results.df,resdata)
  
  res <- microbenchmark( LOOE_BLAST(AD,YD,paral=TRUE),
                         LOOE_BLAST(AD,YD,paral=FALSE),
                         times = repet, unit = "s")
  
  resdata <- as.data.frame(summary(res)[, c(1:7)])
  resdata <- cbind(resdata, tipusop = '',tipusdades = 'Delayed', bloc_size = 0, 
                   M=dim(A)[1], K=dim(A)[2], N=length(Y), ncores, repet)
  results.df <- rbind(results.df,resdata)
  
}


write.csv(results.df,"./doc/benchmark/looei5.csv")

# Grafiquem els resultats
p <- ggplot(results.df, aes( x = results.df$M, 
                             y = results.df$mean, 
                             group= interaction(expr))) + 
  geom_line(aes(color=interaction(expr))) +
  xlab('mida matriu(nxn)') +
  ylab('temps (s)')

print(p)



