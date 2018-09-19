X <- matrix(rnorm(1000), 10, 100)
Y <- crosspord(X)

system.time(ans1 <-  ff(Y))
system.time(ans2 <- solve(Y))

identical(ans1, ans2)
max(ans1 - ans2)


library(Rcpp)
library(DelayedArray)
library("microbenchmark")

# PROVES MATRIU ENTERS
val <- matrix(rep(1:200, 1:200), ncol = 100)
da_mat <- DelayedArray(seed = val) ; da_mat

results <- microbenchmark(transR <-  t(da_mat),
                          transCppEigen <-  BDtrans_numeric(da_mat),
                          transCpp <-  BDtrans_hdf5(da_mat),
                          times = 50L)

print(summary(results)[, c(1:7)],digits=1)

# PROVES MATRIU REALS
provesdelay <- matrix(runif(100000), ncol=200, nrow=500)
X <- DelayedArray(blah)

results <- microbenchmark(transR <- t(X),
                          transCppEigen <- BDtrans_numeric(X),
                          transCpp <- BDtrans_hdf5(X),
                          times = 50L)
print(summary(results)[, c(1:7)],digits=1)

# PROVES MATRIU CARACTERS
delaystring <- matrix(sample(LETTERS, 10000, TRUE), ncol=200, nrow=50)
X <- DelayedArray(delaystring)


results <- microbenchmark(transR <- t(X),
                          #transCppEigen <- BDtrans_numeric(X),
                          transCpp <- BDtrans_hdf5(X),
                          times = 50L)
print(summary(results)[, c(1:7)],digits=1)