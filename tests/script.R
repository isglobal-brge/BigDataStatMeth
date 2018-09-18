X <- matrix(rnorm(1000), 10, 100)
Y <- crosspord(X)

system.time(ans1 <-  ff(Y))
system.time(ans2 <- solve(Y))

identical(ans1, ans2)
max(ans1 - ans2)



library(DelayedArray)
library("microbenchmark")

blah <- matrix(runif(10000), ncol=200, nrow=50)
X <- DelayedArray(blah)

# 
results <- microbenchmark(transR = t(X),
                          transCpp = BDtrans_numeric(X),
                          transCpp = BDtrans_hdf5(X),
                          times = 50L)
