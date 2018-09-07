X <- matrix(rnorm(1000), 10, 100)
Y <- crosspord(X)

system.time(ans1 <-  ff(Y))
system.time(ans2 <- solve(Y))

identical(ans1, ans2)
max(ans1 - ans2)