# La idea es resolver un modelo de regresión Lasso ....  esto se puede hacer con la función `glmnet` del paquete homónimo
# Generamos datos
# con 3 coeficientes distintos de 0 (el 1, el 2 y el 5)

n <- 500
p <- 100
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
MM <- tcrossprod(M)
Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]

# Estimamos los betas de las 100 variables con `glmnet`:

ans1 <- glmnet::cv.glmnet(M, Y)

# el resultado es:
mm <- which(ans1$lambda==ans1$lambda.min)

ans1$glmnet.fit$beta[1:10,mm]

# V1         V2         V3         V4         V5         V6         V7         V8         V9        V10 
# 2.3423592  1.5403282  0.0000000  0.0000000 -0.3444285  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000

# ahora veamos nuestra implementación:

sol <- LOOE(tcrossprod(M),Y)
ans2 <- solveEigen(M, Y, lambda=sol$lambda.min)
ans2[1:10]

# [1]  2.399937e+00  1.599961e+00 -1.200198e-06  6.537065e-08  -3.999911e-01 -3.781633e-06  5.707634e-07 -1.428479e-06
# [9]  1.132573e-06 -4.996911e-06

# con la función de `rfunctions` obtendríamos lo mismo (como es obvio):

ans3 <- rfunctions::cgls(M,Y, lambda=sol$lambda.min)
ans3$x[1:10]

# [1]  2.399937e+00  1.599960e+00 -4.407301e-06  6.936192e-06  -3.999929e-01 -3.114494e-06  1.153641e-06 -5.488184e-06
# [9] -2.234088e-06 -3.076770e-06```