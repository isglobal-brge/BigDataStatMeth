#include "include/looe.h"


Eigen::MatrixXd rcppinversecpp (Eigen::MatrixXd X, double lambda, bool eigen, Eigen::VectorXd Lambda, Eigen::MatrixXd Q )
{
  /* TOT AIXÒ SERIA INNECESSARI?? JA TENIM LA DESCOMPOSICIÓ FETA ANTERIORMENT I ESTALVIEM FEINA I TEMPS !!! 
  svd ee;
 
  if(eigen)
  {
    Rcpp::Rcout<<"Ara haurem de fer la descomposició...\n";
    if( X.isZero() == false )
    {
      Rcpp::Rcout<<"Sembla ser que no és 0.......\n";
      // Eigen::MatrixXd eX = Rcpp::as<Eigen::MatrixXd > (rcpp_parallel_tCrossProd(Rcpp::wrap(X)));
      ee = RcppBDsvd(X, int(), int(), false);
      Rcpp::Rcout<<"Ha arribat a fer la descomposició\n";
    } else {
      Rcpp::stop("X matrix should be provided. \n");
    }
      
    Lambdas.resize(ee.d.size());
    
    for (size_t i=0; i< ee.d.size(); i++)
    {
      if( ee.d[i]<0 || Rcpp::NumericVector::is_na(ee.d[i]) ) Lambda[i] = 0;
      else Lambda[i] = ee.d[i];
    }
      
    Eigen::MatrixXd Q = ee.v;
      
  } else {*/
    
  if( Lambda.isZero() || Q.isZero())
    Rcpp::stop("SVD results should be provided. \n");

  Eigen::VectorXd W;
  if (lambda == 1){
    W = 1/(Lambda.array());
  }else {
    W = 1/(Lambda.array() + lambda);
  }

  Eigen::MatrixXd Ginv = Rcpp::as<Eigen::MatrixXd > (rcpp_xwxt( Rcpp::wrap(Q), Rcpp::wrap(W))); 

  return(Ginv);

}




double rcpplooei( double lambda, Eigen::VectorXd Lambda, Eigen::MatrixXd Q, Eigen::VectorXd Y)
{

  Eigen::VectorXd ans;
  Eigen::MatrixXd nullmat;

  // Rcpp::Rcout<<"Inici rcppinversecpp\n";
  Eigen::MatrixXd Ginv = rcppinversecpp(nullmat, lambda, false, Lambda, Q);
  // Rcpp::Rcout<<"Fi rcppinversecpp - Inici cte\n";
  Eigen::VectorXd cte = Ginv * Y;   // Canviar-ho per una operació en paral·lel??
  // Rcpp::Rcout<<"Fi cte ->  Inici Ginvcte\n";
  ans = cte.binaryExpr(Ginv.diagonal(), cteGinvDiagonal<double>());
  // Rcpp::Rcout<<"Fi Ginvcte\n";
  return (ans.sum());
    
}


// [[Rcpp::export]]
Rcpp::RObject LOOE(Rcpp::RObject X, Rcpp::RObject Y, Rcpp::Nullable<double> nl = R_NilValue,
                   Rcpp::Nullable<double> ml  = R_NilValue,
                   Rcpp::Nullable<Rcpp::RObject> l  = R_NilValue)
{
  
  // Variable declaration
  double nlambdas, maxlambda;
  Eigen::MatrixXd eX;
  Eigen::VectorXd lambdas, eY;
  

  // Variable initialization
  if(nl.isNotNull())  nlambdas = Rcpp::as<double> (nl);
  else    nlambdas = 100;

  if(nl.isNotNull())    maxlambda = Rcpp::as<double> (ml);
  else    maxlambda = 1;
  
  if(l.isNull())  lambdas = generate_seq( 0.01, maxlambda, (maxlambda/nlambdas));
  else  lambdas = Rcpp::as<Eigen::VectorXd> (l);

  if ( X.isS4() == true)    eX = read_DelayedArray(X);
  else{
    try{  eX = Rcpp::as<Eigen::MatrixXd >(X);  }
    catch(std::exception &ex) {  eX = Rcpp::as<Eigen::VectorXd >(X);   }
  }

  if ( Y.isS4() == true)    eY = read_DelayedArray(Y);
  else  eY = Rcpp::as<Eigen::VectorXd >(Y);

  eX = Rcpp::as<Eigen::MatrixXd > (rcpp_parallel_tCrossProd(Rcpp::wrap(eX)));

  svd ee = RcppBDsvd(eX, int(), int(), false);

  Eigen::VectorXd Lambda = ee.d;
  // Length K = k+1
  Lambda.conservativeResize(ee.d.size()+1);
  Lambda[Lambda.size()] = 0;
  Lambda = Lambda.unaryExpr(islambdazero<double>());

  Eigen::MatrixXd Q = ee.v;
  // add column ; columns  K = k+1
  Q.conservativeResize(Q.rows(), Q.cols()+1);
  Q.col(Q.cols()-1) = Eigen::VectorXd::Zero(Q.rows());
  
  Eigen::VectorXd looei = Eigen::VectorXd::Zero(lambdas.size());
  double lambdamin = 0, looemin = 99999999;
  
  for(size_t i=0; i<lambdas.size(); i = i+1)
  {
    looei[i] = rcpplooei(lambdas[i], Lambda, Q, eY);
    
    if(looei[i]<looemin)  {
      looemin = looei[i];
      lambdamin = lambdas[i];
    }
  }


  
  Eigen::MatrixXd Ginv = rcppinversecpp(eX, lambdamin, false, Lambda, Q);
  Rcpp::NumericMatrix coef = rcpp_parallel_XYProd( Rcpp::wrap(Ginv), Rcpp::wrap(X));

  coef =  rcpp_parallel_Xy( Rcpp::transpose(coef), Rcpp::wrap(eY));

  return Rcpp::List::create(Rcpp::Named("coef") = coef,
                            Rcpp::Named("Ginv") = Ginv,
                            Rcpp::Named("lambda.min") = lambdamin,
                            Rcpp::Named("lambdas") = lambdas,
                            Rcpp::Named("looe") = looei);

}





/*** R

library(rfunctions)
library(biocSingular)
library(rbenchmark)

BiocManager::install("BiocSingular", version = "3.8")

n <- 25
p <- 25
M <- matrix(rnorm(n*p), nrow=n, ncol=p)

  Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]

sol <- LOOE_orig(M,Y)
ans2 <- solveEigen(M, Y, lambda=sol$lambda.min)
ans2[1:10]



### PROVES ACTUALS
n <- 300
p <- 100
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
#  DM <- DelayedArray(M)
Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]


res <- benchmark(solprodparal <- Prodparal::LOOE(M,Y),
                 sol3 <- LOOE_orig(tcrossprod(M),Y),
                 ans3 <- solveEigen_orig(M, as.matrix(Y), lambda=sol3$lambda.min),
                 order="relative", replications = c(3))
res[,1:4] 

solprodparal <- Prodparal::LOOE(M,Y)
solprodparal$coef[1:10]
solprodparal$lambda.min

sol3 <- LOOE_orig(tcrossprod(M),Y)
ans3 <- solveEigen_orig(M, as.matrix(Y), lambda=sol3$lambda.min)
sol3$lambda.min
ans3[1:10]


RcppParallel::setThreadOptions(numThreads = 40)

res <- benchmark(solprodparal <- Prodparal::LOOE(M,Y),
                 sol3 <- LOOE_orig(tcrossprod(M),Y),
                 ans3 <- solveEigen_orig(M, as.matrix(Y), lambda=sol3$lambda.min),
                 order="relative", replications = c(1))
res[,1:4] 

###


aa <- BDsvd(M, 10,100,FALSE)
aa$`d$`[1:10]

dd <- BDsvd(M, 99, 100, TRUE)
dd$`d$`[1:10]^2

bb <- svd(tcrossprod(M))
bb$d[1:10]

cc <- eigen(tcrossprod(M))
cc$values[1:10]

gg <- runSVD(tcrossprod(M), k=500)
gg$d[1:10]


RcppParallel::setThreadOptions(numThreads = 40)
setThreadOptions(numThreads = defaultNumThreads() / 2)
  res <- benchmark(solprodparal <- Prodparal::LOOE(M,Y),
                   sol3 <- LOOE_orig(tcrossprod(M),Y),
                   ans3 <- solveEigen_orig(M, as.matrix(Y), lambda=sol3$lambda.min),
      order="relative", replications = c(1))
  res[,1:4] 


  
*/
