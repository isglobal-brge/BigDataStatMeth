#include "include/looe.h"
#include "include/optimizedproduct.h"





// Eigen::MatrixXd rcppinversecpp_eig ( double& lambda, const Eigen::Map<Eigen::VectorXd>* Lambda, const Eigen::Map<Eigen::MatrixXd>* Q )
void rcppinversecpp_eig ( double& lambda, const Eigen::Map<Eigen::VectorXd>* Lambda, 
                          const Eigen::Map<Eigen::MatrixXd>* Q, Eigen::Map<Eigen::MatrixXd>* Ginv  )
{

  Eigen::VectorXd wd(Lambda->size());
  int diff=0;
  
  if( Lambda->isZero() || Q->isZero())
    Rcpp::stop("SVD results should be provided. \n");

  
  if (lambda!=1) 
    diff=lambda;
  
  if (lambda == 1){
    wd = 1/(Lambda->array());
  }else {
    wd = 1/(Lambda->array() + lambda);
  }
  

  Eigen::Map<Eigen::VectorXd> w(wd.data(), wd.size());
  
  rcpp_xwxt_eig( Q, &w, Ginv );
}


Rcpp::NumericMatrix rcppinversecpp ( double& lambda, const Rcpp::NumericVector& Lambda, const Rcpp::NumericMatrix& Q )
{

  if( Lambda.isNULL()|| Q.isNULL())
    Rcpp::stop("SVD results should be provided. \n");
  
  Rcpp::NumericVector W;
  
  if (lambda == 1){
    W = 1/(Lambda);
  }else {
    W = 1/(Lambda + lambda);
  }

  // Rcpp::NumericMatrix Ginv = rcpp_xwxt( Q, W); 
  
  return(rcpp_xwxt( Q, W));
  
}



double rcpplooei_eig( double& lambda, const Eigen::Map<Eigen::VectorXd>* Lambda, const Eigen::Map<Eigen::MatrixXd>* Q,const Eigen::Map<Eigen::VectorXd>* Y)
{

  Eigen::MatrixXd Ginvd = Eigen::MatrixXd::Zero(Q->rows(), Q->rows());
  Eigen::Map<Eigen::MatrixXd> Ginv (Ginvd.data(), Q->rows(), Q->rows());
  
  rcppinversecpp_eig( lambda, Lambda, Q, &Ginv);
  
  Eigen::VectorXd cted = Eigen::VectorXd::Zero(Ginv.rows());
  Eigen::Map<Eigen::VectorXd> cte (cted.data(), Ginv.rows());
  
  //..// Eigen::VectorXd cte = 
  rcpp_parallel_Xy_eigen(&Ginv, Y, &cte);
  
  //..// Rcpp::Rcout<<"\n"<<cte<<"\n";
  
  Eigen::VectorXd ans = cte.binaryExpr(Ginv.diagonal(), cteGinvDiagonal<double>());
  return (ans.sum());
}

double rcpplooei( double& lambda, const Rcpp::NumericVector& Lambda, const Rcpp::NumericMatrix& Q,const Rcpp::NumericVector& Y)
{
  Rcpp::NumericMatrix Ginv = rcppinversecpp( lambda,  Lambda, Q);
  Rcpp::NumericVector cte = rcpp_parallel_Xy(Ginv, Y);
  Rcpp::NumericVector ans = cte/Rcpp::diag(Ginv);
  //ans = parallelpow2(ans);
  return (parallelVectorSum(parallelpow2(ans)));
  
}




// [[Rcpp::export]]
Rcpp::RObject LOOE(Rcpp::RObject& X, Rcpp::RObject& Y, Rcpp::Nullable<double> nl = R_NilValue,
                   Rcpp::Nullable<double> ml  = R_NilValue,
                   Rcpp::Nullable<Rcpp::RObject> l  = R_NilValue)
{
  
  // Variable declaration
  double nlambdas, maxlambda;
  Rcpp::NumericMatrix eX;
  Rcpp::NumericVector eY, lambdas;

  // Variable initialization
  if(nl.isNotNull())  nlambdas = Rcpp::as<double> (nl);
  else    nlambdas = 100;

  if(nl.isNotNull())    maxlambda = Rcpp::as<double> (ml);
  else    maxlambda = 1;
  
  if(l.isNull())  lambdas = generate_seq<Rcpp::NumericVector>( 0.01, maxlambda, (maxlambda/nlambdas));
  else  lambdas = Rcpp::NumericVector(l);
  
  if ( X.isS4() == true)    eX = read_DelayedArray_rcpp(Rcpp::NumericMatrix(X));
  else{
    try{  eX = Rcpp::NumericMatrix(X);   }
    catch(std::exception &ex) { }
  }

  if ( Y.isS4() == true)   
  {
    if(Rf_isVector(Y)) 
      throw std::range_error("Object must be a vector not S4 type");
  }  else {
    eY = Rcpp::NumericVector(Y);
  } 

  // eX = rcpp_parallel_tCrossProd(eX);
  eX = rcpp_parallel_tCrossProd(eX);
  
  svd ee = RcppBDsvd(eX, int(), int(), false);
  
  Rcpp::NumericVector Lambda = ee.d;
  Lambda.push_back(0);
  replace_zero(&Lambda);
  
  Rcpp::NumericMatrix Q = Rcpp::cbind(ee.v,Rcpp::NumericMatrix(ee.v.rows(),1));
  Rcpp::NumericVector looei = Rcpp::NumericVector(lambdas.size());
  
  Rcpp::NumericVector W;
  
  for(size_t i=0; i<lambdas.size(); i = i+1)
  {
    // looei[i] = rcpplooei(lambdas[i], Lambda, Q, eY);
    Rcpp::NumericVector W;
    double ds=0;
    
    if(lambdas[i]!=1) ds = lambdas[i];
    W = 1/(Lambda + ds);
    
    Rcpp::NumericMatrix Ginv = rcpp_xwxt(Q, W);
    Rcpp::NumericVector cte = rcpp_parallel_Xy(Ginv, eY);
    Rcpp::NumericVector ans = cte/Rcpp::diag(Ginv);
    looei[i] = parallelVectorSum(parallelpow2(ans));

  }

  Rcpp::NumericVector::iterator it = std::min_element(looei.begin(), looei.end());

  Rcpp::NumericMatrix Ginv= rcppinversecpp( lambdas[it - looei.begin()], Lambda, Q);
  Rcpp::NumericMatrix coef = rcpp_parallel_XYProd( Ginv, Rcpp::NumericMatrix(X));
  
  
  coef =  rcpp_parallel_Xy( Rcpp::transpose(coef), eY);
  
  return Rcpp::List::create(Rcpp::Named("coef") = coef,
                            Rcpp::Named("Ginv") = Ginv,
                            Rcpp::Named("lambda.min") = lambdas[it - looei.begin()],
                            Rcpp::Named("lambdas") = lambdas,
                            Rcpp::Named("looe") = looei);

}





// [[Rcpp::export]]
Rcpp::RObject LOOE_Eigen(Rcpp::RObject& X, Rcpp::RObject& Y, Rcpp::Nullable<double> nl = R_NilValue,
                   Rcpp::Nullable<double> ml  = R_NilValue,
                   Rcpp::Nullable<Rcpp::RObject> l  = R_NilValue)
{
  
  // Variable declaration
  double nlambdas, maxlambda;
  // Eigen::Map<Eigen::MatrixXd> eX = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  // Eigen::Map<Eigen::MatrixXd> eX = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  //..// Eigen::MatrixXd eX = Rcpp::as<Eigen::MatrixXd>(X);
  Eigen::MatrixXd eX(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X));
  
  Eigen::Map<Eigen::VectorXd> eY = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(Y);
  Eigen::VectorXd  lambdas;

  
  // Variable initialization
  if(nl.isNotNull())  nlambdas = Rcpp::as<double> (nl);
  else    nlambdas = 100;
  
  if(nl.isNotNull())    maxlambda = Rcpp::as<double> (ml);
  else    maxlambda = 1;
  

  if(l.isNull())  lambdas = generate_seq<Eigen::VectorXd>( 0.01, maxlambda, (maxlambda/nlambdas));
  else  lambdas = Rcpp::as<Eigen::VectorXd>(l);
  
      
  
//..//   if ( X.isS4() == true)    eX = read_DelayedArray_rcpp (Rcpp::NumericMatrix(X));
//..//   else{
//..//    try{  eX = Eigen::MatrixXd(X);   }
//..//        catch(std::exception &ex) { }
//..//   }
  
//..//  if ( Y.isS4() == true)   
//..//  {
//..//    if(Rf_isVector(Y)) 
//..//      throw std::range_error("Object must be a vector not S4 type");
//..//  }  else {
//..//    eY = Rcpp::NumericVector(Y);
//..//  } 


  eX = bdcrossproduct(eX);
  
  svdeig ee = RcppBDsvd_eig(eX, int(), int(), false);
  


  Eigen::VectorXd Lambdad(ee.d.size()+1);
  Lambdad<<ee.d,0;

  Eigen::Map<Eigen::VectorXd> Lambda(Lambdad.data(), Lambdad.size());
  setZeros(&Lambda);
  

  Eigen::MatrixXd Qd(ee.v.rows(), ee.v.cols()+1) ;
  Qd << ee.v, Eigen::VectorXd::Zero(Qd.rows());
  Eigen::Map<Eigen::MatrixXd> Q(Qd.data(), Qd.rows(), Qd.cols());
  
  
  Eigen::VectorXd looei (lambdas.size());
  double lambdamin = 0; 

  
  for(size_t i=0; i<lambdas.size(); i = i+1)
  {
     looei[i] = rcpplooei_eig(lambdas[i], &Lambda, &Q, &eY);
  }
  
  lambdamin = lambdas[index_val(looei.minCoeff(), looei)];
  
  // Eigen::MatrixXd Ginv = rcppinversecpp_eig(lambdamin, &Lambda, &Q);
  Eigen::MatrixXd Ginvd = Eigen::MatrixXd::Zero(Q.rows(), Q.rows());
  Eigen::Map<Eigen::MatrixXd> Ginv (Ginvd.data(), Q.rows(), Q.rows());
  
  rcppinversecpp_eig(lambdamin, &Lambda, &Q, &Ginv);
  
  
  Eigen::MatrixXd mcoef = rcpp_parallel_XYProd_eigen(Ginv, Rcpp::as<Eigen::MatrixXd >(X)).transpose();
  Eigen::Map<Eigen::MatrixXd> pmcoef (mcoef.data(), mcoef.rows(), mcoef.cols());
  
  Eigen::VectorXd coefd = Eigen::VectorXd::Zero(mcoef.cols());;
  Eigen::Map<Eigen::VectorXd> coef (coefd.data(), mcoef.cols());
  
  rcpp_parallel_Xy_eigen( &pmcoef, &eY, &coef);
  
  return Rcpp::List::create(Rcpp::Named("coef") = coef,
                            Rcpp::Named("Ginv") = Ginv,
                            Rcpp::Named("lambda.min") = lambdamin,
                            Rcpp::Named("lambdas") = lambdas,
                            Rcpp::Named("looe") = looei);
  
// return (Rcpp::wrap(Q));
}






























































////////////////////////////////////////////////////////////////////////////////////////////////////





// Eigen::MatrixXd rcppinversecpp_eig ( double& lambda, const Eigen::Map<Eigen::VectorXd>* Lambda, const Eigen::Map<Eigen::MatrixXd>* Q )
Eigen::MatrixXd rcppinversecpp_eig_blast ( double lambda, const Eigen::VectorXd& Lambda, 
                                           const Eigen::MatrixXd& Q, bool paral )
{
  
  Eigen::VectorXd w(Lambda.size());
  Eigen::MatrixXd Ginv = Eigen::MatrixXd::Zero(Q.rows(), Q.rows());
  
  if( Lambda.isZero() || Q.isZero())
    Rcpp::stop("SVD results should be provided. \n");
  
  if (lambda == 1){
    w.array() = 1/(Lambda.array());
  }else {
    w.array() = 1/(Lambda.array() + lambda);
  }

  
  if(paral==true)  {
    // Ginv =  block_matrix_mul_parallel(block_matrix_mul_parallel( Q, w.asDiagonal(),256), Q.transpose(), 256);
    Ginv = Q*w.asDiagonal()*Q.transpose();
  } else {
    Ginv = Q*w.asDiagonal()*Q.transpose();
    // Eigen::MatrixXd Ginv = xwxt(Q,w) ; // Opció força mes lenta. Perquè???
  };
  
  return(Ginv);
}



double rcpplooei_eig_blast( double lambda, const Eigen::VectorXd& Lambda, 
                            const Eigen::MatrixXd& Q, const Eigen::VectorXd& Y, 
                            bool paral)
{

  Eigen::MatrixXd mY = Y.col(0);
  Eigen::MatrixXd Ginv = rcppinversecpp_eig_blast( lambda, Lambda, Q, paral);
  
  if(paral==true)  {
    Eigen::MatrixXd cte = block_matrix_mul_parallel(Ginv,Y.col(0),128);
    Eigen::VectorXd dGinv = Ginv.diagonal();
    Eigen::VectorXd rans = (cte.array() / dGinv.array()).array().pow(2);
    
    return (rans.array().sum());
  }else  {
    Eigen::VectorXd cte = Ginv * Y;
    Eigen::VectorXd dGinv = Ginv.diagonal();
    Eigen::VectorXd rans = (cte.array() / dGinv.array()).array().pow(2);
    
    return (rans.array().sum());
  }
  /*
  Eigen::VectorXd dGinv = Ginv.diagonal();
  Eigen::VectorXd rans = (cte.array() / dGinv.array()).array().pow(2);
  
  return (rans.array().sum());
   */
  
}






// [[Rcpp::export]]
Rcpp::RObject LOOE_BLAST(Rcpp::RObject& X, Rcpp::RObject& Y, bool paral, 
                         Rcpp::Nullable<double> nl = R_NilValue,
                         Rcpp::Nullable<double> ml  = R_NilValue,
                         Rcpp::Nullable<Rcpp::RObject> l  = R_NilValue)
{
  try
  {
      // Variable declaration
      double nlambdas, maxlambda;
      Eigen::MatrixXd rX;
      Eigen::VectorXd eY;
      Eigen::VectorXd  lambdas;
      
      
      // Variable initialization
      if(nl.isNotNull())  nlambdas = Rcpp::as<double> (nl);
      else    nlambdas = 100;
      
      if(nl.isNotNull())    maxlambda = Rcpp::as<double> (ml);
      else    maxlambda = 1;
      
      if(l.isNull())  lambdas = generate_seq<Eigen::VectorXd>( 0.01, maxlambda, (maxlambda/nlambdas));
      else  lambdas = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(l);

      // Read DelayedArray's X and Y     
      if ( X.isS4() == true)    
      {
        rX = read_DelayedArray(X);
      } else {
        try{  
          // rX = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
          rX = Rcpp::as<Eigen::MatrixXd>(X);
        }
        catch(std::exception &ex) { }
      }

      if ( Y.isS4() == true)   
      {
        eY = read_DelayedArray(Y);
      }  else {
        // eY = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(Y);
        eY = Rcpp::as<Eigen::VectorXd>(Y);
      } 

      Eigen::MatrixXd eX = bdtcrossproduct(rX);
      
      svdeig ee = RcppBDsvd_eig(eX, int(), int(), false);

      Eigen::VectorXd Lambda(ee.d.size()+1);
      Lambda<<ee.d,0;
      Lambda = equalZero(Lambda);
      
      
      Eigen::MatrixXd Q(ee.v.rows(), ee.v.cols()+1) ;
      Q << ee.v, Eigen::VectorXd::Zero(Q.rows());
      Eigen::VectorXd looei (lambdas.size());
      double lambdamin = 0; 
      
      for(size_t i=0; i<lambdas.size(); i = i+1)
      {
        looei[i] = rcpplooei_eig_blast(lambdas[i], Lambda, Q, eY, paral);
      }
      
      lambdamin = lambdas[index_val(looei.minCoeff(), looei)];

      Eigen::MatrixXd Ginv = rcppinversecpp_eig_blast(lambdamin, Lambda, Q, paral);
      
      Eigen::MatrixXd coef = (Ginv * rX).adjoint() * eY;
      
      return Rcpp::List::create(Rcpp::Named("coef") = coef,
                                Rcpp::Named("Ginv") = Ginv,
                                Rcpp::Named("lambda.min") = lambdamin,
                                Rcpp::Named("lambdas") = lambdas,
                                Rcpp::Named("looe") = looei);
      
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  
}









/*** R


library(microbenchmark)
library(DelayedArray)



n <- 300
p <- 100
M <- matrix(rnorm(n*p), nrow=n, ncol=p)

Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]

MD <- DelayedArray(M)
YD <- DelayedArray(as.matrix(Y))

solprodparal3 <- BigDataStatMeth::LOOE_BLAST(MD,YD, paral = TRUE)
solprodparal3 <- BigDataStatMeth::LOOE_BLAST(M,Y, paral = TRUE)
solprodparal3$coef[1:10]
solprodparal3$looe[1:10]
solprodparal3$lambda.min
solprodparal <- BigDataStatMeth::LOOE_Eigen(M,Y)

temps <- microbenchmark(# solprodparal <- BigDataStatMeth::LOOE(M,Y),
                        solprodparal2 <- BigDataStatMeth::LOOE_BLAST(M,Y),
                                       times = 1L); temps
solprodparal2$coef[1:10]
solprodparal2$looe[1:10]

RcppParallel::setThreadOptions(numThreads = 30)

temps <- microbenchmark(#solparal1 <- BigDataStatMeth::LOOE(M,Y),
                        #solprodparal2 <- BigDataStatMeth::LOOE_Eigen(M,Y),
                        solrcppparal <- BigDataStatMeth::LOOE_BLAST(M,Y,paral=TRUE),
                        solrcpp <- BigDataStatMeth::LOOE_BLAST(M,Y,paral=FALSE),
                        sol.orig <- LOOE_orig(tcrossprod(M),Y),
                        ans.orig <- solveEigen_orig(M, as.matrix(Y), lambda=sol.orig$lambda.min),
                        times = 3L); temps
;

LOOE.all <- function(X,Y)
{
  
}

solrcppparal$coef[1:10]
solrcpp$coef[1:10]
sol.orig[1:10]



solprodparal$looe[1:10]




temps <- microbenchmark(solprodparal <- BigDataStatMeth::LOOE(M,Y),
                        sol.orig <- LOOE.all(tcrossprod(M),Y),
                        times = 1L)
temps
solprodparal$coef[1:10]
solprodparal$lambda.min


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

solprodparal3$coef[1:10]
solprodparal3$lambda.min



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


aa <- BDsvd(M, 99,100,FALSE)
aa$`d$`[1:10]

dd <- BDsvd(M, 99, 100, TRUE)
dd$`d$`[1:10]^2

bb <- svd(tcrossprod(M))
bb$d[1:10]
bb$v

cc <- eigen(tcrossprod(M))
cc$values[1:10]

gg <- runSVD(tcrossprod(M), k=500)
gg$d[1:10]


RcppParallel::setThreadOptions(numThreads = 40)
setThreadOptions(numThreads = defaultNumThreads() / 2)


  res <- benchmark(solprodparal <- BigDataStatMeth::LOOE(M,Y),
                   sol3 <- LOOE_orig(tcrossprod(M),Y),
                   ans3 <- solveEigen_orig(M, as.matrix(Y), lambda=sol3$lambda.min),
      order="relative", replications = c(2))
  res[,1:4] 
  
  temps <- microbenchmark(solprodparal <- BigDataStatMeth::LOOE(M,Y),
                   sol3 <- LOOE_orig(tcrossprod(M),Y),
                   ans3 <- solveEigen_orig(M, as.matrix(Y), lambda=sol3$lambda.min),
                   times = 2L);temps


  
*/
