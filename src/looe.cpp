#include "include/looe.h"
#include "include/optimizedproduct.h"


// Eigen::MatrixXd rcppinversecpp_eig ( double& lambda, const Eigen::Map<Eigen::VectorXd>* Lambda, const Eigen::Map<Eigen::MatrixXd>* Q )
Eigen::MatrixXd rcppinversecpp ( double lambda, const Eigen::VectorXd& Lambda, 
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
    Eigen::MatrixXd tmp = Xwd_parallel(Q,w);
    Ginv  = block_matrix_mul_parallel( tmp, Q.adjoint(), 128);
    // Ginv = Q*w.asDiagonal()*Q.transpose();
  } else {
    Eigen::MatrixXd tmp = Xwd(Q,w);
    Ginv  = block_matrix_mul( tmp, Q.adjoint(), 128);
    // Ginv = Q*w.asDiagonal()*Q.transpose();
    //..// Ginv = Q*w.asDiagonal()*Q.transpose();
  };
  
  return(Ginv);
}



double rcpplooei( double lambda, const Eigen::VectorXd& Lambda, 
                            const Eigen::MatrixXd& Q, const Eigen::VectorXd& Y, 
                            bool paral)
{

  Eigen::MatrixXd mY = Y.col(0);
  Eigen::MatrixXd Ginv = rcppinversecpp( lambda, Lambda, Q, paral);
  
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




//' LOOE
//' 
//' This function performs ...
//' 
//' @param X ......
//' @param Y ...
//' @param paral ...
//' @param nl ...
//' @param ml ...
//' @param l ...
//' @return coef ...
//' @return Ginv  ...
//' @return lambda.min  ...
//' @return lambdas  ...
//' @return looe ...
//' @examples
//' n <- 500
//' p <- 200
//' M <- matrix(rnorm(n*p), nrow=n, ncol=p)
//' Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]
//' 
//' MD <- DelayedArray(M)
//' YD <- DelayedArray(as.matrix(Y))
//' 
//' looe <- LOOE(M,Y,paral=TRUE),
//' looe_delayed <- LOOE(MD,YD,paral=TRUE)
//' 
//' looe$lambda.min
//' looe$lambdas
//' looe_delayed$coef[1:10]
//' @export
// [[Rcpp::export]]
Rcpp::RObject LOOE(Rcpp::RObject& X, Rcpp::RObject& Y, bool paral, 
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
      int tid, chunk=1;
      
      
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
      
      //..//svdeig ee = RcppBDsvd_eig(eX, int(), int(), false);
      svdeig ee = RcppBDsvd(eX, int(), int(), false);

      Eigen::VectorXd Lambda(ee.d.size()+1);
      Lambda<<ee.d,0;
      Lambda = equalZero(Lambda);
      
      
      Eigen::MatrixXd Q(ee.v.rows(), ee.v.cols()+1) ;
      Q << ee.v, Eigen::VectorXd::Zero(Q.rows());
      Eigen::VectorXd looei (lambdas.size());
      double lambdamin = 0; 
      
      /*
      for(size_t i=0; i<lambdas.size(); i = i+1)
      {
        looei[i] = rcpplooei(lambdas[i], Lambda, Q, eY, paral);
      }
       */
      
      #pragma omp parallel shared(looei,lambdas, Lambda, Q, eY, paral) // private (tid)
      {
        // tid = omp_get_thread_num();
        //només per fer proves dels threads i saber que està paralelitzant, sinó no cal tenir-ho descomentat
        /*
        if (tid == 0) {
        Rcpp::Rcout << "Number of threads: " << omp_get_num_threads() << "\n";
        }
        */
        // #pragma omp for schedule (static, chunk)
        #pragma omp for schedule (dynamic)
        for(size_t i=0; i<lambdas.size(); i = i+1)
        {
          looei[i] = rcpplooei(lambdas[i], Lambda, Q, eY, paral);
        }
      }

      lambdamin = lambdas[index_val(looei.minCoeff(), looei)];

      Eigen::MatrixXd Ginv = rcppinversecpp(lambdamin, Lambda, Q, paral);
      
      //Eigen::MatrixXd coef = (Ginv * rX).adjoint() * eY;
       Eigen::MatrixXd coef = block_matrix_mul_parallel(Ginv,rX,128).adjoint() * eY;
      
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
library(BigDataStatMeth)



n <- 500
p <- 200
M <- matrix(rnorm(n*p), nrow=n, ncol=p)

Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]

MD <- DelayedArray(M)
YD <- DelayedArray(as.matrix(Y))

#looe1 <- LOOE(M,Y,paral=TRUE)
res <- microbenchmark( looe1 <- LOOE(M,Y,paral=TRUE),
                       looe3 <- LOOE(MD,YD,paral=TRUE),
                        times = 2L, unit = "s")

print(summary(res)[, c(1:7)],digits=3)

looe1$coef[1:10]
looe1$lambda.min




n <- 1800
p <- 700
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
Y <- 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5]
MD <- DelayedArray(M)
YD <- DelayedArray(as.matrix(Y))

res <- microbenchmark( looe1 <- LOOE(M,Y,paral=TRUE),
                       # looe2 <- LOOE(M,Y,paral=FALSE),
                       looe3 <- LOOE(MD,YD,paral=TRUE),
                       # looe4 <- LOOE(MD,YD,paral=FALSE),
                       # looe5 <- LOOE.all(M,Y),
                       ans1 <- glmnet::cv.glmnet(M, Y),
                       times = 2L, unit = "s")

print(summary(res)[, c(1:7)],digits=3)

looe1$coef[1:10]
# looe2$coef[1:10]
looe3$coef[1:10]
# looe4$coef[1:10]
looe5[1:10]
# el resultado es:
mm <- which(ans1$lambda==ans1$lambda.min)
ans1$glmnet.fit$beta[1:10,mm]






  
*/
