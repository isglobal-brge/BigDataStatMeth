#include "include/pca.h"



Eigen::VectorXd cumsum(Eigen::VectorXd x){
  // initialize an accumulator variable
  double acc = 0;
  // initialize the result vector
  Eigen::VectorXd res = Eigen::VectorXd::Zero(x.size());
  for(int i = 0; i < x.size(); i++){
    acc += x[i];
    res[i] = acc;
  }
  return res;
}

/*

// 
//' PCA Descomposition
//' 
//' Compute PCA
//' 
//' @param matrix or DelayedArray 
//' @return 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdPCA(const Rcpp::RObject & x, int k)
{
  try
  {
    
    Eigen::MatrixXd X;
    Eigen::MatrixXd U, V, Lambda;
    Eigen::VectorXd lambda;
    Eigen::MatrixXd C, D, varcoord;
    bool bkeepon = true;
    // svdeig sing;
    
    // Rcpp::Rcout<<"\nPrintem el que anem fent : \n";
    // Read DelayedArray's X and Y     
    if ( x.isS4() == true)    
    {
      X = read_DelayedArray(x);
    } else {
      try{  
        // rX = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X);
        X = Rcpp::as<Eigen::MatrixXd>(x);
      }
      catch(std::exception &ex) { }
    }
    
    Rcpp::Rcout<<"\n Lectura dades OK : \n";
    X = RcppNormalize_Data(X);
    X = X*(1/sqrt(X.rows()-1));  // Tipify
    
    // Rcpp::Rcout<<"\n Tipificació dades OK : \n";
    
    {
      Rcpp::Rcout<<"\n Calculem crossproduct....: \n";
      Eigen::MatrixXd covX = bdcrossproduct(X);
      Rcpp::Rcout<<"\n Crossproduct OK : \n";
      
      Rcpp::Rcout<<"\n Calculem SVD....: \n";
      svdeig sing = RcppbdSVD(covX, k, int(), false);  
      Rcpp::Rcout<<"\n SVD OK : \n";
      
      if(sing.bokuv == true && sing.bokd==true)
      {
        U = sing.u;
        V = sing.v;
        lambda = sing.d;
      }else {
        bkeepon = false;
        throw(Rcpp::exception("Error with svd decomposition","pca.cpp",4));
      }
    }
    
    if( bkeepon == true)
    {
      Rcpp::Rcout<<"\n Descomposició +  assignació  OK : \n";
      
      {
        Eigen::VectorXd prov = lambda.array().sqrt();
        D = prov.asDiagonal();
      }
      
      Rcpp::Rcout<<"\n D assignada  OK : \n";
      
      //..// Lambda = lambda.asDiagonal();
      
      
      // Cumulative variance
      //..// Eigen::VectorXd pvac = cumsum(lambda/lambda.array().sum()); // percentatge de vari`ancia acumulada pels components pvac
      
      // Get correlations
      {
        Rcpp::Rcout<<"\n Calculem produc per a C....: \n";
        Eigen::MatrixXd mprov = block_matrix_mul_parallel( Eigen::MatrixXd::Identity(V.rows(), V.cols()), V, 512 );
        Rcpp::Rcout<<"\n Product OK : \n";
        C = block_matrix_mul_parallel( mprov, D, 512 );
      }
      
      Rcpp::Rcout<<"\n Correlacions realitzades -  OK : \n";
      
      // Get coordinates
      {
        svdeig singX = RcppbdSVD(X, k, int(), false);
        if(singX.bokuv == true && singX.bokd==true)
        {
          varcoord = block_matrix_mul_parallel(X.adjoint(), singX.u, 512 );  
        } else {
          throw(Rcpp::exception("Error with svd decomposition","pca.cpp",4));
        }
      }
      Rcpp::Rcout<<"\n Càlcul coordenades realitzades -  OK : \n";
      
      Lambda = lambda.asDiagonal(); 
    }

    
    
    return Rcpp::List::create(Rcpp::Named("U") = U,
                              Rcpp::Named("V") = V,
                              Rcpp::Named("D") = D,
                              Rcpp::Named("Lambda") = Lambda,
                              Rcpp::Named("lambda") = lambda,
                              Rcpp::Named("pvac") = cumsum(lambda/lambda.array().sum()), // percentatge de vari`ancia acumulada pels components pvac,
                              Rcpp::Named("var.contr") = V.pow(2), // Variable contrib.
                              Rcpp::Named("C") = C,
                              Rcpp::Named("var.coord") = varcoord,
                              Rcpp::Named("var.cos2") = C.pow(2)
    );
    
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  
}


*/
/*** R

*/
