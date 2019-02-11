#include "include/svdDecomposition.h"
#include "include/ReadDelayedData.h"



svd RcppBDsvd ( Eigen::MatrixXd X, int k, int nev, bool normalize )
{

  svd retsvd;
  if( k==0 )    k = (std::min(X.rows(), X.cols()))-1;
  if(nev == 0)  nev = k + 1 ;
  
  Eigen::MatrixXd Xtcp;
  if(normalize ==true )  {
    Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(RcppNormalize_Data(X))));
  }else {
    Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(X)));
  }
  
  Spectra::DenseSymMatProd<double> op(Xtcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, nev);
  
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  
  if(eigs.info() == Spectra::SUCCESSFUL)
  {
    retsvd.d = eigs.eigenvalues().cwiseSqrt();
    retsvd.u = eigs.eigenvectors();
  }
  
  Eigen::MatrixXd Xcp;
  if(normalize ==true )  {
    Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(RcppNormalize_Data(X))));  
  }else {
    Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(X)));
  }
  
  
  Spectra::DenseSymMatProd<double> opv(Xcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigsv(&opv, k, nev);
  
  // Initialize and compute
  eigsv.init();
  nconv = eigsv.compute();

  // Retrieve results
  if(eigsv.info() == Spectra::SUCCESSFUL)
  {
    retsvd.v = eigsv.eigenvectors();
  }
  
  return retsvd;
  
}



// [[Rcpp::export]]
Rcpp::RObject BDsvd (const Rcpp::RObject & x, int k, int nev, bool normalize )
{

  auto dmtype = beachmat::find_sexp_type(x);

  size_t ncols = 0, nrows=0;
  Eigen::MatrixXd X;
  Rcpp::List ret;
  
  if ( dmtype == INTSXP || dmtype==REALSXP ) {
    if ( x.isS4() == true){
      X = read_DelayedArray(x);
    }else {
      try{
        X = Rcpp::as<Eigen::MatrixXd >(x);
      }catch(std::exception &ex) {
        X = Rcpp::as<Eigen::VectorXd >(x);
      }
    }
    
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
  
  Eigen::MatrixXd Xtcp;

  if( normalize == true)  {
    Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(RcppNormalize_Data(X))));  
  } else {
    Xtcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_tCrossProd( Rcpp::wrap(X)));  
  }
  
  
  Spectra::DenseSymMatProd<double> op(Xtcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, k, nev);
  
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  

  if(eigs.info() == Spectra::SUCCESSFUL)
  {
    ret["d$"] = eigs.eigenvalues();
    ret["u$"] = eigs.eigenvectors();
  }
  
  Eigen::MatrixXd Xcp;
  
  if(normalize==true )  {
    //Rcpp::Rcout<<"Normalitzem\n";
    Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(RcppNormalize_Data(X))));
  } else
  {
    //Rcpp::Rcout<<"NO Normalitzem\n";
    Xcp =  Rcpp::as<Eigen::MatrixXd> (rcpp_parallel_CrossProd( Rcpp::wrap(X)));
  }
  //. OPCIÓ BONA.// Eigen::MatrixXd Xcp = CrossProduct(Normalize_Data(Rcpp::wrap(X)),false); // Normalitzem matriu
  // Eigen::MatrixXd Xcp = CrossProduct_eig(X,false); // No normalitzem matriu - matriu normalitzada prèviament
  

  
  Spectra::DenseSymMatProd<double> opv(Xcp);
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigsv(&opv, k, nev);
  
  // Initialize and compute
  eigsv.init();
  nconv = eigsv.compute();
  

  // Retrieve results
  if(eigsv.info() == Spectra::SUCCESSFUL)
  {
    ret["v$"] = eigsv.eigenvectors();
  }

  return Rcpp::wrap(ret);
  
}