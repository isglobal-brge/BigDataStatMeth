#include "include/utilities.h"


struct Zeros : public RcppParallel::Worker  {
  
  // output matrix to write to

  // input matrix
  Eigen::Map<Eigen::VectorXd> v;
  
  // other variables
  std::size_t numcol;
  
  // Constructor
  Zeros( Eigen::Map<Eigen::VectorXd>* v)
    : v(*v) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    // size_t ncmat = mat.ncol();
    for (std::size_t i = begin; i < end; i++) 
    {
      if(v[i]<0 || Rcpp::NumericVector::is_na(v[i])==true) v[i]=0;
    }
  }
};


void setZeros(Eigen::Map<Eigen::VectorXd>* lambda) {
  
  
  try
  {
     try{
        
       Zeros zeros(lambda);
       RcppParallel::parallelFor(0, lambda->size(), zeros);
        
     }catch(std::exception &ex) {	
        forward_exception_to_r(ex);
     }
      
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }  
 
}



Eigen::VectorXd equalZero(Eigen::VectorXd v) 
{
  
  try
  {
    for(size_t i=0; i<v.size(); i++)
    {
      if(v[i]<0 || Rcpp::NumericVector::is_na(v[i])==true) v[i]=0;
    }

  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }  
  return(v);  
}





bool double_equals(double a, double b )
{
  return std::abs(a - b) < 0.001;
}


