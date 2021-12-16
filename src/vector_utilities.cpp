#include "include/vector_utilities.h"


void replace_zero(Rcpp::NumericVector* v)
{
  for(size_t i=0; i < v->size(); i = i+1)
  {
    if ( v->at(i)<0 || Rcpp::NumericVector::is_na(v->at(i))==true) 
      v->at(i) = 0;
  }
}



int index_val( double val, Eigen::VectorXd& vector)
{
  size_t i;
  for( i = 0; i<vector.size() && vector[i]!=val; i = i+1 )
  {  }
  return(i);
}



struct Sum : public RcppParallel::Worker
{   
  // source
  const RcppParallel::RVector<double> input;
  
  // accumulated value
  double value;
  
  // constructors
  Sum(const Rcpp::NumericVector input) : input(input), value(0) {}
  Sum(const Sum& sum, RcppParallel::Split) : input(sum.input), value(0) {}
  
  // accumulate elements of range
  void operator()(std::size_t begin, std::size_t end) {
    value += std::accumulate(input.begin() + begin, input.begin() + end, 0.0);
  }
  
  // join values 
  void join(const Sum& rhs) { 
    value += rhs.value; 
  }
};



//' Sumarize vector
//' 
//' This function sumarize the elements of a vector
//' 
//' @param x numerical vector
//' @examples
//' library(BigDataStatMeth)
//' 
//' n <- 100 
//' x <- rnorm(n)
//' 
//' # with numeric matrix
//' res <- bdparallelVectorSum(x)
//' 
//' @return none value returned, data are stored in a dataset inside an hdf5 data file.
//' 
//' @export
// [[Rcpp::export]]
double bdparallelVectorSum(Rcpp::NumericVector x) {
  
  // declare instance 
  Sum sum(x);
  
  // start the work
  parallelReduce(0, x.length(), sum);
  
  return sum.value;
}


struct Pow : public RcppParallel::Worker
{   
  // source
  const RcppParallel::RVector<double> input;
  
  // output vector to write to
  RcppParallel::RVector<double> rv;
  
  
  // constructors
  Pow(const Rcpp::NumericVector input, Rcpp::NumericVector rv) 
    : input(input), rv(rv) {}
  
  // pow2
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) 
    {
      rv[i] = input[i] * input[i];
    }
  }
  
};

//' Pow vector
//' 
//' Gets pow2 vector
//' 
//' @param x numerical vector
//' @examples
//' library(BigDataStatMeth)
//' 
//' n <- 100 
//' x <- rnorm(n)
//' 
//' # with numeric matrix
//' res <- bdparallelpow2(x)
//' 
//' @return Numeric Vector
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector bdparallelpow2(Rcpp::NumericVector x) {
  
  // Return vector
  Rcpp::NumericVector rv(x.length());
  
  // declare instance 
  Pow pow(x, rv);
  
  // start the work
  parallelFor(0, x.length(), pow);
  
  return rv;
}