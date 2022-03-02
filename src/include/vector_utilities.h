#ifndef vector_utilities
#define vector_utilities

  #include <RcppEigen.h>
  #include <RcppParallel.h>

  template <class T>
  T generate_seq (double start, double end, double inc)
  {
    int ilen = std::ceil((end - start)/inc) + 1;
    T seq(ilen);
    double val = start;
    
    seq[0] = start;
    for (size_t i= 1; val<=end; i++)
    {
      val = val + inc;
      seq[i] = val;
    }
    return(seq);
  }


  void replace_zero(Rcpp::NumericVector* v);
  int index_val( double val, Eigen::VectorXd& vector);
  double bdparallelVectorSum(Rcpp::NumericVector x);
  Rcpp::NumericVector bdparallelpow2(Rcpp::NumericVector x);
  template <class T> T generate_seq (double start, double end, double inc);
  


#endif
