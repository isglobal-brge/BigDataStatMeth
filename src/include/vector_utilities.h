#ifndef vector_utilities
#define vector_utilities

  #include <RcppEigen.h>
  #include <RcppParallel.h>
    #include <vector>

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


    // Return index where value val is found in array
    template<class _InputIterator, class T>
    std::vector<_InputIterator>
    find_all(_InputIterator begin, _InputIterator end, const T& val)
    {
        std::vector<_InputIterator> matches;
        while(begin != end)
        {
            if((*begin) == val)
                matches.push_back(begin);
            ++begin;
        }
        
        return matches;
    }
    
    
  void replace_zero(Rcpp::NumericVector* v);
  int index_val( double val, Eigen::VectorXd& vector);
  double bdparallelVectorSum(Rcpp::NumericVector x);
  Rcpp::NumericVector bdparallelpow2(Rcpp::NumericVector x);
  template <class T> T generate_seq (double start, double end, double inc);
  


#endif
