#include "include/vector_utilities.h"


Eigen::VectorXd generate_seq( double start, double end, double inc)
{
  int ilen = std::ceil((end - start)/inc) + 1;
  Eigen::VectorXd seq(ilen);
  double val = start;
  
  seq[0] = start;
  for (size_t i= 1; val<=end; i++)
  {
    val = val + inc;
    seq[i] = val;
  }
  return(seq);
}
