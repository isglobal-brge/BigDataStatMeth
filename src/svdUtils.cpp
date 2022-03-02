#include "include/svdutils.h"



Rcpp::IntegerVector getInitialPosition(bool transp, int desp )
{
  Rcpp::IntegerVector voffset(2);
  
  if(transp == true)
  {
    voffset[0] = 0;
    voffset[1] = desp;
  } else {
    voffset[0] = desp;
    voffset[1] = 0;
  }
  
  return(voffset);
}


Rcpp::IntegerVector getSizetoRead(bool transp, int count, int rows, int cols )
{
  Rcpp::IntegerVector vcount(2);
  
  if(transp == true)
  {
    vcount[0] = rows;
    vcount[1] = count;
    
  } else {
    vcount[0] = count;
    vcount[1] = cols;
  }
  
  return(vcount);
}
