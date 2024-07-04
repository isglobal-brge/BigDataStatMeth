#ifndef BIGDATASTATMETH_OMICS_UTILS_HPP
#define BIGDATASTATMETH_OMICS_UTILS_HPP


namespace BigDataStatMeth {


    double calc_freq(Rcpp::NumericVector x)
    {
        
        int len = x.size();
        
        std::vector<double> xc = Rcpp::as<std::vector<double> >(x);
        
        int n0 = std::count (xc.begin(), xc.end(), 0);
        int n1 = std::count (xc.begin(), xc.end(), 1);
        
        double maf = (double(n0)/len) + 0.5*(double(n1)/len);
        
        if( maf > 0.5 ) { 
            maf = 1 - maf;
        }
        
        return maf;
        
    }

}

#endif // BIGDATASTATMETH_OMICS_UTILS_HPP