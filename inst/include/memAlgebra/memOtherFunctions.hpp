#ifndef BIGDATASTATMETH_ALGEBRA_MEM_OTHER_FUNCTIONS_HPP
#define BIGDATASTATMETH_ALGEBRA_MEM_OTHER_FUNCTIONS_HPP

// #include <RcppEigen.h>
// #include "Utilities/openme-utils.hpp"
// #include <thread>

namespace BigDataStatMeth {


    extern inline void getBlockPositionsSizes( hsize_t maxPosition, hsize_t blockSize, std::vector<hsize_t>& starts, std::vector<hsize_t>& sizes ){
        
        hsize_t isize = blockSize + 1;
        
        for (hsize_t ii = 0; ii < maxPosition; ii += blockSize)
        {
            if( ii + blockSize > maxPosition ) {
                isize = maxPosition - ii; }
            
            hsize_t sizetoRead = getOptimBlockSize( maxPosition, blockSize, ii, isize);
            
            starts.push_back(ii);
            sizes.push_back(sizetoRead);
            
            if( ii + blockSize > maxPosition ) isize = blockSize + 1;
            if( sizetoRead > blockSize ) {
                ii = ii - blockSize + sizetoRead; }
        }
        
    }
    
    
    extern inline Eigen::VectorXd cumsum(Eigen::VectorXd x)
    {
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

}

#endif // BIGDATASTATMETH_ALGEBRA_MEM_OTHER_FUNCTIONS_HPP