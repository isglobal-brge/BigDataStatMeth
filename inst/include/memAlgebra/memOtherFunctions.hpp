#ifndef BIGDATASTATMETH_ALGEBRA_MEM_OTHER_FUNCTIONS_HPP
#define BIGDATASTATMETH_ALGEBRA_MEM_OTHER_FUNCTIONS_HPP

#include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include <thread>

namespace BigDataStatMeth {

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