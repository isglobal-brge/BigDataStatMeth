#ifndef BIGDATASTATMETH_ALGEBRA_MEM_OPTIMIZED_PRODS_HPP
#define BIGDATASTATMETH_ALGEBRA_MEM_OPTIMIZED_PRODS_HPP

// #include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
// #include <thread>

namespace BigDataStatMeth {

// extern inline Eigen::MatrixXd block_matrix_mul( Eigen::MatrixXd& A, Eigen::MatrixXd& B, int block_size);
// extern inline Eigen::MatrixXd Bblock_matrix_mul_parallel( Eigen::MatrixXd& A, Eigen::MatrixXd& B, 
//                                                                           int block_size, Rcpp::Nullable<int> threads  = R_NilValue);


// extern inline Eigen::MatrixXd bdcrossproduct (Eigen::MatrixXd& mat);
// extern inline Eigen::MatrixXd bdtcrossproduct (Eigen::MatrixXd& mat);
template< typename T> extern inline Eigen::MatrixXd bdcrossproduct ( T X );
template< typename T> extern inline Eigen::MatrixXd bdtcrossproduct ( T X );

extern inline Eigen::MatrixXd xwxt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w);
extern inline Eigen::MatrixXd xtwx(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w);
extern inline Eigen::MatrixXd Xwd(const Eigen::MatrixXd& X, const Eigen::VectorXd& w);
extern inline Eigen::MatrixXd Xwd_parallel(const Eigen::MatrixXd& X, const Eigen::VectorXd& w, Rcpp::Nullable<int> threads);
extern inline Eigen::MatrixXd wdX_parallel(const Eigen::MatrixXd& X, const Eigen::VectorXd& w, Rcpp::Nullable<int> threads);


// Computes CrossProduct X'X
// extern inline Eigen::MatrixXd bdcrossproduct (Eigen::MatrixXd& X)
// {
//     size_t nc(X.cols());
//     Eigen::MatrixXd XtX(Eigen::MatrixXd(nc, nc).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint()));
//     return(XtX);
// }

template< typename T>
extern inline Eigen::MatrixXd bdcrossproduct ( T X )
{
    
    static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
                  std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                  std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value ||
                  std::is_same<T, Eigen::Transpose<Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> > >::value ||
                  std::is_same<T, Eigen::Transpose<Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> > >::value,
                  "Error - type not allowed");
    
    Eigen::MatrixXd Xem = X;
    size_t nc(Xem.cols());
    Eigen::MatrixXd XtX(Eigen::MatrixXd(nc, nc).setZero().selfadjointView<Eigen::Lower>().rankUpdate(Xem.adjoint()));
    return(XtX);
}




// Compute tCrossProduct XX'
// extern inline Eigen::MatrixXd bdtcrossproduct (Eigen::MatrixXd& X)
// {
//     
//     size_t nr(X.rows());
//     Eigen::MatrixXd XXt(Eigen::MatrixXd(nr, nr).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X));
//     return(XXt);
// }


template< typename T>
extern inline Eigen::MatrixXd bdtcrossproduct ( T X )
{
    
    static_assert(std::is_same<T, Eigen::MatrixXd >::value || 
                  std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                  std::is_same<T, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value ||
                  std::is_same<T, Eigen::Transpose<Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> > >::value ||
                  std::is_same<T, Eigen::Transpose<Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> > >::value,
                  "Error - type not allowed");
    
    Eigen::MatrixXd Xem = X;
    
    size_t nr(Xem.rows());
    Eigen::MatrixXd XXt(Eigen::MatrixXd(nr, nr).setZero().selfadjointView<Eigen::Lower>().rankUpdate(Xem));
    return(XXt);
}



// Compute weighted crossproduct XwX'
extern inline Eigen::MatrixXd xwxt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w) 
{
    const int n(X.rows());
    Eigen::MatrixXd XwXt(Eigen::MatrixXd(n, n).setZero().
                             selfadjointView<Eigen::Lower>().rankUpdate(X * w.array().sqrt().matrix().asDiagonal()));
    return (XwXt);
}


// Compute transposed weighted crossproduct X'wX 
extern inline Eigen::MatrixXd xtwx(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w)
{
    const int n(X.cols());
    Eigen::MatrixXd XtwX(Eigen::MatrixXd(n, n).setZero().
                             selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint() * w.array().sqrt().matrix().asDiagonal()));
    return (XtwX);
}


// Compute weighted crossproduct Xw
extern inline Eigen::MatrixXd Xw(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w) 
{
    Eigen::MatrixXd Xw = X * w.array().matrix().asDiagonal();
    return (Xw);
}


// Compute weighted crossproduct Xw
extern inline Eigen::MatrixXd wX(const Eigen::MatrixXd& X, const Eigen::MatrixXd& w) 
{
    Eigen::MatrixXd wX = w.array().matrix().asDiagonal()*X;
    return (wX);
}



// Matirx - vector as diagonal matrix Multiplication
extern inline Eigen::MatrixXd Xwd(const Eigen::MatrixXd& X, const Eigen::VectorXd& w)
{
    int n = X.rows();
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n,X.cols()) ; 
    
    for (int i=0; i<n; i++) {
        C.col(i) = X.col(i)*w(i);
    }
    return(C);
}


//vector as diagonal matrix -  Matirx Multiplication
extern inline Eigen::MatrixXd wdX(const Eigen::MatrixXd& X, const Eigen::VectorXd& w)
{
    int n = X.cols();
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(X.rows(),n) ; 
    
    for (int i=0; i<n; i++) {
        C.row(i) = w(i)*X.row(i);
    }
    return(C);
}


// Matirx - vector as diagonal matrix Multiplication (parallel)
extern inline Eigen::MatrixXd Xwd_parallel(const Eigen::MatrixXd& X, const Eigen::VectorXd& w, Rcpp::Nullable<int> threads = R_NilValue)
{
    int n = X.rows();
    unsigned int ithreads;
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n,X.cols()) ; 
    
    if(threads.isNotNull()) {
        if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
            ithreads = Rcpp::as<int> (threads);
        } else {
            ithreads = getDTthreads(0, true);
            //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;}
        }
    } else {
        ithreads = getDTthreads(0, true);
        //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;
    }
    
    //.OpenMP.//omp_set_num_threads(ithreads);
    
    //.OpenMP.//#pragma omp parallel shared(X, w, C) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(X, w, C) 
{
#pragma omp for schedule (dynamic)
    for (int i=0; i<n; i++)
    {
        C.col(i) = X.col(i)*w(i);
    }  
}
return(C);
}


//vector as diagonal matrix -  Matirx Multiplication (parallel)
extern inline Eigen::MatrixXd wdX_parallel(const Eigen::MatrixXd& X, const Eigen::VectorXd& w, Rcpp::Nullable<int> threads = R_NilValue)
{
    int n = X.cols();
    unsigned int ithreads;
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(X.rows(),n);
    
    if(threads.isNotNull()) {
        if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
            ithreads = Rcpp::as<int> (threads);
        } else {
            ithreads = getDTthreads(0, true);
            //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;}
        }
    } else {
        ithreads = getDTthreads(0, true);
        //.11-04-2022.// ithreads = std::thread::hardware_concurrency()/2;
    }
    
    //.OpenMP.// omp_set_num_threads(ithreads);
    
    //.OpenMP.//#pragma omp parallel shared(X, w, C) 
#pragma omp parallel num_threads(getDTthreads(ithreads, true)) shared(X, w, C) 
{
#pragma omp for schedule (dynamic)
    for (int i=0; i<n; i++)
    {
        C.row(i) = w(i)*X.row(i);
    }
}
return(C);
}

}

#endif // BIGDATASTATMETH_ALGEBRA_MEM_OPTIMIZED_PRODS_HPP