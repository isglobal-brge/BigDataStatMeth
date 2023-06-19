#ifndef BIGDATASTATMETH_HDF5_MATRIXSNORMALIZATION_HPP
#define BIGDATASTATMETH_HDF5_MATRIXSNORMALIZATION_HPP

#include <RcppEigen.h>
#include "Utilities/openme-utils.hpp"
#include "hdf5Utilities/hdf5Utilities.hpp"
#include "H5Cpp.h"

namespace BigDataStatMeth {



    // Internal call - 
    //   To be used when we have SD and Mean computed for each row/column and we 
    //   use this data to compute normalized matrix
    
    //.. ORIGINAL name ..// Eigen::MatrixXd RcppNormalize_Data_hdf5 ( Eigen::MatrixXd  X, bool bc, bool bs, bool btransp, Eigen::MatrixXd normdata )
    extern inline Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X, bool bc, bool bs, bool btransp, Eigen::MatrixXd normdata )
    {
        Eigen::MatrixXd rX;

        if( btransp == true)
        {
            if( bc==true && bs==true )  {

                rX = (X.colwise() - normdata.row(0).transpose() ).array().colwise() / normdata.row(1).transpose().array();

            }   else if (bc == true  && bs==false)   {

                Eigen::VectorXd mean = X.rowwise().mean();
                rX = (X.colwise() - normdata.row(0).transpose());

            }  else if ( bc == false && bs == true)   {
                rX = X.array().colwise() / normdata.row(1).transpose().array();
            }

        } else {

            if( bc==true && bs==true )  {

                Eigen::RowVectorXd mean = X.colwise().mean();
                Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
                rX = (X.rowwise() - mean).array().rowwise() / std.array();

            }   else if (bc == true  && bs==false)   {

                Eigen::RowVectorXd mean = X.colwise().mean();
                rX = (X.rowwise() - mean);

            }  else if ( bc == false && bs == true)   {

                Eigen::RowVectorXd mean = X.colwise().mean();
                Eigen::RowVectorXd std = (X.array().square().colwise().sum() / (X.rows() - 1)).sqrt();
                rX = X.array().rowwise() / std.array();
            }
        }

        return(rX);
    }


    // Internal call - 
    //   To be used when we don't have SD and Mean computed and we 
    //   to compute this data to get normalized matrix
    template< typename M>
    extern inline M RcppNormalize_Data ( M  X, bool bc, bool bs )
    {
        static_assert(std::is_same<M, Eigen::MatrixXd >::value || 
                      std::is_same<M, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> >::value || 
                      std::is_same<M, Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> >::value,
                      "Error - type not allowed");
        
        Eigen::MatrixXd rX = X; // use data as colmajor
        
        if( bc==true && bs==true )  {
            
            Eigen::RowVectorXd mean = X.colwise().mean();
            Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
            rX = (X.rowwise() - mean).array().rowwise() / std.array();
            
        }   else if (bc == true  && bs==false)   {
            
            Eigen::RowVectorXd mean = X.colwise().mean();
            rX = (X.rowwise() - mean);
            
        }  else if ( bc == false && bs == true)   {
            
            Eigen::RowVectorXd mean = X.colwise().mean();
            Eigen::RowVectorXd std = (X.array().square().colwise().sum() / (X.rows() - 1)).sqrt();
            rX = X.array().rowwise() / std.array();
        } 
        
        X = rX;
        return(X); // Return data as initial type
    };


    // // Internal call - 
    // //   To be used when we don't have SD and Mean computed and we 
    // //   to compute this data to get normalized matrix
    // Eigen::MatrixXd RcppNormalize_Data ( Eigen::MatrixXd  X, bool bc, bool bs )
    // {
    //     Eigen::MatrixXd rX;
    //     
    //     if( bc==true && bs==true )  {
    //         
    //         Eigen::RowVectorXd mean = X.colwise().mean();
    //         Eigen::RowVectorXd std = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
    //         rX = (X.rowwise() - mean).array().rowwise() / std.array();
    //         
    //     }   else if (bc == true  && bs==false)   {
    //         
    //         Eigen::RowVectorXd mean = X.colwise().mean();
    //         rX = (X.rowwise() - mean);
    //         
    //     }  else if ( bc == false && bs == true)   {
    //         
    //         Eigen::RowVectorXd mean = X.colwise().mean();
    //         Eigen::RowVectorXd std = (X.array().square().colwise().sum() / (X.rows() - 1)).sqrt();
    //         rX = X.array().rowwise() / std.array();
    //     } 
    //     
    //     return(rX);
    // }


    extern inline Eigen::MatrixXd RcppNormalize_Data_R_hdf5( Eigen::MatrixXd  X, bool bc, bool bs,
                                               bool btransp, Eigen::MatrixXd normdata)
    {
        Eigen::MatrixXd rX;

        if( btransp == true) {
            if( bc==true && bs==true )  {
                rX = (X.colwise() - normdata.row(0).transpose() ).array().colwise() / normdata.row(1).transpose().array();
            }   else if (bc == true  && bs==false)   {
                rX = (X.colwise() - normdata.row(0).transpose());
            }  else if ( bc == false && bs == true)   {
                rX = X.array().colwise() / normdata.row(1).transpose().array();
            }
        } else {
            if( bc==true && bs==true )  {
                rX = (X.rowwise() - normdata.row(0)).array().rowwise() / normdata.row(1).array();
            } else if (bc == true  && bs==false) {
                rX = (X.rowwise() - normdata.row(0));
            }  else if ( bc == false && bs == true)   {
                rX = X.array().rowwise() / normdata.row(1).array();
            }
        }

        return(rX);
    }



}

#endif // BIGDATASTATMETH_HDF5_MATRIXSNORMALIZATION_HPP

