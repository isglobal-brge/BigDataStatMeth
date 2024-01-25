#include <BigDataStatMeth.hpp>
#include "memAlgebra/memOptimizedProducts.hpp"
#include "hdf5Algebra/matrixQR.hpp"
#include "hdf5Algebra/multiplication.hpp"
#include "hdf5Algebra/matrixInvCholesky.hpp"
#include "hdf5Algebra/matrixEquationSolver.hpp"
#include "hdf5Algebra/matrixSdMean.hpp"

//' Apply function to different datasets inside a group
//'
//' Apply function to different datasets inside a group
//' 
//' @param filename, Character array, indicating the name of the file to create
//' @param group, Character array, indicating the input group where the data set
//' to be imputed is. 
//' @param datasets, Character array, indicating the input datasets to be used
//' @param outgroup, Character, array, indicating group where the data set will 
//' be saved after imputation if `outgroup` is NULL, output dataset is stored 
//' in the same input group. 
//' @param func, Character array, function to be applyed : 
//' QR to apply bdQR() function to datasets
//' CrossProd to apply bdCrossprod() function to datasets
//' tCrossProd to apply bdtCrossprod() function to datasets
//' invChol to apply bdInvCholesky() function to datasets
//' blockmult to apply matrix multiplication, in that case, we need the datasets 
//' to be used defined in b_datasets variable, datasets and b_datasets must be 
//' of the same lenght, in that case, the operation is performed according to 
//' index, for example, if we have datasets = {"A1", "A2", "A3} and 
//' b_datasets = {"B1", "B2", "B3}, the functions performs : A1%*%B1, 
//' A2%*%B2 and A3%*%B3 
//' CrossProd_double to  performs crossprod using two matrices, see blockmult 
//' tCrossProd_double to  performs transposed crossprod using two matrices, 
//' see blockmult 
//' solve to solve matrix equation system, see blockmult for parametrization 
//' sdmean to get sd and mean from de datasets by cols or rows
//' @param b_group, optional Character array indicating the input group where 
//' data are stored when we need a second dataset to operate, for example in 
//' functions like matrix multiplication
//' @param b_datasets, optional Character array indicating the input datasets 
//' to be used when we need a second dataset in functions like matrix 
//' multiplication
//' @param force, optional Boolean if true, previous results in same location 
//' inside hdf5 will be overwritten, by default force = false, data was not 
//' overwritten.
//' @param transp_dataset optional parameter. Boolean if true we use the 
//' transposed dataframe to perform calculus. By default transp_dataset = false, 
//' we use the original dataset stored in hdf5 data file. Currently this option 
//' is only valid with "blockmult", "CrossProd_double" and "tCrossProd_double"
//' @param transp_bdataset optional parameter. Boolean if true we use the 
//' transposed dataframe to perform calculus.By default transp_bdataset = false, 
//' we use the original dataset stored in hdf5 data file. Currently this option 
//' is only valid with "blockmult", "CrossProd_double" and "tCrossProd_double"
//' @param fullMatrix boolean, optional parameter used in Inverse Cholesky, by 
//' default false. If fullMatrix = true, in the hdf5 file the complete matrix 
//' is stored. If false, only the lower triangular matrix is stored
//' @param byrows boolean, optional parameter used in sd and mean calculus, by 
//' default false. If byrows = true, the sd and mean is computed by columns. 
//' If false, sd and mean is computed by rows.
//' @param threads optional parameter. Integer with numbers of threads to be used
//' @return Original hdf5 data file with results after apply function to 
//' different datasets
//' @export
// [[Rcpp::export]]
void bdapply_Function_hdf5( std::string filename, 
                            std::string group, 
                            Rcpp::StringVector datasets, 
                            std::string outgroup, 
                            std::string func, 
                            Rcpp::Nullable<std::string> b_group = R_NilValue, 
                            Rcpp::Nullable<Rcpp::StringVector> b_datasets = R_NilValue,
                            Rcpp::Nullable<bool> force = false,
                            Rcpp::Nullable<bool> transp_dataset = false,
                            Rcpp::Nullable<bool> transp_bdataset = false,
                            Rcpp::Nullable<bool> fullMatrix = false,
                            Rcpp::Nullable<bool> byrows = false,
                            Rcpp::Nullable<int> threads = 2 )
{
    
    Rcpp::StringVector str_bdatasets;
    std::string str_bgroup, outputdataset;
    bool btransdataA, btransdataB, bfullMatrix, bbyrows;
    long dElementsBlock = MAXELEMSINBLOCK; 
    Rcpp::NumericVector oper = {0, 1, 2, 3, 4, 11, 22, 5, 6, 7, 8};
    oper.names() = Rcpp::CharacterVector({"QR", "CrossProd", "tCrossProd",
               "invChol", "blockmult", "CrossProd_double", "tCrossProd_double",
               "solve", "normalize", "sdmean", "descChol"});
    
    try
    {
        
        bool bforce;
        
        std::vector<hsize_t> stride = {1, 1},
                             block = {1, 1};
        
        if(force.isNull()) { bforce = false; } 
        else {   bforce = Rcpp::as<bool>(force); }
        
        if(transp_dataset.isNull()) { btransdataA = false; } 
        else {   btransdataA = Rcpp::as<bool>(transp_dataset); }
        
        if(transp_bdataset.isNull()) { btransdataB = false; } 
        else {   btransdataB = Rcpp::as<bool>(transp_bdataset); }
        
        if(fullMatrix.isNull()) { bfullMatrix = false; } 
        else {   bfullMatrix = Rcpp::as<bool>(fullMatrix); }
        
        if(byrows.isNull()) { bbyrows = false; } 
        else {   bbyrows = Rcpp::as<bool>(byrows); }
        
        if(b_group.isNull()) { str_bgroup = group; } 
        else {   str_bgroup = Rcpp::as<std::string>(b_group); }
        
        if( b_datasets.isNotNull() &&  ( oper(oper.findName(func)) == 1 ||  
            oper(oper.findName(func)) == 2 || oper(oper.findName(func)) == 4 || 
            oper(oper.findName(func)) == 5 || oper(oper.findName(func)) == 6) ) {
            
            str_bdatasets = Rcpp::as<Rcpp::StringVector>(b_datasets);
            
            if( oper.findName( func ) == 1) {
                func = "CrossProd_double";
            } else if( oper.findName( func ) == 2) {
                func = "tCrossProd_double";
            }
        }
        
        if( b_group.isNull()) { 
            str_bgroup = group; 
        } else {   
            str_bgroup = Rcpp::as<std::string>(b_group); 
        }
        
        // Seek all datasets to perform calculus
        for( int i=0; i < datasets.size(); i++ ) 
        {
            
            BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, group, Rcpp::as<std::string>(datasets(i)), false);
            dsA->openDataset();

            if( oper(oper.findName( func )) == 0) {
                // ==> QR Decomposition 
                
                BigDataStatMeth::hdf5Dataset* dsQ = new BigDataStatMeth::hdf5Dataset(filename, outgroup, Rcpp::as<std::string>(datasets(i)) + ".Q", bforce);
                BigDataStatMeth::hdf5Dataset* dsR = new BigDataStatMeth::hdf5Dataset(filename, outgroup, Rcpp::as<std::string>(datasets(i)) + ".R", bforce);
               
                RcppQRHdf5(dsA, dsQ, dsR, true, R_NilValue, threads);
                
                delete dsA;
                delete dsQ;
                delete dsR;
                
            } else if( oper(oper.findName( func )) == 1 || oper(oper.findName( func )) == 2) {
                // ==> CrossProd and transposed CrossProd
                
                //  PREPARAR LA FUNCIÓ PER A QUE TREBALLI EN MEMÒRIA SI DADES PETITES
                //  I PER BLOCKS DIRECTAMETN DES DE DE FITXER SI DADES GRANS ???!!!!
                //  TENIR-HO EN TCOMPTE QUAN TOT OK
                
                hsize_t* dims_out = dsA->dim();
                
                std::vector<double> vdA( dims_out[0] * dims_out[1] ); 
                dsA->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, stride, block, vdA.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> original (vdA.data(), dims_out[1], dims_out[0] );
                
                delete dsA;
                
                Eigen::MatrixXd results;
                BigDataStatMeth::hdf5Dataset* dsOut = new BigDataStatMeth::hdf5Dataset(filename, outgroup, Rcpp::as<std::string>(datasets(i)) , bforce);
                
                if(  oper(oper.findName( func )) == 1 ) {
                    results = BigDataStatMeth::bdcrossproduct(original);
                } else {
                     results = BigDataStatMeth::bdtcrossproduct(original);
                }
                
                dsOut->createDataset(results.rows(), results.cols(), "numeric"); 
                dsOut->openDataset(); 
                dsOut->writeDataset(Rcpp::wrap(results));
                
                delete dsOut;
                
            } else if( oper(oper.findName( func )) == 3 || oper(oper.findName( func )) == 8) {
                // ==> Inverse Cholesky and Cholesky decomposition

                int nrows = dsA->nrows();
                int ncols = dsA->ncols();
                
                if(nrows == ncols) {
                    
                    BigDataStatMeth::hdf5DatasetInternal* dsOut = new BigDataStatMeth::hdf5DatasetInternal(filename, outgroup, Rcpp::as<std::string>(datasets(i)) , bforce);
                    dsOut->createDataset(nrows, ncols, "real");
                    
                    if(oper(oper.findName( func )) == 3 ) {
                        BigDataStatMeth::Rcpp_InvCholesky_hdf5( dsA, dsOut, bfullMatrix, dElementsBlock, threads);    
                    } else {
                        int res = BigDataStatMeth::Cholesky_decomposition_hdf5( dsA, dsOut,  nrows, ncols, dElementsBlock, threads);
                    }
                    
                    delete dsA;
                    delete dsOut;
                    
                } else {
                    delete dsA;
                    Rcpp::Rcout<<"\n Can't get inverse matrix for "<<Rcpp::as<std::string>(datasets(i))<<" dataset using Cholesky decomposition - not an square matrix\n";
                    return void();
                }
                

            } else if( oper(oper.findName( func )) == 4 ||  oper(oper.findName( func )) == 11 ||  oper(oper.findName( func )) == 22) {
                // ==> blockmult, CrossProd Double, tCrossProd Double"
                
                BigDataStatMeth::hdf5Dataset* dsB = new BigDataStatMeth::hdf5Dataset(filename, str_bgroup, Rcpp::as<std::string>(str_bdatasets(i)), false);
                dsB->openDataset();
                
                // Real data set dimension
                hsize_t* dims_outB = dsB->dim();
                hsize_t* dims_out = dsA->dim();
                
                std::vector<double> vdA( dims_out[0] * dims_out[1] ); 
                dsA->readDatasetBlock( {0, 0}, {dims_out[0], dims_out[1]}, stride, block, vdA.data() );
                
                std::vector<double> vdB( dims_outB[0] * dims_outB[1] ); 
                dsB->readDatasetBlock( {0, 0}, {dims_outB[0], dims_outB[1]}, stride, block, vdB.data() );
                
                Eigen::MatrixXd original;
                Eigen::MatrixXd originalB;
                
                if( oper(oper.findName( func )) == 4 ) {
                    original = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> (vdA.data(), dims_out[1], dims_out[0] );
                    originalB = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> (vdB.data(), dims_outB[1], dims_outB[0] );
                    outputdataset = Rcpp::as<std::string>(datasets(i)) + "_" + str_bdatasets(i);
                } else if  (oper(oper.findName( func )) == 11) {
                    original = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdA.data(), dims_out[0], dims_out[1] );
                    originalB = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> (vdB.data(), dims_outB[1], dims_outB[0] );
                    outputdataset = "Cross_" + datasets(i) + "_" + str_bdatasets(i);
                } else if ( oper(oper.findName( func )) == 22) {
                    original = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> (vdA.data(), dims_out[1], dims_out[0] );
                    originalB = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdB.data(), dims_outB[0], dims_outB[1] );
                    outputdataset =  "tCross_" + datasets(i) + "_" + str_bdatasets(i);
                }

                Eigen::MatrixXd results;

                if (btransdataA == false && btransdataB == false) {
                    results = BigDataStatMeth::Bblock_matrix_mul_parallel(original, originalB, 2048, R_NilValue);
                } else if (btransdataA == true && btransdataB == false) {
                    results = BigDataStatMeth::Bblock_matrix_mul_parallel(original.transpose(), originalB, 2048, R_NilValue);
                }else if (btransdataA == false && btransdataB == true) {
                    results = BigDataStatMeth::Bblock_matrix_mul_parallel(original, originalB.transpose(), 2048, R_NilValue);
                } else {
                    results = BigDataStatMeth::Bblock_matrix_mul_parallel(original.transpose(), originalB.transpose(), 2048, R_NilValue);
                }
                
                if( results != Eigen::MatrixXd::Zero(original.rows(),originalB.cols()) ) {
                    BigDataStatMeth::hdf5Dataset* dsOut = new BigDataStatMeth::hdf5Dataset(filename, outgroup, outputdataset , bforce);
                    dsOut->createDataset(results.rows(), results.rows(), "real");
                    dsOut->writeDataset(Rcpp::wrap(results));
                    
                    delete dsOut;    
                    
                } else {
                    Rcpp::Rcout<<"Multiplication: "<< group<<"/"<< Rcpp::as<std::string>(datasets(i))<< " x "<< str_bgroup<<"/"<< Rcpp::as<std::string>(str_bdatasets(i)) <<" can not be computed \n";
                }
                
                delete dsA;
                delete dsB;
                
            } else if( oper(oper.findName( func )) == 5) {
                // ==> Solve matrix equation Ax = B
 
                BigDataStatMeth::hdf5Dataset* dsB = new BigDataStatMeth::hdf5Dataset(filename, str_bgroup, Rcpp::as<std::string>(str_bdatasets(i)), false);
                dsB->openDataset();
                
                BigDataStatMeth::hdf5DatasetInternal* dsOut = new BigDataStatMeth::hdf5DatasetInternal(filename, outgroup, Rcpp::as<std::string>(datasets(i)) + "_eq_" + Rcpp::as<std::string>(str_bdatasets(i)) , bforce);
                dsOut->createDataset( dsB->nrows(), dsB->ncols(), "real" );
                
                RcppSolveHdf5(dsA, dsB, dsOut );
                
                delete dsOut;
                delete dsB;
                delete dsA;
                

            } else if( oper(oper.findName( func )) == 7) {
                // ==> Compute sd and mean by rows or columns
                
                Eigen::MatrixXd datanormal;
                hsize_t* dims_out = dsA->dim();
                
                if( bbyrows == false) {
                    datanormal = Eigen::MatrixXd::Zero( 2, (int)dims_out[0]);
                    get_HDF5_mean_sd_by_column( dsA, datanormal, R_NilValue );

                } else {
                    datanormal = Eigen::MatrixXd::Zero( 2, (int)dims_out[1]);
                    get_HDF5_mean_sd_by_row( dsA, datanormal, R_NilValue );

                }
                
                BigDataStatMeth::hdf5Dataset* dsmean = new BigDataStatMeth::hdf5Dataset(filename, outgroup, "mean." + datasets(i) , bforce);
                dsmean->createDataset(1, datanormal.cols(), "real");
                dsmean->writeDataset(Rcpp::wrap(datanormal.row(0)));
                
                BigDataStatMeth::hdf5Dataset* dssd = new BigDataStatMeth::hdf5Dataset(filename, outgroup, "sd." + datasets(i) , bforce);
                dssd->createDataset(1, datanormal.cols(), "real");
                dssd->writeDataset(Rcpp::wrap(datanormal.row(1)));
                
                delete dssd;
                delete dsmean;
                delete dsA;
            
            } else {
                delete dsA;
                Rcpp::Rcout<< "Function does not exists, please use one of the following : \"QR\", \"CrossProd\","<<
                        " \"tCrossProd\", \"invChol\", \"blockmult\", \"CrossProd_double\", \"tCrossProd_double\","<<
                        " \"solve\", \"normalize\", \"sdmean\", \"descChol\" ";
                return void();
            }
        }
        
    } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
        ::Rf_error( "c++ exception bdapply_Function_hdf5 (File IException)" );
        return void();
    }
    
    // Rcpp::Rcout<< func <<" function has been computed in all blocks\n";  
    return void();
}

