#include <BigDataStatMeth.hpp>
#include "hdf5Algebra/matrixEquationSolver.hpp"
#include "hdf5Algebra/matrixPseudoinverse.hpp"
#include "Utilities/Utilities.hpp"



//' Solve matrix equations
//' 
//' This function solve matrix equations 
//'  \code{A * X = B } 
//' where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//' 
//' @param A numerical matrix. 
//' @param B numerical matrix.
//' @return X numerical matrix. 
//' @examples
//' 
//' library(BigDataStatMeth)
//' 
//' n <- 500
//' m <- 500
//' 
//' # R Object
//' 
//' A <- matrix(runif(n*m), nrow = n, ncol = m)
//' B <- matrix(runif(n), nrow = n)
//' AS <- A%*%t(A)
//'       
//' X <- bdSolve(A, B)
//' XR <- solve(A,B)
//'       
//' all.equal(X, XR, check.attributes=FALSE)
//'   
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdSolve(const Rcpp::RObject A, const Rcpp::RObject B) 
{
    
    try {
        
        Eigen::MatrixXd a, b;
        
        char Uchar = 'U';
        int info = 0;
        
        a = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(A);
        b = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(B);

        // Declare matrix variables
        int n = a.rows();
        int nrhs = b.cols();
        int lwork = std::max( 1, n );
        int lda = std::max( 1, n );
        int ldb = std::max( 1, n );
        std::vector<int> ipiv(n);
        std::vector<double> work(lwork);
        
        // Solve matrix equation
        if( a == a.transpose()  )
        {
            
            // dsysv_( char* UPLO, int* N , int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, double* WORK, int* LWORK, int* INFO);
            BigDataStatMeth::dsysv_( & Uchar, &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, work.data(), &lwork, &info);
        } else {
            
            // dgesv( int N, int NRHS, double A, int LDA, int IPIV, double B, int LDB, int INFO);
            BigDataStatMeth::dgesv_( &n, &nrhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, &info );
        }
    
        return(Rcpp::wrap(b));
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return Rcpp::wrap(-1);
    }
    
}



/********************+
 * 
 * ESTIC AQUÍ !!! 
 * 
 * HE DE PROGRAMAR EL CÀLCUL UTILITZANT FITXERS HDF5!!!!!
 * FALTA TESTEJAR LA SOLUCIÓ EN MEMÒRIA
 * 
 * ******************/


//' Solve matrix equations
//' 
//' This function solve matrix equations 
//'  \code{A * X = B } 
//' where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//' 
//' @param A numerical matrix. 
//' @param B numerical matrix.
//' @return X numerical matrix. 
//' @examples
//' 
//' library(BigDataStatMeth)
//' 
//' n <- 500
//' m <- 500
//' 
//' # R Object
//' 
//' A <- matrix(runif(n*m), nrow = n, ncol = m)
//' B <- matrix(runif(n), nrow = n)
//' AS <- A%*%t(A)
//'       
//' X <- bdSolve(A, B)
//' XR <- solve(A,B)
//'       
//' all.equal(X, XR, check.attributes=FALSE)
//'   
//' @export
// [[Rcpp::export]]
void bdSolve_hdf5(std::string filename, std::string groupA, std::string datasetA,
                      std::string groupB, std::string datasetB,
                      Rcpp::Nullable<std::string> outgroup = R_NilValue, 
                      Rcpp::Nullable<std::string> outdataset = R_NilValue, 
                      Rcpp::Nullable<bool> overwrite = R_NilValue) 
{
    
    try {
        
        std::string strOutgroup, strOutdataset;
        bool bforce;
        
        if(outgroup.isNull()) { strOutgroup = "Solved"; } 
        else {   strOutgroup = Rcpp::as<std::string>(outgroup); }
        
        if(outdataset.isNull()) { strOutdataset = datasetA+"_"+datasetB; } 
        else { strOutdataset = Rcpp::as<std::string>(outdataset); }
        
        if(overwrite.isNull()) { bforce = false ; }
        else { bforce = Rcpp::as<bool>(overwrite); }
        
        BigDataStatMeth::hdf5Dataset* dsA = new BigDataStatMeth::hdf5Dataset(filename, groupA, datasetA, false);
        dsA->openDataset();
        
        BigDataStatMeth::hdf5Dataset* dsB = new BigDataStatMeth::hdf5Dataset(filename, groupB, datasetB, false);
        dsB->openDataset();
        
        BigDataStatMeth::hdf5Dataset* dsRes = new BigDataStatMeth::hdf5Dataset(filename, strOutgroup, strOutdataset, bforce);
        dsRes->createDataset( dsB->nrows(), dsB->ncols(), "real" );
        
        RcppSolveHdf5(dsA, dsB, dsRes );
        
        delete dsRes;
        delete dsB;
        delete dsA;
        
    } catch(std::exception &ex) {
        Rcpp::Rcout<< ex.what();
        return void();
    }
    return void();
     
 }