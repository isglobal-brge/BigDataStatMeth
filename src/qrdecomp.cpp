#include "include/qrdecomp.h"


using namespace std;


// 
//' QR Decomposition 
//' 
//' This function compute QR decomposition (also called a QR factorization) 
//' of a matrix \code{A} into a product \code{A = QR} of an 
//' orthogonal matrix Q and an upper triangular matrix R.
//' 
//' @param A a real square matrix 
//' @param boolean thin, if thin = true returns Q thin  decomposition else returns Q full decomposition, default thin = false
//' @return List with orthogonal matrix \code{Q}  and upper triangular matrix \codi{R}
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdQR(Eigen::MatrixXd xl, Rcpp::Nullable<bool> thin = R_NilValue)
{
  
  int m = xl.rows(), n = xl.cols();
  int irank;
  bool bthin;
  Eigen::MatrixXd R;
  Eigen::MatrixXd Q;
  
  Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(xl);
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(xl);
  
  
  qr.compute(xl);
  irank = lu_decomp.rank();
  
  if (irank == m + 1 || irank == n + 1 )
  {
    R = qr.matrixQR().template triangularView<Eigen::Upper>();
  } else {
    R = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>(); 
  }
  
  if( thin.isNull()) {
    bthin = false;
  } else {
    bthin = Rcpp::as<bool> (thin);
  }
  
  if (bthin == false)
  {
    Q =  qr.householderQ();       // Full decomposition
  } else {
    
    Q = Eigen::MatrixXd::Identity(m,n);
    Q = qr.householderQ() * Q;    // Thin decomposition
  }
  
  return Rcpp::List::create(Rcpp::Named("Q") = Q,
                            Rcpp::Named("R") = R
  );
}



// [[Rcpp::export]]
Rcpp::RObject review_decomposition(Eigen::MatrixXd R, int n)
{
  
  // if R is smaller than nxn
  while( R.cols() < n )
  {
    R.conservativeResize(R.rows()+1, R.cols()+1);
    R.col(R.cols()-1) = Eigen::VectorXd::Zero(R.rows());
    R.row(R.rows()-1) = Eigen::VectorXd::Zero(R.cols());
    R(R.rows()-1, R.cols()-1) = 0.0000000000000001;
  }
  
  
  if( (R.colwise().sum().array() == 0.0).any())
  {
    Eigen::VectorXd colsum = R.colwise().sum();
    for (int i = 0; i < n; i++){
      if( colsum(i) == 0 )   R(i, i) = 0.0000000000000001;
    }
  }
  
  return(Rcpp::wrap(R));
  
}





//' Solves matrix equations : A*X = B
//' 
//' 
//' 
//' @param R numerical or Delayed Array matrix. 
//' @param Z numerical or Delayed Array matrix.
//' @return X numerical matrix. 
//' @examples
//' 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bddtrsm(Rcpp::RObject R, Rcpp::RObject Z) 
{
  char Nchar='N';
  char Uchar='U';
  char Lchar='L'; //Nova def.
  double done = 1.0;
  
  
  Eigen::MatrixXd A = Rcpp::as<Eigen::MatrixXd> (R);
  Eigen::FullPivLU<Eigen::MatrixXd>matinv(A);
  Eigen::MatrixXd B = Rcpp::as<Eigen::MatrixXd> (Z);
  
  int m = B.rows(), n = B.cols();
  int lda = std::max( 1, m );
  int ldb = std::max( 1, m );
  
  if(matinv.isInvertible()==true)
  {
    
    // dtrsm( char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double A, int LDA, double B, int LDB )
    dtrsm_(&Lchar, &Uchar, &Nchar, &Nchar, &m, &n, &done, A.data(), &lda, B.data(), &ldb ); // La solució (X) es sobrescriu a (B) --> B es converteix en X
    return(Rcpp::wrap(B));
    
  } else {

    int block_size = 128;
    if( block_size > std::min(A.rows(),A.cols()) || block_size > std::min(B.rows(),B.cols()) )  {
      block_size = std::min(  std::min(A.rows(),A.cols()), std::min(B.rows(),B.cols()));
    }
    // Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd>( bdpseudoinv(A) )* B;
    Eigen::MatrixXd X = block_matrix_mul_parallel( Rcpp::as<Eigen::MatrixXd>(bdpseudoinv(A)), B, block_size);
    return(Rcpp::wrap(X));
  }
  
  return(0);
  
}






/***R
library(microbenchmark)
library(lapack)
library(BigDataStatMeth)
library(gdsfmt)
library(SNPRelate)

# Preparant dades .... 
X1 <- read.table( "data/colorectalnna.txt", header = TRUE )  # (1)
head(X1)
Y <- as.data.frame(X1[,c(6)]);   # Variable resposta bmi  --> Desprès es pot provar d'augmentar a bmi, meat i chol
colnames(Y) <- colnames(X1[6])
head(Y);
XL <- X1[,c(3:5)];  # Efectes fixats : sexe, edat i fumador
#write.table(XL, "data/colorectal2.txt", sep="\t", row.names = FALSE)
####

snpgdsBED2GDS("data/colorectal.bed", "data/colorectal.fam", 
              "data/colorectal.bim", "data/test.gds")

sortida <- Ols_Grid5("data/colorectal2.txt", "data/test.gds", as.matrix(Y), 5000, 10, 1,1)

sortida$beta[50:100]
sortida$est[50:100]
sortida$pvalue[50:100]
sortida$est

sortida$est

genofile <- openfn.gds("data/test.gds")




estimat <- get_pvalue( sortida$beta, genofile, Y )

estimat[1:10]

closefn.gds(genofile)

which(is.na(sortida$beta), arr.ind=TRUE)
which(is.na(Y), arr.ind=TRUE)
sortida$
  
  
  print(which(is.na(sortida$beta), arr.ind=TRUE)) # Detectar valors NA



BBikj[ rfrom:rto, 1:tbmaxsize] <- bddtrsm( RBRi[[k]], as.matrix(ZBij[ rfrom:rto, cfrom:cto ]))


as.matrix(sortida$ZBij[ 101:120, 1:1 ])
bdpseudoinv(sortida$RBRi[[7]])
solve( sortida$RBRi[[7]])

# Comprovacions



genofile <- openfn.gds("data/test.gds")

snps <- gdsfmt::objdesp.gdsn(index.gdsn(genofile, "genotype") )$dim[2]
n <- gdsfmt::objdesp.gdsn(index.gdsn(genofile, "genotype") )$dim[1]

xi <- read.gdsn( index.gdsn(genofile, "genotype"),
                 start = c( 1, 1 ),
                 count = c( n, snps ) )


est_vf <- BigDataStatMeth::blockmult(matrix(unlist(XL), ncol = 3, byrow = FALSE), 
                                     as.matrix( sortida$beta[ 1 : (dim(as.matrix(sortida$beta))[1] - dim(xi)[2]) ]), paral = TRUE)

cbind(Y,est_vf)[1:20,]

estl2[1:10,]

sorti <- xi[1:200,1:200] %*% as.matrix( sortida$beta[ 1:200])
sorti[1:10]

estr <- xi %*% as.matrix( sortida$beta[ (dim(as.matrix(sortida$beta))[1] + 1 - dim(xi)[2]) : dim(as.matrix(sortida$beta))[1]]) # Tenim estimacions

which(is.na(sortida$beta), arr.ind=TRUE)


est <- estl + estr
# Tenim originals

# print(which(is.na(xi), arr.ind=TRUE)) # Detectar valors NA

residuals <- as.matrix(Y) - est

RSS <- sum(residuals^2)

se2 <- 1 / sum( est - mean(as.matrix(Y)))^2


t <- betas / sqrt( se2 )
pvalue <- 2*pt(-abs(t), df=n-1)





bdQR()






BigDataStatMeth::bdInvCholesky(sortida$NoInversa, triangular = TRUE)[1:10,1:10]


genofile <- openfn.gds("data/test.gds")
gdsfmt::objdesp.gdsn(index.gdsn(genofile, "genotype") )$dim

read.gdsn( index.gdsn(genofile, "genotype"),
           start = c( 1, 1 ),
           count = c( 2133, 1000 ) )

closefn.gds(genofile)






*/