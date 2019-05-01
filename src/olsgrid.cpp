#include "include/olsgrid.h"


/* DE MOMENT AQUESTA FUNCIÓ NO S'UTILITZA */
int bdQR(Eigen::MatrixXd xl, Eigen::MatrixXd ql, Eigen::MatrixXd rtl)
{
   int inrows = xl.rows(), incols = xl.rows();
   int irank;
   
   Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(xl);
   Eigen::HouseholderQR<Eigen::MatrixXd> qr(xl);
   
   qr.compute(xl);
   irank = lu_decomp.rank();
   
   if (irank == inrows + 1 || irank == incols + 1 )
   {
      rtl = qr.matrixQR().template triangularView<Eigen::Upper>();
   }else {
      rtl = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>();  
   }
   // Només la part amb dades : ql = qr.householderQ().setLength(qr.nonzeroPivots()) 
   ql =  qr.householderQ();
   
   return(0);
}


// t: number of phenotypic traits in populations of n individuals
// tb: number of t elements in slabs. Com particionar de forma adequada t? quina és la longitud adequada???

// [[Rcpp::export]]
Rcpp::RObject ols_grid(Rcpp::RObject X, Rcpp::RObject y, int l, int r, int m, int t, int mb, int tb) //, Rcpp::RObject Yj, Rcpp::RObject Xr)
{
   char Lchar='L'; 
   char Nchar='N';
   char Tchar='T';
   char Uchar='U';
   
   int ione = 1;
   double done = 1.0;
   double dzero = 0.0;
   
   // 1a Part : 
   //    Codi : 
   //       Load XL  -> S'haurà de llegir del fitxer, de moment particiono la matriu que es passa per paràmetre
   //                   tenint en compte el valor de L (number of fixed covariates (sex, age, height, ....) )
   //       {QL,RTL} = qr(XL) 
   
      Eigen::MatrixXd XL = Rcpp::as<Eigen::MatrixXd>(X).block( 0, 0, Rcpp::as<Eigen::MatrixXd>(X).rows(), l);
      
      // NOTA : < La matriu XR en realitat no s'hauria de llegir aquí sinó que s'hauria de llegir per parts al pas 9 >
      //< La matriu XR està formada per m matrius de tamany mb >
      Eigen::MatrixXd XR = Rcpp::as<Eigen::MatrixXd>(X).block( 0, l, Rcpp::as<Eigen::MatrixXd>(X).rows(), Rcpp::as<Eigen::MatrixXd>(X).cols()-l);
      Eigen::MatrixXd QL, RTL;
      
      Rcpp::Rcout<<"\nMatrius XL i XR : \nXL\n"<<XL<<"\nXR\n"<<XR<<"\n";
      
      // bdQR(XL, QL, RTL);
      /* Es podria passar a funció  bdQR */
      int inrows = XL.rows(), incols = XL.cols();
      int irank;

      Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(XL);
      Eigen::HouseholderQR<Eigen::MatrixXd> qr(XL);
      // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(XL);
      
      qr.compute(XL);
      irank = lu_decomp.rank();
      
      if (irank == inrows + 1 || irank == incols + 1 ) {
         RTL = qr.matrixQR().template triangularView<Eigen::Upper>();
      }else {
         RTL = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>();
      }
      // QL =  qr.householderQ() * QL; // thin QL
      QL = qr.householderQ(); // QL Completa

   /* FI 1a part, es podria passar a funció  bdQR */
   
   Rcpp::Rcout<<"\nDescomposició QR de XL : \nQL\n"<<QL<<"\nRTL\n"<<RTL<<"\n";
   
   
   /* 2a part
         Codi (3-7): 
    
            for j:=1 to t/tb do
               Load slab Yˆj 
               Zˆ_T_j := Q^T_L * Yˆ_j                     (GEMM)
               Kˆ_j := R^−1_TL * Zˆ_T_j                   (TRSM)
            end for
   */
   
   Eigen::MatrixXd Y = Rcpp::as<Eigen::MatrixXd>(y); // Matriu n x t on hi ha emmagatzemats t marcadors phenotipics per n individus a la mostra
                                                      // Llegeixo tota la matriu passada per paràmetre però desprès s'haurà ade llegir
                                                      // només a trossos des del fitxer. S'haurà de llegir per blocs dins el bucle directament
                                                      // segons j i t/tb a Yj directament

   Rcpp::Rcout<<"\nMatriu Y : \n"<<Y<<"\n";
   
   int tbmaxsize = tb;
   Eigen::MatrixXd invRTL = RTL;
   int n = invRTL.rows(), info = 0;
   
   // dtrtri( char UPLO, char DIAG, int N, double A, int LDA, int INFO )	
   dtrtri_(&Uchar, &Nchar, &n, invRTL.data(), &n, &info ); // Càlcula la inversa de RTL i la funció retorna solució a invRTL directament
   
   Rcpp::Rcout<<"\nMatriu inversa RTL : \n"<<invRTL<<"\n";   
   
   for(int j=0; j<t/tb; j++) // Cada tb és un vector de Y
   {
        
      if( t < ((j+1)*tb) ) tbmaxsize = t - ( j * tb );
      
      int mm = QL.rows(), mk = QL.cols(), mn = tbmaxsize; 
      
      Rcpp::Rcout<<"\n tbmaxsize : "<<tbmaxsize<<"\n";   
      
      // Yj =  Rcpp::as<Eigen::MatrixXd>(Y).block(i*tb,tbmaxsize) ; // Load slab Yj
      Rcpp::Rcout<<"\nTamany matriu a llegir i posicions : (0, "<< j*tb<<" , "<<QL.cols()<<" , "<<tbmaxsize<<" )";
      Eigen::MatrixXd Yj = Y.block( 0, j*tb, QL.cols(), tbmaxsize );

      Eigen::MatrixXd ZTj (mm, mn);
      
      Rcpp::Rcout<<"\n(4) Bloc llegit matriu Yj : \n"<<Yj<<"\n amb "<<mn<<" vectors\n";

      // get ZT_j
      // dgemm_( char TRANSA, char TRANSB, int M, int N, int K, double ALPHA, double* A, int LDA, double* B, int LDB, double BETA, double* C, int LDC )	
      dgemm_(&Tchar, &Nchar, &mk, &mn, &mk, &done, QL.data(), &mk, Yj.data(), &mm, &dzero, ZTj.data(), &mk);   // mk enlloc mm perquè transposa(QL)
      
      Rcpp::Rcout<<"\n (5) dgemm càlcul Ztj : \n"<<ZTj<<"\n";
      
      Eigen::MatrixXd Kj = ZTj;
      mm = Kj.rows(), mn = Kj.cols(), mk = invRTL.rows();
      
      if( mk < mm ) {
         int insert = mm-mk;
         invRTL.conservativeResizeLike( Eigen::MatrixXd::Identity(invRTL.rows() + insert, invRTL.cols()+insert));
      }
      
      // get Kˆj
      // dtrsm( char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double A, int LDA, double B, int LDB )
      dtrsm_(&Lchar, &Uchar, &Nchar, &Nchar, &mm, &mn, &done, invRTL.data(), &mm, Kj.data(), &mm ); // La solució (X) es sobrescriu a (B) --> B es converteix en X

      Rcpp::Rcout<<"\n (6) trsm càlcul Kj : \n"<<Kj<<"\n";
   }
   /* FI 2a part, es podria passar a funció  .... */
   
   
   //..// ----------------------------------  FINS AQUÍ OK !!!!  ------------------------------------ //..//
   
   
   
   /* 3a part :  Codi (8-25): 
    
            for i:=1 to m/mb do
               Load slab XˆRi
               Rˆ_TRi : = QTL * Xˆ_Ri			            (GEMM)
               Tˆi :=Xˆ_Ri - Q_L * Rˆ _TRi 	            (GEMM)
               Hˆi :=Rˆ-1_TL * RˆTRi 			            (TRSM)
               for k:=1 to mb do
                  {QRik , RBRik } := qr(Tik )            (QR)
               end for
               for j := 1 to t/tb do
                  Load slab Yˆj
                  Zˆ_Bij : = QˆT_Ri * Yˆj                (GEMM)
                  for k:=1 to mb do
                     Bˆ_Bikj : = R^-1_BRik * Qˆ_Bikj     (TRSM)
                     B_Tikj := Kˆ_j -  Hˆ_ik * Bˆ_Bikj   (GEMM)
                  end for
                  Store slab Bˆij
               end for
            end for
   */
   
   
   // m : Nombre de matrius
   // mb : matriu de m
   // slab XRi : bloc de mb matrius de m
   
   int mbmaxsize = mb;
   
   Rcpp::Rcout<<"\nm i mb valen : "<<m<<" - "<<mb<<"\n";
   int xrc = XR.cols(), xrr = XR.rows();
   
   for(int i=0; i<m/mb; i++)  // mb(number of matrix)
   {
      int mm = QL.cols(), mk = QL.rows();
      
      if( m < ((i+1)*mb) ) mbmaxsize = m - ( i * mb );
      
      // Yj =  Rcpp::as<Eigen::MatrixXd>(Y).block(i*tb,tbmaxsize) ; // Load slab Yj
      Eigen::MatrixXd XRi = XR.block( 0, i*mb, xrr, mbmaxsize);
      
      int mn = XRi.cols();
      
      Rcpp::Rcout<<"\n(9) Bloc llegit matriu XRi : \n"<<XRi<<"\n Amb columnes = "<<mbmaxsize<<"\n";   
      
      
      
      Eigen::MatrixXd RTRi (mm, mn);
   
      // get RTRi
      dgemm_(&Tchar, &Nchar, &mm, &mn, &mk, &done, QL.data(), &mm, XRi.data(), &mk, &dzero, RTRi.data(), &mm);
      
      Rcpp::Rcout<<"\n (10) dgemm càlcul RTRi : \n"<<RTRi<<"\n";
      
      mm = QL.rows(); mk = QL.cols(); mn = RTRi.cols();
      Eigen::MatrixXd QL_RTRi (mm, mn);

      // get Ti
      dgemm_(&Nchar, &Nchar, &mm, &mn, &mk, &done, QL.data(), &mm, RTRi.data(), &mk, &dzero, QL_RTRi.data(), &mm);   
      Eigen::MatrixXd Ti = XRi - QL_RTRi;
      
      Rcpp::Rcout<<"\n (11*) dgemm càlcul QL_RTRi : \n"<<QL_RTRi<<"\n";
      Rcpp::Rcout<<"\n (11) dgemm càlcul Ti = XRi - QL_RTRi : \n"<<Ti<<"\n";

      Eigen::MatrixXd Hi = RTRi;
      // get Hi
      dtrsm_(&Lchar, &Uchar, &Nchar, &Nchar, &mm, &mn, &done, invRTL.data(), &mm, Hi.data(), &mm ); // La solució (X) es sobrescriu a (B) --> B es convertei
      
      Rcpp::Rcout<<"\n (12) dgemm càlcul Hi : \n"<<Hi<<"\n";
      
      Eigen::MatrixXd QRik, RBRik;
      int mbprov = Ti.cols();
      
      Eigen::MatrixXd QR (xrc, xrc);  // XR.cols() * Xr.cols()
      
      for( int k=0; k<mbprov; k++ )
      {
         // COM HE DE PARTICIONAR Ti ???
         Eigen::VectorXd Tik = Ti.col(k);
         Eigen::FullPivLU<Eigen::MatrixXd>lu_decomp(Tik);
         Eigen::HouseholderQR<Eigen::MatrixXd> qr(Tik);
         
         qr.compute(Tik);
         irank = lu_decomp.rank();
         
         if (irank == inrows + 1 || irank == incols + 1 ) {
            RBRik = qr.matrixQR().template triangularView<Eigen::Upper>();
         }else {
            RBRik = qr.matrixQR().topLeftCorner(irank, irank).template triangularView<Eigen::Upper>();  
         }
         // Només la part amb dades : ql = qr.householderQ().setLength(qr.nonzeroPivots()) 
         QRik =  qr.householderQ();
         
         QR.block(  0, i*mb, xrr, mbmaxsize) = QRik;
         
         Rcpp::Rcout<<"\n (14) Descomposició Tik : \n"<<Tik<<"\nen QRik : \n"<<QRik<<"\n RBRik : \b"<<RBRik<<"\n";
         Rcpp::Rcout<<"\n La matriu QR es va construint : \n"<<QR<<"n";

      }
      /**** AQUÍ !!! MIRANT COM FER EL DE LA LÍNIA VERDA !!, 
       * 
       *COM MONTAR LA MATRIU QRi A PARTIR DE les operacions QRik ???
       si es pot montar la matriu ja estarà soluconat però falta saber com fer-ho !!!!
       
       -- Revisar si hi ha informació a internet del càlcul QR per blocs !!!!
      */
      /*
      for(int j=0; j<t/tb; j++)
      {
         
         if( t < ((j+1)*tb) ) tbmaxsize = t - ( j * tb );
         
         // Yj =  Rcpp::as<Eigen::MatrixXd>(Y).block(i*tb,tbmaxsize) ; // Load slab Yj
         Rcpp::Rcout<<"\nTamany matriu a llegir i posicions : (0, "<< j*tb<<" , "<<QL.cols()<<" , "<<tbmaxsize<<" )";
         Eigen::MatrixXd Yj = Y.block( 0, j*tb, QRik.cols(), tbmaxsize );
      }
       */
      
   }
   
   
   return Rcpp::List::create(Rcpp::Named("QL") = QL,
                             Rcpp::Named("RTL") = RTL
                           );
   
   
}



// [[Rcpp::export]]
Rcpp::RObject matv_gemv(Rcpp::RObject M, Rcpp::RObject V) //, Rcpp::RObject Yj, Rcpp::RObject Xr)
{
   
   Eigen::MatrixXd mat = Rcpp::as<Eigen::MatrixXd> (M);
   Eigen::VectorXd vect = Rcpp::as<Eigen::VectorXd> (V);
   int imrows = mat.rows(), imcols = mat.cols();
   Eigen::VectorXd result(imrows);

   char Nchar='N';
   char Tchar='T';
   int ione = 1;
   double done = 1.0;
   double dzero = 0.0;

   dgemv_(&Nchar, &imrows, &imcols, &done, mat.data(), &imrows, vect.data(), &ione, &dzero, result.data(), &ione);
                      
   return(Rcpp::wrap(result));
   
}


// [[Rcpp::export]]
Rcpp::RObject matv_gemm(Rcpp::RObject a, Rcpp::RObject b) //, Rcpp::RObject Yj, Rcpp::RObject Xr)
{
   
   Eigen::MatrixXd A = Rcpp::as<Eigen::MatrixXd> (a);
   Eigen::MatrixXd B = Rcpp::as<Eigen::MatrixXd> (b);
   int m = A.rows(), k = A.cols(), n = B.cols();
   
   Eigen::MatrixXd C(m, n);
   
   char Nchar='N';
   int ione = 1;
   double done = 1.0, dzero = 0.0;
   
   //dgemm(char TRANSA,char TRANSB,int M,int N,int K,double ALPHA,double A,int LDA,double B,int LDB,double BETA,double C,int LDC)
   
   dgemm_(&Nchar, &Nchar, &m, &n, &k, &done, A.data(), &m, B.data(), &k, &dzero, C.data(), &m);
   
   return(Rcpp::wrap(C));
   
}


// [[Rcpp::export]]
Rcpp::RObject matv_dtrsv(Rcpp::RObject M, Rcpp::RObject v) 
{
   
   Eigen::MatrixXd mat = Rcpp::as<Eigen::MatrixXd> (M);
   Eigen::VectorXd vect = Rcpp::as<Eigen::VectorXd> (v);
   int imrows = mat.rows(), imcols = mat.cols();
   
   char Nchar='N';
   char Tchar='T'; 
   char Uchar='U'; //Nova def.
   int ione = 1;
   double done = 1.0;
   double dzero = 0.0;

   // dtrsv	( char UPLO, char TRANS, char DIAG, int N, double precision dimension(lda,*) A, int LDA, double precision dimension(*) X, int INCX )
   dtrsv_(&Uchar, &Nchar, &Nchar, &imcols, mat.data(), &imrows, vect.data(), &ione); // La solució es sobrescriu a vect - Ax=b --> b es converteix en x

   return(Rcpp::wrap(vect));

}


// [[Rcpp::export]]
Rcpp::RObject matv_dtrsm(Rcpp::RObject R, Rcpp::RObject Z) 
{
   char Nchar='N';
   char Rchar='R';
   char Tchar='T';
   char Uchar='U';
   char Lchar='L'; //Nova def.
   int ione = 1;
   double done = 1.0;
   double dzero = 0.0;
   
   Eigen::MatrixXd A = Rcpp::as<Eigen::MatrixXd> (R);
   Eigen::VectorXd B = Rcpp::as<Eigen::VectorXd> (Z);
   int m = B.rows(), n = B.cols(), k = A.rows();
   
   if( k < m ) {
      int insert = m-k;
      A.conservativeResizeLike( Eigen::MatrixXd::Identity(A.rows() + insert, A.cols()+insert));
   }
   
   
   // dtrsm( char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double A, int LDA, double B, int LDB )
   dtrsm_(&Lchar, &Uchar, &Nchar, &Nchar, &m, &n, &done, A.data(), &m, B.data(), &m ); // La solució (X) es sobrescriu a (B) --> B es converteix en X
   
   return(Rcpp::wrap(B));
   
}

// [[Rcpp::export]]
Rcpp::RObject matm_dtrtri(Rcpp::RObject R )
{
   char Nchar='N';
   char Uchar='U';

   Eigen::MatrixXd A = Rcpp::as<Eigen::MatrixXd> (R);
   int n = A.cols();
   int info=0;

   // dtrtri( char UPLO, char DIAG, int N, double A, int LDA, int INFO )	
   dtrtri_(&Uchar, &Nchar, &n, A.data(), &n, &info ); // La solució sobrescriu (A)
   
   Rcpp::Rcout<<"\nInfo val : "<<info<<"\n";
   
   return(Rcpp::wrap(A));
   
}





/***R
library(microbenchmark)

# Proves olds-grid
Xl <- matrix(c(1.0,-1.0, 4.0, 1.0, 4.0, -2.0, 1.0, 4.0, 2.0, 1.0, -1.0, 0.0), byrow=TRUE, nrow = 4)
Xr <- matrix(c(1.0, 4.0, 4.0, 1.0, 2.0, -1.0, 0.0, 3.0), byrow=TRUE, nrow = 4)

Xm <- matrix(c(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0,
               1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,
               1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0), byrow = TRUE, nrow = 4)

X <- cbind(Xl,Xr,Xm) ; X

Yt <- matrix(c(1.0, 0.0, 1.0, 1.0, 0.0, 1.0,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
               0.0, 1.0, 1.0, 0.0, 1.0, 1.0,
               1.0, 0.0, 1.0, 1.0, 0.0, 1.0), byrow = TRUE, nrow = 4)



Y <- Yt; Y

l <- ncol(Xl) # Part fixa
r <- ncol(Xr)
n <- nrow(X)   # individus
Y <- cbind(Yt)  # Matriu trets fenotípics
t <- ncol(Yt)  # Trets fenotípics
m <- r + ncol(Xm) # Grups marcadors + part variable

tb <- 1
mb <- 2


# m groups of genetic markers and t phenotypic traits in populations of n individuals.
# each trait is represented by a vector yj containing the trait measurements (one per individual)
# each matrix Xi = [XL|XRi ] is composed of a set of l fixed covariates such as sex, age,... (XL)
# and one of the groups of r = p -l XRi

# n individuals : numero de files
# m grups de marcadors genètic : m*nombre marcadors columnes ?
# t phenotypic : t columnes adicionals ?
# p : nombre total columnes
# l : nombre fixat de covariables.
# r : nombre variable de covariables

# ols_grid(X, Y, l, r, m, t, mb, tb) 
ols_grid(X, Y, l, r, m, t, mb, tb) 

# Comprova funcionament primera descomposició QR (linea 2 (previ primer for))
qrstr <- qr(Xl)   # dim(x) == c(n,p)
qr.Q(qrstr)
qr.R(qrstr)







## A PARTIR D'AQUÍ COMPROVACIONS DIFERENTS FUNCION LAPACK-BLAS

# Proves dtrtri_
## upper triangular matrix 'r':
r <- matrix(c(1,2,3,0,1,1,0,0,2), byrow = TRUE, nrow = 3)

matm_dtrtri(r)



# Proves drstm_
## upper triangular matrix 'r':
r <- rbind(c(1,2,3),
           c(0,1,1),
           c(0,0,2))
x <- c(8,4,2,4)

backsolve(r, x) # 8 -12 -5
matv_dtrsm(r,x)

results <- microbenchmark(rsol <- backsolve(r, x),
                          dtrsm <- matv_dtrsm(r,x),
                          times = 2L)  # Proves multiplicacions x blocs
print(summary(results)[, c(1:7)],digits=3)
stopifnot(all.equal(rsol, dtrsm))

# Proves dgemm_
A <- matrix(c(1.0,-1.0, 4.0, 1.0, 4.0, -2.0, 1.0, 4.0, 2.0, 1.0, -1.0, 0.0), byrow=TRUE, nrow = 4)
B <- matrix(c(1,2,3,0,1,1,0,0,2), byrow=TRUE, nrow=3)


n <- 1000; p <- 500
A <- matrix(rnorm(n*p), nrow=n, ncol=p)

n <- 500 ; p <- 1000
B <- matrix(rnorm(n*p), nrow=n, ncol=p)

dgemm <- matv_gemm(A,B);dgemm
A%*%B

results <- microbenchmark(rsol <- A%*%B,
                          dgemm <- matv_gemm(A,B),
                          times = 2L)  # Proves multiplicacions x blocs

print(summary(results)[, c(1:7)],digits=3)

stopifnot(all.equal(rsol, dgemm))


# Proves dtsrv_

A <- rbind(c(1,2,3),c(0,1,1),c(0,0,2))
b <- c(8,4,2)
y <- backsolve(A, b) ; y


## x1 x2 x3 
##  2  3 -1
sol <- matv_dtrsv(A,b)
sol

# Proves dgemv_
X <- matrix(c(1.0,-1.0, 4.0, 1.0, 4.0, -2.0, 1.0, 4.0, 2.0, 1.0, -1.0, 0.0), byrow=TRUE, nrow = 4)
Y <- matrix(c(3.0, 2.0, 4.0), byrow=TRUE, nrow = 3)
X

dec <- ols_grid( X, Y, 2)
dec$Q
dec$R

results <- microbenchmark(rsol <- X%*%Y,
                          blastsol <- matv_gemv(X,Y),
                          times = 2L)  # Proves multiplicacions x blocs

print(summary(results)[, c(1:7)],digits=3)





library(microbenchmark)
library(DelayedArray)
test_LAPACK()

X <- matrix(c(-1.0,-8.0,0.0,-1.0,1.0,-5.0,3.0,0.0,2.0), byrow=TRUE, nrow = 3)
n <- 100
X <- matrix(rnorm(n*n), nrow=n, ncol=n)

results <- microbenchmark(svdLP1 <- test_LAPACK_eig(X),
                          svd <- svd(X),
                          times = 2L)  # Proves multiplicacions x blocs

print(summary(results)[, c(1:7)],digits=3)
svdLP1$d[1:10]
svd$d[1:10]




Q2 <- cbind(qr.Q(qrstr), c(-0.5, -0.5, 0.5, 0.5))
# yj <- matrix(c( 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0),byrow=TRUE, nrow = 4)
yj <- matrix(c( 1.0, 1.0, 0.0, 1.0),byrow=TRUE, nrow = 4)
yj <- matrix(c( 0.0, 1.0, 1.0, 0.0),byrow=TRUE, nrow = 4) # Segona columna
t(Q2)%*%yj

# En aquesta comprovació !!!
invrtl <- solve(qr.R(qrstr))
library(matlib)

backsolve(invrtl, t(Q2)%*%yj) # Això és el que hauria de donar però que no donar Kj
matv_dtrsm(invrtl,t(Q2)%*%yj)
matv_dtrsv(invrtl,t(Q2)%*%yj)

qr.R(qrstr)
solve(qr.R(qrstr))

invR <- inv(qr.R(qrstr))

XRi <-  matrix(c(1,4,4,1,2,-1,0,3), byrow = TRUE, nrow = 4)
RTRi <- matrix(c(-3.5,-3.5,-2.5,3.5,0.5,0.5,-1.5,-1.5), byrow = TRUE, nrow = 4)
t(Q2)%*%RTRi
XRi - Q2%*%RTRi

invR <- inv(qr.R(qrstr))
invR <- cbind(invR, c(0.0, 0.0, 0.0))
invR <- rbind(invR, c(0.0, 0.0, 0.0,1.0))

backsolve(invR, RTRi)
invR%*%RTRi



*/