#include "include/pca_hdf5.h"



Eigen::VectorXd cumsum_hdf5(Eigen::VectorXd x){
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


// Get's variance and cumulative variance from svd decomposition
// and write results to hdf5 file
int get_HDF5_PCA_var_ptr(  H5File* file, std::string strdataset)
{
  try
  {
    
    Exception::dontPrint();
    
    DataSet dataset = file->openDataSet(strdataset);
    
    // Hem de llegir el vector lambda (d)
    // Un cop llegit hem de calcular : 
    //        var <- res$d^2/sum(res$d^2)
    //        pervar <- cumsum(var) 
    //        
    // Escriure els resultats al fitxer
    
    /***
     * FALTA COMPROVAR COM S'ESCRIUEN LES MATRIU V I D AL FITXER !!! SEMBLA QUE V NO S'ESCRIVIA GENS BÉ !!!!
     * La modificació que he fet ha estat crear una funció que en principi hauria d'escriure els valors per files enlloc de per
     * columnes.... bé, la veritat es que només l'he creat i he canviat les dimensions, faltarà veure qeu tal funciona !!!! '
     */
    
  }catch( FileIException error ) {
    error.printErrorStack();
    return -1;
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    error.printErrorStack();
    return -1;
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    error.printErrorStack();
    return -1;
  } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
    error.printErrorStack();
    return -1;
  }
  
  
}



// 
//' PCA Descomposition
//' 
//' Compute PCA
//' 
//' @param matrix or DelayedArray 
//' @return 
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdPCA_hdf5(std::string filename, std::string strsubgroup, std::string strdataset)
{
  try
  {
    
    Eigen::MatrixXd X;
    Eigen::MatrixXd U, V, Lambda;
    Eigen::VectorXd lambda;
    Eigen::MatrixXd C, D, varcoord;
    
    pcaeig pcaRes;
    
    // Open an existing file and dataset.
    H5File* file = new H5File( filename, H5F_ACC_TRUNC );
    //..// 
    
    // Look up for svd decomposition in file
    if( exists_HDF5_element_ptr(file,"/SVD/" + strdataset )!=0)
    {
      Rcpp::Rcout<<"\n Doncs ara ens toca calcular .... Fatal !!!!\n";
      DataSet dataset = file->openDataSet(strsubgroup + "/" + strdataset);
      svdeig retsvd = RcppbdSVD_hdf5( filename, strsubgroup, strdataset, 4, 1, 0, true, true );
      dataset.close();
    }
    
    //..// DataSet dataset = file->openDataSet("/SVD/"+strdataset+"/d");
    
    get_HDF5_PCA_var_ptr(file,"/SVD/"+strdataset+"/d" );
    
    
      
      
      //..// ans <- list(varcoord=var.coord, Y=Y, var = svdX$d^2, percvar=variance, components = svdX$v, method = method)
    
    /****

      Rcpp::Rcout<<"\n Descomposició +  assignació  OK : \n";
      
      {
        Eigen::VectorXd prov = lambda.array().sqrt();
        D = prov.asDiagonal();
      }
      
      Rcpp::Rcout<<"\n D assignada  OK : \n";
      
      //..// Lambda = lambda.asDiagonal();
      
      
      // Cumulative variance
      //..// Eigen::VectorXd pvac = cumsum(lambda/lambda.array().sum()); // percentatge de vari`ancia acumulada pels components pvac
      
      // Get correlations
      {
        Rcpp::Rcout<<"\n Calculem produc per a C....: \n";
        Eigen::MatrixXd mprov = block_matrix_mul_parallel( Eigen::MatrixXd::Identity(V.rows(), V.cols()), V, 512 );
        Rcpp::Rcout<<"\n Product OK : \n";
        C = block_matrix_mul_parallel( mprov, D, 512 );
      }
      
      Rcpp::Rcout<<"\n Correlacions realitzades -  OK : \n";
      
      // Get coordinates
      {
        svdeig singX = RcppbdSVD(X, k, int(), false);
        if(singX.bokuv == true && singX.bokd==true)
        {
          varcoord = block_matrix_mul_parallel(X.adjoint(), singX.u, 512 );  
        } else {
          throw(Rcpp::exception("Error with svd decomposition","pca.cpp",4));
        }
      }
      Rcpp::Rcout<<"\n Càlcul coordenades realitzades -  OK : \n";
      
      Lambda = lambda.asDiagonal(); 
    

    
    
    return Rcpp::List::create(Rcpp::Named("U") = U,
                              Rcpp::Named("V") = V,
                              Rcpp::Named("D") = D,
                              Rcpp::Named("Lambda") = Lambda,
                              Rcpp::Named("lambda") = lambda,
                              Rcpp::Named("pvac") = cumsum(lambda/lambda.array().sum()), // percentatge de vari`ancia acumulada pels components pvac,
                              Rcpp::Named("var.contr") = V.pow(2), // Variable contrib.
                              Rcpp::Named("C") = C,
                              Rcpp::Named("var.coord") = varcoord,
                              Rcpp::Named("var.cos2") = C.pow(2)
    );
    ***/
    return Rcpp::List::create(Rcpp::Named("U") = "U");
    
    
    
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  
}



/*** R

library(BigDataStatMeth)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/UVic/TFM2 - BigDatStatMeth_Mgcca/Analisi/")
bdPCA_hdf5("tmp/testgdsfile.h5","invs","geno")










library(FactoMiner)

data(decathlon, package="FactoMineR")
data<-decathlon

temps <- data[c(1,5,6,10)]
maxt <- apply(temps,2,max)
transf <- -sweep(temps,2,maxt,FUN="-")
data[c(1,5,6,10)] <- transf

basic<-data[-(39:41),-(11:13)]
n <- nrow(basic)

X <- BigDataStatMeth::Normalize_Data(as.matrix(basic)) # OK
X <- X*(1/sqrt(n-1))  # Tipifiquem ¿¿Cal fer-ho sempre?? hi ha gent que ho fa i d'altres no.

# Anàlisi de components :

# Calculem matriu de covariances : Crossprod
covX <- BigDataStatMeth::bdcrossprod(X)
sing <- BigDataStatMeth::bdSVD(covX)      # Fins aquí ho dono OK (hi ha una petita dif que entenc que és perquè tenim k-1 components)


U <- sing$u
V <- sing$v
D <- sqrt(diag(sing$d))
Lambda <- D^2
lambda <- diag(Lambda)  # Es equivalent a sing$d
lambda               # FINS AQUÍ HAURIA DE SER CORRECTE PERQUÈ LES FUNCIONS FUNCIONEN BÉ, 
# LES DADES NO SON CORRECTES PERQUÈ TREBALLEM AMB N-1 COMPONENTS I EN AQUEST CAS N'HI HA MOLT
# POQUES, AIXÒ CAUSA DIFERÈNCIES IMPORTANTS EN AQUEST CAS !!!!!

pvac <- cumsum(lambda/sum(lambda)) # percentatge de vari`ancia acumulada pels components pvac
pvac

var.contr <- V^2*100 # en %
pvac


# funció programada BDSM
res <- bdPCA(as.matrix(basic), 500)

res$var.contr
res$pvac
res$var.coord

rownames(res$var.coord) <- colnames(basic)[1:nrow(res$var.coord)]
par(cex=.7)
plot(res$var.coord[,1:2],xlim=c(-1,1),ylim=c(-1,1),
     type="n",xlab="PC1 (33%)",ylab="PC2 (15.5%)") # type="n" -no punts- 
abline(h=0); abline(v=0)
noms<-abbreviate(rownames(res$var.coord))
text(res$var.coord[,1:2],labels=noms) 
arrows(rep(0,10),rep(0,10),res$var.coord[,1], res$var.coord[,2],col="red") 
curve(sqrt(1-x^2),-1,1,add=T)
curve(-sqrt(1-x^2),-1,1,add=T)



par(cex=.7)
plot(res$var.coord[,1:2],xlim=c(-4,6),ylim=c(-4,4),
     xlab="Dim1 (33,04%)",ylab="Dim2 (15.47%)", pch=19)
abline(h=0); abline(v=0)
noms<-abbreviate(rownames(res$var.coord))
text(res$var.coord[,1:2]+0.1,labels=noms)





############# GRAFIQUILLES EVITABLES
###     Gràfica amb valors propis
plot(lambda, xlab="Ordre valors propis", ylab="Valor propi",
     xaxt="n", yaxt="n", col.lab="blue", cex.lab=0.8)
axis(1, at = seq(1, 10, by = 1), las=1, cex.axis=0.6)
axis(2, at = seq(0, 3, by = 0.5), las=1, cex.axis=0.6)
lines(lambda,lwd=2)
abline(h=1,col="red",lty=2)
###     Gràfica variància explicada acumulada
plot(pvac*100, xlab="Ordre valors propis", ylab="Variància acumulada (%)", xaxt="n", yaxt="n", col.lab="blue", cex.lab=0.8)
axis(1, at = seq(1, 10, by = 1), las=1, cex.axis=0.6)
axis(2, at = seq(0, 100, by = 10), las=1, cex.axis=0.6)
lines(pvac*100,lwd=2)
abline(h=70,col="red",lty=2)
########## FI GRAFIQUILLES EVITABLES

## FINS AQUÍ TOT OK

## analisi per columnes
# contrib de variables en comp
var.contr <- V^2*100 # en %
rownames(var.contr) <- colnames(basic) colnames(var.contr)<-paste("PC",1:ncol(var.contr),sep="") round(var.contr,2)







Xr[1:10,1:10]
X[1:10,1:10]

### Fet seguint els apunts per comprovar
mitjanesbasiques <- apply(basic,2,mean)
sdbasiques <- apply(basic,2,sd)

Xr <- scale(basic,center = TRUE,scale = sdbasiques)
Xr <- Xr*(1/sqrt(n-1)) 
singr <- svd(t(Xr) %*% Xr)

Ur <- singr$u
Vr <- singr$v
Dr <- sqrt(diag(singr$d))
Lambdar <- Dr^2
lambdar <- diag(Lambdar)  # Es equivalent a sing$d
lambdar
pvacr <- cumsum(lambdar/sum(lambdar))


res$pvac
pvacr

# Grafica amb Valors propis
plot(lambdar, xlab="Ordre valors propis", ylab="Valor propi",
     xaxt="n", yaxt="n", col.lab="blue", cex.lab=0.8)
axis(1, at = seq(1, 10, by = 1), las=1, cex.axis=0.6)
axis(2, at = seq(0, 3, by = 0.5), las=1, cex.axis=0.6)
lines(lambdar,lwd=2)
abline(h=1,col="red",lty=2)



par(cex=.7)
plot(res$var.coord[,1:2],xlim=c(-1,1),ylim=c(-1,1),
     type="n",xlab="PC1 (33%)",ylab="PC2 (15.5%)") # type="n" -no punts- 
abline(h=0); abline(v=0)
noms<-abbreviate(rownames(res$var.coord))
text(res$var.coord[,1:2],labels=noms) 
arrows(rep(0,10),rep(0,10),res$var.coord[,1], res$var.coord[,2],col="red") 
curve(sqrt(1-x^2),-1,1,add=T)
curve(-sqrt(1-x^2),-1,1,add=T)

*/
