#ifndef BIGDATASTATMETH_UTILITIES_HPP
#define BIGDATASTATMETH_UTILITIES_HPP

#include <RcppEigen.h>


namespace BigDataStatMeth {

    // Structs
    struct fullpath {
        std::string path ;
        std::string filename ;
    };


    struct svdeig {
        Eigen::VectorXd d;
        Eigen::MatrixXd u;
        Eigen::MatrixXd v;
        bool bokuv = false;
        bool bokd = false;
        std::string hdf5file = "";
    };



    // Functions
    
    // Get path and filename from full route
    extern inline fullpath SplitElementName (std::string str)
    {
        fullpath currentpath;
        std::size_t found = str.find_last_of("/\\");
        
        if( found< str.length() ) {
            currentpath.filename =  str.substr(found+1);
            currentpath.path = str.substr(0,found);
        }else {
            currentpath.filename = str;
            currentpath.path = "";
        }
        return(currentpath);
    }
    
    
    // Get RObject data type
    extern inline std::string getObjecDataType(Rcpp::RObject obj) 
    {
        
        std::string strtype = "";
        
        try 
        {
            if( Rcpp::is<Rcpp::NumericVector>(obj) ) {
                strtype = "numeric";
            } else if( Rcpp::is<Rcpp::IntegerVector>(obj) ) {
                strtype = "int";
            } else if( Rcpp::is<Rcpp::CharacterVector>(obj) ) {
                strtype = "char";
            } else if( Rcpp::is<Rcpp::LogicalVector>(obj) ) {
                strtype = "logic";
            } else if( Rcpp::is<Rcpp::DataFrame>(obj) ) {
                strtype = "dataframe";
            } else if( Rcpp::is<Rcpp::List>(obj) ) {
                strtype = "list";
            } else if( obj.isS4() ) {
                strtype = "S4";
            } else if( obj.isNULL() ) {
                strtype = "NULL";
            } else {
                strtype = "unknown";
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getObjecDataType: "<< ex.what() << "\n";
        }
        
        return(strtype);
    }
    
    
    
    // Get RObject dims
    extern inline Rcpp::IntegerVector getObjectDims(Rcpp::RObject obj, std::string strtype) 
    {
        
        Rcpp::IntegerVector dims(2);
        
        try 
        {
            if(strtype =="") {
                strtype = getObjecDataType(obj);    
            }
            
            if( strtype == "numeric"  || strtype == "int" || strtype == "factor" ){
                if( Rf_isMatrix(obj)) {
                    dims[0] = Rcpp::as<Rcpp::NumericMatrix>(obj).rows();
                    dims[1] = Rcpp::as<Rcpp::NumericMatrix>(obj).cols();
                } else {
                    dims[0] = 1;
                    dims[1] = Rcpp::as<Rcpp::NumericVector>(obj).length();
                }
            }else if( strtype == "logic" ) {
                if( Rf_isMatrix(obj)) {
                    dims[0] = Rcpp::as<Rcpp::LogicalMatrix>(obj).rows();
                    dims[1] = Rcpp::as<Rcpp::LogicalMatrix>(obj).cols();
                } else {
                    dims[0] = 1;
                    dims[1] = Rcpp::as<Rcpp::LogicalVector>(obj).length();
                }
            } else if( strtype == "char" ){
                if( Rf_isMatrix(obj)) {
                    dims[0] = Rcpp::as<Rcpp::CharacterMatrix>(obj).rows();
                    dims[1] = Rcpp::as<Rcpp::CharacterMatrix>(obj).cols();
                } else {
                    dims[0] = 1;
                    dims[1] = Rcpp::as<Rcpp::CharacterVector>(obj).length();
                }
            } else if(strtype == "dataframe"){
                dims[0] = Rcpp::as<Rcpp::DataFrame>(obj).nrows();
                dims[1] = Rcpp::as<Rcpp::DataFrame>(obj).length();
                
            } else {
                dims[0] = 0;
                dims[1] = 0;
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getObjecDataType: "<<ex.what()<< " \n";
        }
        
        return(dims);
    }
    
    
    
    // Get RObject dims
    extern inline int getMaxBlockSize ( int nRowsA, int nColsA, int nRowsB, int nColsB, Rcpp::Nullable<int> block_size = R_NilValue) 
    {
        
        int iblock_size;
        try
        {
            iblock_size = std::min( std::min( nRowsA, nColsA), std::min( nRowsB, nColsB) );
            
            if (block_size.isNotNull()) {
                if( Rcpp::as<int> (block_size) < iblock_size ){   
                    iblock_size = Rcpp::as<int> (block_size); }
            } else {
                //..// iblock_size = std::min(  std::min(dsA->nrows(),dsA->ncols()),  std::min(dsB->nrows(), dsB->ncols()));
                if (iblock_size>1024)
                    iblock_size = 1024;
            }
            
                
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getObjecDataType: "<<ex.what()<< " \n";
        }
        
        return(iblock_size);
    }
    
    
    extern inline size_t getOptimBlockSize( size_t fullSize, size_t blockSize, size_t iDesp, size_t currentSize ) 
    {
        
        try
        {
            if( iDesp + blockSize == fullSize - 1) {
                currentSize = blockSize + 1;
            } else if( iDesp + blockSize > fullSize ) { 
                currentSize = fullSize - iDesp; 
            } else {
                currentSize = blockSize;
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getObjecDataType: "<<ex.what()<< " \n";
        }
        
        return(currentSize);
    }
    
}

#endif // BIGDATASTATMETH_UTILITIES_HPP