#ifndef BIGDATASTATMETH_UTILITIES_HPP
#define BIGDATASTATMETH_UTILITIES_HPP

// #include <RcppEigen.h>
#include <BigDataStatMeth.hpp>


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

    struct strQR {
        Eigen::MatrixXd Q;
        Eigen::MatrixXd R;
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
    
    
    
    // Get maximum block size 
    extern inline int getMaxBlockSize ( int nRowsA, int nColsA, int nRowsB, int nColsB, int ifactor, Rcpp::Nullable<int> block_size = R_NilValue) 
    {
        
        int iblock_size;
        
        try
        {
            
            iblock_size = std::min( std::min( nRowsA, nColsA), std::min( nRowsB, nColsB) );

            if (block_size.isNotNull()) {
                // if( Rcpp::as<int> (block_size) < iblock_size ) {
                //     iblock_size = Rcpp::as<int> (block_size); }
                iblock_size = Rcpp::as<int> (block_size);
                if( iblock_size > (MAXBLOCKSIZE / ifactor) ) {
                    Rcpp::warning("Warning: block size %i is bigger than the maximum recomended %i.", iblock_size, (MAXBLOCKSIZE / ifactor));
                }
            } else {
                //..// iblock_size = std::min(  std::min(dsA->nrows(),dsA->ncols()),  std::min(dsB->nrows(), dsB->ncols()));
                if (iblock_size > (MAXBLOCKSIZE / ifactor))
                    iblock_size = MAXBLOCKSIZE / ifactor;
            }
                
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getObjecDataType: "<<ex.what()<< " \n";
        }
        
        return(iblock_size);
    }
    
    
    // Compute the maximum block size tacking in to accout the last column / row to
    // avoiding single rows or columns and perform an extra step
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
            Rcpp::Rcout<< "c++ exception getOptimBlockSize: "<<ex.what()<< " \n";
        }
        
        return(currentSize);
    }
    
    
    
    // Get the number of rows to read taking in to account the maximum elements per block
    // util when we have rectangular matrices, especially in omics data where we have
    // few samples and thousands of variables
    extern inline std::vector<hsize_t> getMatrixBlockSize( int nrows, int ncols ) 
    {
        size_t  maxRows = nrows,
                maxCols = ncols;
        
        std::vector<hsize_t> blockSize = {0, 0};
        
        try
        {
            // Calculem el mínim de files
            if( nrows < ncols ) {
                if( maxRows < MAXBLOCKSIZE ){
                    maxRows = nrows;
                } else{
                    maxRows = MAXBLOCKSIZE;
                }
                
                maxCols = std::floor( MAXELEMSINBLOCK / maxRows );
                if( maxCols> ncols || maxCols + 1 == ncols) {
                    maxCols = ncols;
                }
            } else {
                if( maxCols < MAXBLOCKSIZE ){
                    maxCols = ncols;
                } else{
                    maxCols = MAXBLOCKSIZE;
                }
                maxRows = std::floor( MAXELEMSINBLOCK / maxCols );
                if( maxRows> nrows || maxRows + 1 == nrows) {
                    maxRows = nrows;
                }
                
            }    
            blockSize[0] = maxRows;
            blockSize[1] = maxCols;
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getMatrixBlockSize: "<<ex.what()<< " \n";
        }
        
        return(blockSize);
    }
    
    
    // Get the number of rows to read taking in to account the maximum elements per block
    // util when we have rectangular matrices, especially in omics data where we have
    // few samples and thousands of variables
    extern inline hsize_t getVectorBlockSize( int maxSize ) 
    {
        // size_t  maxSize = nrows * ncols;
        hsize_t blockSize = 0;
        
        try
        {
            if( maxSize > MAXELEMSINBLOCK) {
                blockSize = MAXELEMSINBLOCK;
            } else {
                blockSize = maxSize;
            }
            
        } catch(std::exception& ex) {
            Rcpp::Rcout<< "c++ exception getVectorBlockSize: "<<ex.what()<< " \n";
        }
        
        return(blockSize);
    }
    
    
    
    
    // extern inline Rcpp::IntegerVector getInitialPosition(bool transp, int desp )
    // {
    //     Rcpp::IntegerVector voffset(2);
    //     
    //     if(transp == true)
    //     {
    //         voffset[0] = 0;
    //         voffset[1] = desp;
    //     } else {
    //         voffset[0] = desp;
    //         voffset[1] = 0;
    //     }
    //     
    //     return(voffset);
    // }
    
    extern inline std::vector<hsize_t> getInitialPosition(bool transp, int desp )
    {
        std::vector<hsize_t> voffset = {0, 0};
        
        if(transp == true) {
            voffset[1] = desp;
        } else {
            voffset[0] = desp;
        }
        
        return(voffset);
    }
    
    
    // extern inline Rcpp::IntegerVector getSizetoRead(bool transp, int count, int rows, int cols )
    // {
    //     Rcpp::IntegerVector vcount(2);
    //     
    //     if(transp == true)
    //     {
    //         vcount[0] = rows;
    //         vcount[1] = count;
    //         
    //     } else {
    //         vcount[0] = count;
    //         vcount[1] = cols;
    //     }
    //     
    //     return(vcount);
    // }
    
    
    
    extern inline std::vector<hsize_t> getSizetoRead(bool transp, int count, int rows, int cols )
    {
        std::vector<hsize_t> vcount = {0, 0};

        if(transp == true)
        {
            vcount[0] = cols;
            vcount[1] = count;
            
        } else {
            vcount[0] = count;
            vcount[1] = cols;
        }
        
        return(vcount);
    }
    
    
    
    // Return the numbers of threads to be used in parallel processes
    extern inline unsigned int get_threads(bool bparal,  Rcpp::Nullable<int> threads  = R_NilValue) 
    {
        unsigned int ithreads;
        
        if(bparal == false) {
            ithreads = 1;
        } else {
            if(threads.isNotNull()) {
                if (Rcpp::as<int> (threads) <= std::thread::hardware_concurrency()){
                    ithreads = Rcpp::as<int> (threads);
                } else {
                    ithreads = getDTthreads(0, true);
                }
            } else {
                ithreads = getDTthreads(0, true);
            }    
        }
        
        return(ithreads);
    }
    
}

#endif // BIGDATASTATMETH_UTILITIES_HPP