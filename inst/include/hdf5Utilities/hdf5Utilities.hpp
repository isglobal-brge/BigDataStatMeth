#ifndef BIGDATASTATMETH_HDF5_UTILITIES_HPP
#define BIGDATASTATMETH_HDF5_UTILITIES_HPP

// #include "BigDataStatMeth.hpp"


namespace BigDataStatMeth {

// // Constants
// 
//     #define MAX_NAME 1024
//     //. ORIGINAL - CORRECTE.//#define MAXSVDBLOCK 1500
//     #define MAXSVDBLOCK 10 //..NOMÉS DEBUG..//
//     const int RANK1 = 1;
//     const int RANK2 = 2;
//     const int RANK3 = 3;
//     const int DIM1 = 1;
//     const int DIM2 = 2;
//     const int DIM3 = 3;
//     const int MAXSTRING = 20;
//     const hsize_t MAXSTRBLOCK = 100000;
//     const hsize_t MAXELEMSINBLOCK = 250000;
//     const int maxElemBlock = 3000000;
//     //.Only test.// const int maxElemBlock = 30;
//     const int EXEC_OK = 0;
//     const int EXEC_ERROR = 1;
//     const int EXEC_WARNING = 2;



// Functions

    extern inline bool pathExists(hid_t id, const std::string& path)
    {
        try {
            return H5Lexists( id, path.c_str(), H5P_DEFAULT ) > 0;    
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception pathExists (File IException)" );
            return false;
        } 
        
    }
    
    // Check if HDF5 element (group or dataset) exists
    extern inline bool exists_HDF5_element(H5::H5File* file, std::string element)
    {
        bool bexists = false;
        try
        {
            // H5::Exception::dontPrint();
         
             if( element.substr(element.length(), element.length()) == "/" ) {
                 element = element.substr( 0, element.length()-1);
             } 
             
            // Search dataset
            if(pathExists( file->getId(), element)) 
                bexists = true;
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            file->close();
            ::Rf_error( "c++ exception exists_HDF5_element (File IException)" );
            return bexists;
        } 
        return bexists;
    }
    
    
    extern inline bool remove_elements(H5::H5File* file, std::string strgroup, Rcpp::StringVector elements)
    {
        
        bool bremok = true;
        
        try
        {
            H5::Exception::dontPrint();
            
            // Remove group
            if(elements.size() == 0) {
                H5std_string element = strgroup;
                
                int result = H5Ldelete(file->getId(), element.c_str(), H5P_DEFAULT);  
                if(result<0) {
                    Rcpp::Rcout<<"\n Error removing group: "<<element<<"\n";
                    bremok = false;
                } 
                
            } else { // Remove datasets
                for (int i=0; i<elements.size(); i++) 
                {
                    H5std_string element = strgroup + "/" + elements[i];
                    
                    int result = H5Ldelete(file->getId(), element.data(), H5P_DEFAULT);  
                    if(result<0) {
                        Rcpp::Rcout<<"\n Error removing : "<<element<<"\n";
                        bremok = false;
                    } 
                }    
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (File IException)" );
            return(bremok);
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (Group IException)" );
            return(bremok);
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (DataSet IException)" );
            return(bremok);
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (DataSpace IException)" );
            return(bremok);
        } 
        
        return(bremok);
    }


    // Remove single element, could be a group or dataset, we need to set the full path in element parameter
    //..// extern inline bool remove_elements(H5::H5File* file, std::string strgroup, Rcpp::StringVector elements)
    extern inline bool remove_elements(H5::H5File* file, H5std_string element)
    {
        
        bool bremok = true;
        
        try
        {
            H5::Exception::dontPrint();
            
            // H5std_string elementtoremove = element;
            
            int result = H5Ldelete(file->getId(), element.data(), H5P_DEFAULT);  
            if(result<0) {
                Rcpp::Rcout<<"\n Error removing : "<<element<<"\n";
                bremok = false;
            } 
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (File IException)" );
            return(bremok);
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (Group IException)" );
            return(bremok);
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (DataSet IException)" );
            return(bremok);
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (DataSpace IException)" );
            return(bremok);
        } 
        
        return(bremok);
    }
    
    
    /////////////////////
    // NOT TESTED !!!! //
    /////////////////////
    extern inline void createHardLink( H5::H5File* file, std::string original, std::string link)
    {
        
        try{
            
            H5::Exception::dontPrint();
            
            const char * charOriginal = original.c_str();
            const char * charLink = link.c_str();
            
            herr_t status = H5Lcreate_hard(file->getId(), charOriginal, file->getId(), charLink, H5P_DEFAULT, H5P_DEFAULT);
            
            if(status<0) {
                ::Rf_error( "c++ exception createHardLink (create_hard IException)" );
            }
            
        } catch(H5::FileIException& error) { 
            ::Rf_error( "c++ exception createHardLink (File IException)" );
        } catch(H5::DataSetIException& error) { 
            ::Rf_error( "c++ exception createHardLink (DataSet IException)" );
        } catch(H5::GroupIException& error) { 
            ::Rf_error( "c++ exception createHardLink (Group IException)" );
        } 
        
        return void();
    }


    extern inline void renameElement( H5::H5File* file, std::string original, std::string link)
    {
        
        try{
            
            H5::Exception::dontPrint();
            
            const char * charOriginal = original.c_str();
            const char * charLink = link.c_str();
            
            // Rcpp::Rcout<<"original: "<< original<<" - Desti: "<<link<<"\n";
            
            herr_t status = H5Lmove(file->getId(), charOriginal, file->getId(), charLink, H5P_DEFAULT, H5P_DEFAULT);
            
            if(status<0) {
                ::Rf_error( "c++ exception renameElement (rename_element IException)" );
            } 
            
            
        } catch(H5::FileIException& error) { 
            ::Rf_error( "c++ exception renameElement (File IException)" );
        } catch(H5::DataSetIException& error) { 
            ::Rf_error( "c++ exception renameElement (DataSet IException)" );
        } catch(H5::GroupIException& error) { 
            ::Rf_error( "c++ exception renameElement (Group IException)" );
        } 
        
        return void();
    }


    extern inline void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
    {
        unsigned int numRows = matrix.rows()-1;
        unsigned int numCols = matrix.cols();
        
        if( rowToRemove < numRows )
            matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove).eval();
        
        matrix.conservativeResize(numRows,numCols);
    }
    
    extern inline void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
    {
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols()-1;
        
        if( colToRemove < numCols )
            matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove).eval();
        
        matrix.conservativeResize(numRows,numCols);
    }
    

}

#endif // BIGDATASTATMETH_HDF5_UTILITIES_HPP

