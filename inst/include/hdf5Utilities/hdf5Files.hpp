#ifndef BIGDATASTATMETH_HDF5_FILES_HPP
#define BIGDATASTATMETH_HDF5_FILES_HPP

#include "BigDataStatMeth.hpp"
#include "Utilities/Utilities.hpp"
#include "hdf5Utilities/hdf5Utilities.hpp"

namespace BigDataStatMeth {

class hdf5File
{
    
public:
    
    hdf5File(std::string route, std::string filen, bool overwrite) 
    {
        filename = filen;
        boverwrite = overwrite;
        
        if( route.substr(route.length(), route.length()) == "/" ) {
            path = route.substr( 0, route.length()-1);
        } else {
            path = route;
        }
        
        if(path == "") {
            fullPath = filename;
        } else {
            fullPath = path + "/" + filename;    
        }
    }
    
    hdf5File(std::string route, std::string filen, H5::H5File* file, bool overwrite) 
    {
        filename = filen;
        boverwrite = overwrite;
        pfile = file;
        
        if( route.substr(route.length(), route.length()) == "/" ) {
            path = route.substr( 0, route.length()-1);
        } else {
            path = route;
        }
        
        if(path == "") {
            fullPath = filename;
        } else {
            fullPath = path + "/" + filename;    
        }
        
    }
    
    
    hdf5File(std::string filen, bool overwrite)
    {
        fullPath = filen;
        fullpath routefile = SplitElementName(filen);
        filename = routefile.filename;
        path = routefile.path;
        boverwrite = overwrite;
    }
    
    
    hdf5File(H5::H5File* file)
    {
        pfile = file;
    }

        
    // Create empty hdf5 data file
    int createFile() 
    {
        int iExec = BigDataStatMeth::EXEC_OK;
        
        try
        {
            
            H5::Exception::dontPrint();
            
            bool bFileExists = ResFileExist_filestream();
            
            if( !bFileExists || ( bFileExists && boverwrite) ) {
                pfile = new H5::H5File( fullPath, H5F_ACC_TRUNC ); 
            } else if ( bFileExists && !boverwrite){
                iExec = BigDataStatMeth::EXEC_WARNING;
            } else {
                Rcpp::Rcout<<"\n File exits, please set force = TRUE";
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception hdf5File (File IException) " );
        } 
        
        return(iExec);
    }
    
    
    // Open hdf5 data file
    //  opentype:
    //      r: read
    //      rw: read/write (default)
    H5::H5File* openFile(std::string opentype)
    {
        try
        {
            H5::Exception::dontPrint();
            
            //..// if( ResFileExist_filestream(fullPath) ) {
            if( ResFileExist_filestream() ) {
                if(opentype == "r") {
                    pfile = new H5::H5File( fullPath, H5F_ACC_RDONLY );
                } else {
                    pfile = new H5::H5File( fullPath, H5F_ACC_RDWR );
                }
            } else {
                Rcpp::Rcout<<"\n File does not exits, please create before open it";
                pfile = nullptr;
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception hdf5File (File IException) " );
        } 
        
        return(pfile);
    }

    H5::H5File* getFileptr() { return(pfile); }  // Return file pointer
    std::string getFilename() { return(filename); }  // Return file name
    std::string getPath() { return(path); }  // Return file path
    std::string getFullPath() { return(fullPath); }  // Return file path
    bool checkFile() { return(ResFileExist_filestream()); } // Return file exists
    
    
    // Destructor
    ~hdf5File(){
        pfile->close();
    }
    
protected:
    // Variables declaration
    H5::H5File* pfile = nullptr;
    std::string filename;
    std::string path;
    
private:
    
    // Variables declaration
    std::string opentype;
    std::string fullPath;
    bool boverwrite;
    
    
    // Function declarations
    
    // Test if file exsits
    bool ResFileExist_filestream() 
    {
        bool exists = true;
        std::fstream fileStream;
        fileStream.open(fullPath);
        
        if (fileStream.fail()) {
            exists = false;
        }
        return(exists);
    }
    
    void close_file() {
        pfile->close();
    }
    
};
}

#endif // BIGDATASTATMETH_HDF5_FILES_HPP