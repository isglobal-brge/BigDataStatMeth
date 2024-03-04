#ifndef BIGDATASTATMETH_HDF5_FILES_HPP
#define BIGDATASTATMETH_HDF5_FILES_HPP

#include "BigDataStatMeth.hpp"
#include "Utilities/Utilities.hpp"
// #include <filesystem>


// // #include "hdf5Utilities/hdf5Utilities.hpp"

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
        int iExec = EXEC_OK;
        
        try
        {
            
            H5::Exception::dontPrint();
            
            bool bFileExists = ResFileExist_filestream();
            
            if( !bFileExists || ( bFileExists && boverwrite) ) {
                pfile = new H5::H5File( fullPath, H5F_ACC_TRUNC ); 
            } else if ( bFileExists && !boverwrite){
                iExec = EXEC_WARNING;
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
            
            bool bFileExists = ResFileExist_filestream();
            
            if( bFileExists ) {
                if(opentype == "r") {
                    pfile = new H5::H5File( fullPath, H5F_ACC_RDONLY );
                } else {
                    pfile = new H5::H5File( fullPath, H5F_ACC_RDWR );
                }
            } else {
                ::Rf_error("\n File does not exits, please create before open it");
                pfile = nullptr;
                return(pfile);
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            pfile = new H5::H5File(fullPath, H5F_ACC_TRUNC);
            // ::Rf_error( "c++ exception hdf5File (File IException) " );
        } 
        
        return(pfile);
    }

    H5::H5File* getFileptr() { return(pfile); }  // Return file pointer
    std::string getFilename() { return(filename); }  // Return file name
    std::string getPath() { return(path); }  // Return file path
    std::string getFullPath() { return(fullPath); }  // Return file path
    bool checkFile() { return(ResFileExist_filestream()); } // Return file exists
    Rcpp::StringVector getDatasetNames( std::string strgroup, std::string strprefix){ 
        return(get_dataset_names_from_group( strgroup, strprefix)); 
    } // Return a dataset name list with all the datasets inside
    
    
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
        
        bool exists = false;
        
        std::fstream fileStream;
        fileStream.open(fullPath);
        
        if (fileStream.good()) {
            exists = true;
        } else {
            exists = false;
        }
        return(exists);
    }

    
    void close_file() {
        pfile->close();
    }
    
    
    // Get dataset names inside the strgroup
    Rcpp::StringVector get_dataset_names_from_group(std::string strgroup, std::string strprefix)
    {
        
        Rcpp::StringVector datasetnames;
        
        try{
            
            H5::Exception::dontPrint();
            
            herr_t err;
            ssize_t len;
            hsize_t nobj;
            int otype;
            char memb_name[MAX_NAME];
            
            // get file id
            H5::Group grp = pfile->openGroup(strgroup);
            hid_t gid = grp.getId();
            
            // get dataset names inside group
            err = H5Gget_num_objs(gid, &nobj);
            if(err<0 ) {
                ::Rf_error( "c++ exception get_dataset_names_from_group (err IException)" );
                return -1;
            } else {
                for (int i = 0; i < nobj; i++) 
                {
                    len = H5Gget_objname_by_idx(gid, (hsize_t)i, memb_name, (size_t)MAX_NAME );
                    
                    if(len == 0) {
                        ::Rf_error( "c++ exception get_dataset_names_from_group (len IException)" );
                        return -1;
                    }
                    
                    otype =  H5Gget_objtype_by_idx(gid, (size_t)i );
                    
                    // 202109
                    if( strprefix.compare("")!=0 ){
                        if(otype == H5G_DATASET && (memb_name[0] == strprefix[0])) {
                            datasetnames.push_back(memb_name);
                        }
                    } else {
                        if(otype == H5G_DATASET ) {
                            datasetnames.push_back(memb_name);
                        }
                    }
                }    
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception get_dataset_names_from_group (File IException)" );
            return -1;
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            ::Rf_error( "c++ exception get_dataset_names_from_group (DataSet IException)" );
            return -1;
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            ::Rf_error( "c++ exception get_dataset_names_from_group (Group IException)" );
            return -1;
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception get_dataset_names_from_group (DataSpace IException)" );
            return -1;
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception get_dataset_names_from_group (Data TypeIException)" );
            return -1;
        }
        return(datasetnames);
    }
    
};
}

#endif // BIGDATASTATMETH_HDF5_FILES_HPP