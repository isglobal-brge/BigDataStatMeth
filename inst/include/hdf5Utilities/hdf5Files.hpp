/**
 * @file hdf5Files.hpp
 * @brief HDF5 file handling class and utilities
 * @details This header file provides a class for managing HDF5 files and related
 * operations. The implementation includes:
 * 
 * Key features:
 * - File creation and opening
 * - File status checking
 * - Dataset management
 * - Error handling
 * - Resource management
 */

#ifndef BIGDATASTATMETH_HDF5_FILES_HPP
#define BIGDATASTATMETH_HDF5_FILES_HPP


#include <RcppEigen.h>
#include "H5Cpp.h"

namespace BigDataStatMeth {

/**
 * @class hdf5File
 * @brief Class for managing HDF5 files
 * @details Provides functionality for creating, opening, and managing HDF5 files.
 * Includes methods for file operations, dataset management, and resource cleanup.
 */
class hdf5File
{
    
public:
    
    /**
     * @brief Constructor with path and filename
     * @param route Directory path
     * @param filen Filename
     * @param overwrite Whether to overwrite existing file
     */
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
    
    /**
     * @brief Constructor with path, filename, and file pointer
     * @param route Directory path
     * @param filen Filename
     * @param file HDF5 file pointer
     * @param overwrite Whether to overwrite existing file
     */
    hdf5File(std::string route, std::string filen, H5::H5File* file, bool overwrite) 
    {
        filename = filen;
        boverwrite = overwrite;
        pfile = file;
        bOwnsFile = false;
        
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
    
    
    /**
     * @brief Constructor with full path
     * @param filen Full path to file
     * @param overwrite Whether to overwrite existing file
     */
    hdf5File(std::string filen, bool overwrite)
    {
        fullPath = filen;
        fullpath routefile = SplitElementName(filen);
        filename = routefile.filename;
        path = routefile.path;
        boverwrite = overwrite;
    }
    
    
    /**
     * @brief Constructor with file pointer
     * @param file HDF5 file pointer
     */
    hdf5File(H5::H5File* file)
    {
        pfile = file;
        bOwnsFile = false;
    }


    /**
     * @brief Create a new HDF5 file
     * @details Creates a new HDF5 file at the specified location. If the file
     * exists and overwrite is true, it will be truncated.
     * 
     * @return EXEC_OK on success, EXEC_ERROR on error, EXEC_WARNING if file exists
     */
    int createFile() 
    {
        int iExec = EXEC_OK;
        
        try
        {
            H5::Exception::dontPrint();
            
            enable_hdf5_locking_once();
            
            bool bFileExists = ResFileExist_filestream();
            bool bInUse = bFileExists ? lockedByOtherProcess() : false;

            
            if( !bFileExists || ( bFileExists && boverwrite) ) {
                
                if (bInUse) {
                    Rf_error("HDF5 file is in use by another process; cannot overwrite.");
                }
                
                pfile = new H5::H5File( fullPath, H5F_ACC_TRUNC ); 
                bOwnsFile = true;
                iExec = EXEC_OK; //.. 2025/08/13 ..//
            } else if ( bFileExists && !boverwrite){
                iExec = EXEC_WARNING;
            } else {
                Rcpp::Rcout<<"\n File exits, please set force = TRUE";
            }    
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            Rf_error("c++ exception hdf5File (File IException) " );
        } 
        
        return(iExec);
    }
    
    
    /**
     * @brief Open an existing HDF5 file
     * @details Opens an HDF5 file in read or read/write mode.
     * 
     * @param opentype Access mode ("r" for read-only, "rw" for read-write)
     * @return Pointer to opened file or nullptr on error
     */
    H5::H5File* openFile(std::string opentype)
    {
        try
        {
            H5::Exception::dontPrint();
            
            enable_hdf5_locking_once();
            
            if( checkHDF5File() ) {
                if(opentype == "r") {
                    pfile = new H5::H5File( fullPath, H5F_ACC_RDONLY );
                    bOwnsFile = true;
                } else {
                    if (lockedByOtherProcess()) {
                        Rf_error("HDF5 file is in use by another process.");
                    }
                    
                    pfile = new H5::H5File( fullPath, H5F_ACC_RDWR );
                    bOwnsFile = true;
                }
            } else {
                
                if (opentype == "r") {        //..2025/08/13..//
                    Rf_error("HDF5 file not found."); //..2025/08/13..//
                }         //..2025/08/13..//
                pfile = new H5::H5File(fullPath, H5F_ACC_TRUNC); //..2025/08/13..//
                bOwnsFile = true;
            }
            
        } catch (const H5::Exception& e) {
            pfile = nullptr;
            std::string error_msg = std::string("openFile HDF5 error: ") + e.getDetailMsg();
            Rf_error("%s", error_msg.c_str());
        } catch (const std::exception& e) {
            pfile = nullptr;
            std::string error_msg = std::string("openFile error: ") + e.what();
            Rf_error("%s", error_msg.c_str());
        }
        
        return(pfile);
    }

    /**
     * @brief Get file pointer
     * @return Pointer to HDF5 file
     */
    H5::H5File* getFileptr() { return(pfile); }  // Return file pointer

    /**
     * @brief Get filename
     * @return Filename without path
     */
    std::string getFilename() { return(filename); }  // Return file name

    /**
     * @brief Get file path
     * @return Directory path without filename
     */
    std::string getPath() { return(path); }  // Return file path

    /**
     * @brief Get full file path
     * @return Complete path including filename
     */
    std::string getFullPath() { return(fullPath); }  // Return file path

    /**
     * @brief Check if file exists
     * @return True if file exists, false otherwise
     */
    bool checkFile() { return(ResFileExist_filestream()); } // Return file exists

    /**
     * @brief Get list of dataset names
     * @param strgroup Group path
     * @param strprefix Prefix filter
     * @param strsufix Suffix filter
     * @return Vector of dataset names
     */
    Rcpp::StringVector getDatasetNames( std::string strgroup, std::string strprefix, std::string strsufix){ 
        return(get_dataset_names_from_group( strgroup, strprefix, strsufix)); 
    } // Return a dataset name list with all the datasets inside
    
    
    /**
     * @brief List datasets in a group, optionally recursing into subgroups.
     *
     * @param strgroup  HDF5 group path. Use "/" for the file root.
     * @param strprefix Only return names starting with this prefix (empty = all).
     * @param recursive If true, recurse into every subgroup and return full paths
     *                  relative to \p strgroup (e.g. "SUB/dataset").
     * @return StringVector of dataset names (or relative paths when recursive).
     */
    Rcpp::StringVector getAllDatasetNames(std::string strgroup,
                                          std::string strprefix,
                                          bool        recursive) {
        if (!recursive) {
            return get_dataset_names_from_group(strgroup, strprefix, "");
        }
        Rcpp::StringVector results;
        collect_datasets_recursive(strgroup, "", strprefix, results);
        return results;
    }
    
    
    
    /**
     * @brief Close file and cleanup resources
     * @details Closes all open objects and the file itself. Used for emergency cleanup.
     */
    void close_file() {
        try {
            
            // Obtener el nÃºmero de objetos abiertos
            ssize_t sObjects = H5Fget_obj_count( pfile->getId(), H5F_OBJ_ALL);
            
            if (sObjects > 0) {
                // Get the Ids of opened objects 
                std::vector<hid_t> vIds(sObjects);
                ssize_t sOpenIds = H5Fget_obj_ids(pfile->getId(), H5F_OBJ_ALL, sObjects, vIds.data());
                
                if (sOpenIds > 0) {
                    for (hid_t OpenId : vIds) {
                        if (H5Iis_valid(OpenId)) {
                            H5Oclose(OpenId); // Cierra el objeto genÃ©rico (puede ser dataset, grupo, etc.)
                        }
                    }
                } 
            }
            
            pfile->close();
        } catch(std::exception& ex) {
            
            Rcpp::Rcerr<< "c++ exception close_file (err FileException)";
        }
        
        return void();
    }
    
    
    /**
     * @brief Test whether an HDF5 file is locked / in use.
     *
     * @param filename Path to the HDF5 file (relative or absolute).
     *
     * @return true if the file appears locked (i.e., cannot be opened in
     *         read/write mode under HDF5 locking); false otherwise. If the
     *         file does not exist, returns false.
     *
     * @details Constructs a lightweight temporary @c hdf5File and delegates
     *          to @c lockedByOtherProcess(). It does not create, truncate,
     *          or keep the file open; no state is persisted.
     *
     * @note On systems without HDF5 file locking (or if disabled), the
     *       result may be conservative. Lack of write permissions may also
     *       yield @c true (indistinguishable from "locked").
     *
     * @since 0.99.0
     */
    static bool isLocked(const std::string& filename) noexcept {
        hdf5File tmp(filename, /*overwrite=*/false);
        return tmp.lockedByOtherProcess();
    }
    
    
    
    /**
     * @brief Destructor
     * @details Closes the file and releases resources
     */
    ~hdf5File(){
        try {
            if(pfile != nullptr && bOwnsFile) {
                pfile->close();
                pfile = nullptr;
            }
        } catch (...) {}
    }
    
protected:
    // Variables declaration
    H5::H5File* pfile = nullptr;
    std::string filename;
    std::string path;
    bool bOwnsFile = false;
    
private:
    
    // Variables declaration
    std::string opentype;
    std::string fullPath;
    bool boverwrite;
    
    // hdf5 lock file
    
    #ifdef _WIN32
    #include <cstdlib>
        static inline void enable_hdf5_locking_once() {
            static bool done = (_putenv_s("HDF5_USE_FILE_LOCKING","TRUE"), true);
            (void)done;
        }
    #else
    #include <cstdlib>
        static inline void enable_hdf5_locking_once() {
            static bool done = (setenv("HDF5_USE_FILE_LOCKING","TRUE",1), true);
            (void)done;
        }
    #endif
    
    // Function
    
    #if __cplusplus >= 201703L // C++17 and later 
    #include <string_view>
        
        /**
         * @brief Check if string ends with suffix
         * @param str String to check
         * @param suffix Suffix to match
         * @return True if string ends with suffix
         */
        static bool ends_with(std::string_view str, std::string_view suffix)
        {
            return str.size() >= suffix.size() && str.compare(str.size()-suffix.size(), suffix.size(), suffix) == 0;
        }
        
        /**
         * @brief Check if string starts with prefix
         * @param str String to check
         * @param prefix Prefix to match
         * @return True if string starts with prefix
         */
        static bool starts_with(std::string_view str, std::string_view prefix)
        {
            return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
        }
        
    #endif // C++17
    
    
    /**
     * @brief Check if file exists using file stream
     * @return True if file exists and is accessible
     */
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

    
    /**
     * @brief Check if file is corrupt, open, accessible or has_valid_structure
     * @return True if file exists and is accessible
     */
    bool checkHDF5File() {
        
        bool is_accessible = false;
        std::string error_message = "";
        
        try {
            // Turn off automatic error printing
            H5::Exception::dontPrint();
            
            // Method 1: Check if file is accessible
            try {
                // Try to check if file exists and is HDF5 format
                if (H5::H5File::isHdf5(fullPath)) {
                    is_accessible = true;
                } else {
                    Rf_error("c++ exception File is not in HDF5 format" );
                }
            } catch (const H5::FileIException& e) {
                error_message = "c++ exception File access error: " + std::string(e.getCDetailMsg());
                Rf_error("%s", error_message.c_str() );
            } catch (const H5::Exception& e) {
                error_message = "c++ exception HDF5 Exception during accessibility check: " + std::string(e.getCDetailMsg());
                Rf_error("%s", error_message.c_str() );
            }
            
            // Method 2: Try to open the file if accessible
            if (is_accessible) {
                try {
                    H5::H5File* file = new H5::H5File(fullPath, H5F_ACC_RDONLY);
                    // is_open = true;
                    
                    // Method 3: Validate file structure
                    try {
                        // Try to access root group
                        H5::Group root_group = file->openGroup("/");
                        
                        root_group.close();
                        
                    } catch (const H5::GroupIException& e) {
                        error_message =  "c++ exception (checkHDF5File) Root group access failed: " + std::string(e.getCDetailMsg());
                        Rf_error("%s", error_message.c_str() );
                    } catch (const H5::Exception& e) {
                        error_message =  "c++ exception (checkHDF5File) Structure validation failed: " + std::string(e.getCDetailMsg() );
                        Rf_error("%s", error_message.c_str() );
                    }
                    
                    // Close the file
                    file->close();
                    delete file;
                    
                } catch (const H5::FileIException& e) {
                    error_message = "c++ exception (checkHDF5File) Cannot open file: " + std::string(e.getCDetailMsg());
                    Rf_error("%s", error_message.c_str() );
                } catch (const H5::Exception& e) {
                    error_message ="c++ exception HDF5 Exception during file opening: " + std::string(e.getCDetailMsg() );
                    Rf_error("%s", error_message.c_str() );
                }
            }
            
        } catch (const std::exception& e) {
            error_message = "c++ exception (checkHDF5File): " + std::string(e.what());
            Rf_error("%s", error_message.c_str() );
            // ::Rf_error( error_message.c_str() );
            
            // error_message = "Standard exception: " + std::string(e.what());
            // is_corrupt = true;
        } catch (...) {
            error_message = "c++ exception (checkHDF5File): Unknown exception occurred" ;
            Rf_error("%s", error_message.c_str() );
        }
        
        return(true);
        
    }
    
    
    /**
     * @brief Return true if existing HDF5 file appears locked/busy.
     * @note Requires HDF5 file locking enabled (env var set above).
     */
    bool lockedByOtherProcess() {
        
        if (!ResFileExist_filestream()) {
            return false;    
        }
        
        H5::Exception::dontPrint();
        enable_hdf5_locking_once();
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        
        #if H5_VERSION_GE(1,12,0)
                H5Pset_file_locking(fapl, 1 , 0 );
        #endif
                
        hid_t fid = H5Fopen(fullPath.c_str(), H5F_ACC_RDWR, fapl);
        H5Pclose(fapl);
        
        if (fid < 0) {
            return true;       
        }
        
        H5Fclose(fid);
        return false;
    }
    
    
    
    
    
    

    /**
     * @brief Check if HDF5 file is already open
     * @return True if file is open
     */
    bool isHDF5FileOpen() 
    {
        
        bool bCorrupt = false;
        
        hid_t file_id = H5Fopen(fullPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        
        if (file_id < 0) {
            bCorrupt = true; 
        } else {
            // opend objects
            int num_open_objects = H5Fget_obj_count(file_id, H5F_OBJ_ALL);
            
            if(num_open_objects > 1 ) 
                bCorrupt = true; 
            
            H5Fclose(file_id);
        }
        
        return bCorrupt;
        
    }
    
    
    /**
     * @brief Recursive helper: walk \p group_path and collect dataset paths.
     *
     * @param group_path   Current group being visited.
     * @param rel_prefix   Path prefix accumulated so far (empty at root call).
     * @param strprefix    User-supplied name filter (empty = accept all).
     * @param results      Accumulator vector — paths appended in-place.
     */
    void collect_datasets_recursive(const std::string&  group_path,
                                    const std::string&  rel_prefix,
                                    const std::string&  strprefix,
                                    Rcpp::StringVector& results)
    {
        try {
            H5::Exception::dontPrint();
            
            hsize_t nobj = 0;
            char    memb_name[MAX_NAME];
            
            H5::Group grp = pfile->openGroup(group_path);
            hid_t     gid = grp.getId();
            
            H5Gget_num_objs(gid, &nobj);
            
            for (hsize_t i = 0; i < nobj; ++i) {
                ssize_t len = H5Gget_objname_by_idx(gid, i, memb_name,
                                                    (size_t)MAX_NAME);
                if (len <= 0) continue;
                
                int otype = H5Gget_objtype_by_idx(gid, (size_t)i);
                std::string name(memb_name);
                std::string rel_path = rel_prefix.empty()
                    ? name
                : rel_prefix + "/" + name;
                
                if (otype == H5G_DATASET) {
                    if (strprefix.empty() || starts_with(name, strprefix)) {
                        results.push_back(rel_path);
                    }
                } else if (otype == H5G_GROUP) {
                    std::string child_path = (group_path == "/")
                    ? "/" + name
                    : group_path + "/" + name;
                    collect_datasets_recursive(child_path, rel_path,
                                               strprefix, results);
                }
            }
        } catch (...) { /* skip unreadable groups */ }
    }
    
    
    /**
     * @brief Get dataset names from group
     * @param strgroup Group path
     * @param strprefix Prefix filter
     * @param strsufix Suffix filter
     * @return Vector of dataset names
     */
    Rcpp::StringVector get_dataset_names_from_group(std::string strgroup, std::string strprefix, std::string strsufix)
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
                Rcpp::Rcerr<<"\nc++ exception get_dataset_names_from_group (err IException)\n";
                return -1;
            } else {
                for (unsigned int i = 0; i < nobj; i++) 
                {
                    len = H5Gget_objname_by_idx(gid, (hsize_t)i, memb_name, (size_t)MAX_NAME );
                    
                    if(len == 0) {
                        Rcpp::Rcerr<<"c++ exception get_dataset_names_from_group (len IException)\n";
                        return -1;
                    }
                    
                    otype =  H5Gget_objtype_by_idx(gid, (size_t)i );
                    
                    // 202505
                    if( strprefix.compare("")!=0 && strsufix.compare("")==0 ){
                        if(otype == H5G_DATASET && (starts_with(memb_name, strprefix)) ) {
                            datasetnames.push_back(memb_name);
                        }
                    } else if( strprefix.compare("")==0 && strsufix.compare("")!=0 ){
                        if(otype == H5G_DATASET && (ends_with(memb_name, strsufix)) ) {
                            datasetnames.push_back(memb_name);
                        }
                    } else if( strprefix.compare("")!=0 && strsufix.compare("")!=0 ) {
                        if(otype == H5G_DATASET && starts_with(memb_name, strprefix) && ends_with(memb_name, strsufix)  ) {
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
            Rcpp::Rcerr<<"\nc++ exception get_dataset_names_from_group (File IException)\n";
            return -1;
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception get_dataset_names_from_group (DataSet IException)\n";
            return -1;
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            Rcpp::Rcerr<<"\nc++ exception get_dataset_names_from_group (Group IException)\n";
            return -1;
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception get_dataset_names_from_group (DataSpace IException)\n";
            return -1;
        } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception get_dataset_names_from_group (Data TypeIException)\n";
            return -1;
        }
        return(datasetnames);
    }
    
};
}

#endif // BIGDATASTATMETH_HDF5_FILES_HPP