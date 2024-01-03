#ifndef BIGDATASTATMETH_HDF5_GROUPS_HPP
#define BIGDATASTATMETH_HDF5_GROUPS_HPP

#include "BigDataStatMeth.hpp"
#include "hdf5Utilities/hdf5Utilities.hpp"
#include "Utilities/Utilities.hpp"


namespace BigDataStatMeth {

class hdf5Group : public hdf5File
{
    
public:
    
    hdf5Group(std::string filename, std::string group) :
    hdf5File(filename, false)
    {
        openFile("rw");
        groupname = group;
    }
    
    
    hdf5Group(H5::H5File* file, std::string group) : 
    hdf5File(file)
    {
        openFile("rw");
        groupname = group;
    }
    
    hdf5Group(BigDataStatMeth::hdf5File* objFile, std::string group, bool forceGroup) : 
    hdf5File(objFile->getPath() , objFile->getFilename(), objFile->getFileptr(), false)
    {
        if( pfile ){
            openFile("rw");
        } else {
            ::Rf_error( "c++ exception Please create file before proceed" );
        }
        
        if( exists_HDF5_element(pfile, group) ) {
            if( forceGroup == true) {
                remove_elements(pfile, getGroupName(), {}); 
            } else {
                ::Rf_error( "c++ exception Please create file before proceed" );
            }
            
        }
        create_HDF5_groups(group);    
        
        groupname = group;
    }
    
    
    hdf5Group(BigDataStatMeth::hdf5File* objFile, std::string group) : 
    hdf5File(objFile->getPath() , objFile->getFilename(), objFile->getFileptr(), false)
    {
        if( pfile ){
            openFile("rw");
        } else {
            ::Rf_error( "c++ exception Please create file before proceed" );
        }
        
        if( !exists_HDF5_element(pfile, group) ) {
            create_HDF5_groups(group);    
        }
        
        groupname = group;
    }
    
    
    // Create multiple group in hdf5 data file, groups must be separated by "/"
    void create_HDF5_groups( H5std_string mGroup)
    {
        try
        {
            H5::Exception::dontPrint();
            
            char * pch;
            std::string strgroup = mGroup;
            char*  cpgroup = &strgroup[0];
            std::string results = "";
            
            pch = strtok(cpgroup, "/"); 
            
            while (pch != NULL)  
            {  
                if( results.compare("") == 0 ) {
                    results = pch;
                } else {
                    results = results + "/" + pch;
                }
                
                if(!pathExists( pfile->getId(), results )) {
                    pfile->createGroup(results);
                }
                pch = strtok (NULL, "/");  
            }  
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            // pfile->close();
            ::Rf_error( "c++ exception create_HDF5_groups_ptr (File IException)" );
            return void();
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            // pfile->close();
            ::Rf_error( "c++ exception create_HDF5_groups_ptr (Group IException)" );
            return void();
        } 
        
        groupname = mGroup;
        return void();
    }
    
    std::string getGroupName() { return(groupname); }  // Return group name
    
    
    // Destructor
    ~hdf5Group(){
    }
    
protected:
    // Variables declaration
    std::string groupname;
    
private:
    
    // Variables declaration
    
    // Function declarations
    
    
    
};
}

#endif // BIGDATASTATMETH_HDF5_GROUPS_HPP