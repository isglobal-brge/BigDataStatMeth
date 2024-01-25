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
    
    
    
    // // Get dataset name inside a group
    // Rcpp::StringVector get_dataset_list( std::string strprefix)
    // {
    //     
    //     
    //     Rcpp::StringVector datasetnames;
    //     try{
    //         
    //         Rcpp::Rcout<<"Dins la funció Rcpp => del get_dataset_list\n";
    //         
    //         H5::Exception::dontPrint();
    //         
    //         H5::Group grp;
    //         herr_t err;
    //         ssize_t len;
    //         hsize_t nobj;
    //         int otype;
    //         char memb_name[MAX_NAME];
    //         
    //         
    //         if( exists_HDF5_element(pfile, groupname) ) {
    //             Rcpp::Rcout<<"Pas - 1\n";
    //             grp = pfile->openGroup(groupname);
    //         } else {
    //             ::Rf_error( "c++ exception group does not exists, please review group and file names" );
    //         }
    //         
    //         // get file id
    //         //..mogut....// Rcpp::Rcout<<"Pas - 1\n";
    //         //..mogut....// H5::Group grp = pfile->openGroup(groupname);
    //         Rcpp::Rcout<<"Pas - 2\n";
    //         hid_t gid = grp.getId();
    //         
    //         Rcpp::Rcout<<"Pas - 3\n";
    //         // get dataset names inside group
    //         err = H5Gget_num_objs(gid, &nobj);
    //         for (int i = 0; i < nobj; i++) 
    //         {
    //             Rcpp::Rcout<<"Pas - 4\n";
    //             len = H5Gget_objname_by_idx(gid, (hsize_t)i, memb_name, (size_t)MAX_NAME );
    //             
    //             Rcpp::Rcout<<"Pas - 5\n";
    //             otype =  H5Gget_objtype_by_idx(gid, (size_t)i );
    //             
    //             if( strprefix.compare("")!=0 ){
    //                 Rcpp::Rcout<<"Pas - 5.1\n";
    //                 if(otype == H5G_DATASET && (memb_name[0] == strprefix[0])) {
    //                     datasetnames.push_back(memb_name);
    //                 }
    //             } else {
    //                 Rcpp::Rcout<<"Pas - 5.2\n";
    //                 if(otype == H5G_DATASET ) {
    //                     datasetnames.push_back(memb_name);
    //                 }
    //             }
    //             
    //         }
    //     } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
    //         ::Rf_error( "c++ exception get_dataset_list (File IException)" );
    //         datasetnames = R_NilValue;
    //         return(datasetnames);
    //     } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
    //         ::Rf_error( "c++ exception get_dataset_list (Group IException)" );
    //         datasetnames = R_NilValue;
    //         return(datasetnames);
    //     } 
    //     
    //     return(datasetnames);
    //     
    // }
    
    
    
    
    
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