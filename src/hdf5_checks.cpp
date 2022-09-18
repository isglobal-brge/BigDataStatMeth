#include "include/hdf5_checks.h"

bool exist_FileGroupDataset(std::string filename, std::string group, std::string dataset)
{
  
  H5File* file = nullptr;
  // DataSet* pdataset = nullptr;
  
  try {
    
      if( ResFileExist_filestream(filename) ) {
        file = new H5File( filename, H5F_ACC_RDONLY ); 
      } else {
        Rcpp::Rcout<<"\nFile not exits, create file before split dataset";
        return false;
      }
    
      if( group.compare("") != 0 ) {
        
        if( exists_HDF5_element_ptr(file, group) != 0 ) {
            
          if( dataset.compare("") != 0 ) {
            if( exists_HDF5_element_ptr(file, group + "/" + dataset ) == 0 ) {
                  file->close();
                  Rcpp::Rcout<<"Group not exists, create the input dataset before proceed";
                  return false;
            }
          }
          
        }  else { 
      
          file->close();
          Rcpp::Rcout<<"Group not exists, create the group and dataset before proceed";
          return false;
          
        }
      }
    
  } catch( FileIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception exist_FileGroupDataset (File IException)" );
    return false;
  } catch( GroupIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception exist_FileGroupDataset (Group IException)" );
    return false;
  } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
    file->close();
    ::Rf_error( "c++ exception exist_FileGroupDataset (DataSet IException)" );
    return false;
  } catch(std::exception& ex) {
    Rcpp::Rcout<< ex.what();
    return false;
  }

  file->close();
  return true;
  
}



double prepare_outGroup(H5File* file, std::string outGroup, bool bforce)
{
  
  double res = 0;
  
  try {
    
    if( exists_HDF5_element_ptr(file, outGroup) == 0 ) {
      res = create_HDF5_groups_ptr(file, outGroup );
    } /*  Commented to prevent to remove complete group
          maybe there are other interesting datasets inside the group !!!
 
    else if (exists_HDF5_element_ptr(file, outGroup) != 0 && bforce == true)
    {
      res = remove_HDF5_element_ptr(file, outGroup);
      res = create_HDF5_group_ptr(file, outGroup );
    } else {
      throw std::range_error("Group exists, please set force = true if you want to rewrite data");
    }*/
    
    
  } catch( FileIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception prepare_outGroup (File IException)" );
    return (res);
  } catch( GroupIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception prepare_outGroup (Group IException)" );
    return (res);
  } catch(std::exception& ex) {
    Rcpp::Rcout<< ex.what();
    return (res);
  }
  
  return(res);
  
  
}


double prepare_outDataset(H5File* file, std::string outDataset, bool bforce)
{
    
    double res = 0;
    
    try {
        
        if( exists_HDF5_element_ptr(file, outDataset) != 0 && bforce == false) {
            Rcpp::Rcout<<"Output dataset exists, please set force = true if you want to rewrite data";
            //..// throw std::range_error("Output dataset exists, please set force = true if you want to rewrite data");
        } else if ( exists_HDF5_element_ptr(file, outDataset) !=0 && bforce == true) {
          
            res = remove_HDF5_element_ptr(file, outDataset);
        }
        

    } catch( FileIException& error ) { // catch failure caused by the H5File operations
        file->close();
        ::Rf_error( "c++ exception prepare_outGroup (File IException)" );
        return (res);
    } catch( GroupIException& error ) { // catch failure caused by the H5File operations
        file->close();
        ::Rf_error( "c++ exception prepare_outGroup (Group IException)" );
        return (res);
    } catch(std::exception& ex) {
        Rcpp::Rcout<< ex.what();
        return (res);
    }
    
    return(res);
    
    
}

