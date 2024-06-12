#ifndef BIGDATASTATMETH_UTIL_REMOVE_ELEMENT_HPP
#define BIGDATASTATMETH_UTIL_REMOVE_ELEMENT_HPP



namespace BigDataStatMeth {


    void RcppRemove_hdf5_elements(BigDataStatMeth::hdf5File* file, std::vector<std::string> elements)
    {
        
        try
        {
            H5::Exception::dontPrint();
            
            // Remove group
            if(elements.size() == 0) {
                std::string strmessage = "Nothing to be removed removed";
                Rcpp::message(Rcpp::wrap(strmessage));
            } else { // Remove datasets
                for (int i=0; i<elements.size(); i++) 
                {
                    H5std_string element = "" + elements[i];
                    if(exists_HDF5_element( file->getFileptr(), element))
                    {
                        int result = H5Ldelete( (file->getFileptr())->getId(), element.data(), H5P_DEFAULT);  
                        if(result<0) {
                            Rcpp::Rcout<<"\n Error removing : "<<elements[i]<<"\n";
                        } 
                    } else {
                        Rcpp::Rcout<<"\n Element: "<<elements[i]<<" does not exists \n";
                    }
                }    
            }
            
        } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (File IException)" );
            return void();
        } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (Group IException)" );
            return void();
        } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (DataSet IException)" );
            return void();
        } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception remove_HDF5_multiple_elements_ptr (DataSpace IException)" );
            return void();
        } 
        return void();
    }

}

#endif // BIGDATASTATMETH_UTIL_REMOVE_ELEMENT_HPP