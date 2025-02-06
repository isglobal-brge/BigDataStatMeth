#ifndef BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP
#define BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP


namespace BigDataStatMeth {

    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } 
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"\nc++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"\nC++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                         BigDataStatMeth::hdf5Dataset* ds2) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } 
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"\nc++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"\nC++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                         BigDataStatMeth::hdf5Dataset* ds2, 
                         BigDataStatMeth::hdf5Dataset* ds3) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"\nc++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"\nC++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                         BigDataStatMeth::hdf5Dataset* ds2, 
                         BigDataStatMeth::hdf5Dataset* ds3, 
                         BigDataStatMeth::hdf5Dataset* ds4) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } 
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"\nc++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"\nC++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    
    
    
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                                BigDataStatMeth::hdf5Dataset* ds2, 
                                BigDataStatMeth::hdf5Dataset* ds3, 
                                BigDataStatMeth::hdf5Dataset* ds4, 
                                BigDataStatMeth::hdf5Dataset* ds5) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } else if(ds5->getDatasetptr() !=nullptr){
                ds5->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"\nc++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"\nC++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }



    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                                       BigDataStatMeth::hdf5Dataset* ds2, 
                                       BigDataStatMeth::hdf5Dataset* ds3, 
                                       BigDataStatMeth::hdf5Dataset* ds4, 
                                       BigDataStatMeth::hdf5Dataset* ds5,
                                       BigDataStatMeth::hdf5Dataset* ds6) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } else if(ds5->getDatasetptr() !=nullptr){
                ds5->close_file();
            } else if(ds6->getDatasetptr() !=nullptr){
                ds6->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"\nc++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"\nC++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }
    
    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                                       BigDataStatMeth::hdf5Dataset* ds2, 
                                       BigDataStatMeth::hdf5Dataset* ds3, 
                                       BigDataStatMeth::hdf5Dataset* ds4, 
                                       BigDataStatMeth::hdf5Dataset* ds5,
                                       BigDataStatMeth::hdf5Dataset* ds6,
                                       BigDataStatMeth::hdf5Dataset* ds7) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } else if(ds5->getDatasetptr() !=nullptr){
                ds5->close_file();
            } else if(ds6->getDatasetptr() !=nullptr){
                ds6->close_file();
            } else if(ds7->getDatasetptr() !=nullptr){
                ds7->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"\nc++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"\nC++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }


    extern inline void checkClose_file(BigDataStatMeth::hdf5Dataset* ds1, 
                                       BigDataStatMeth::hdf5Dataset* ds2, 
                                       BigDataStatMeth::hdf5Dataset* ds3, 
                                       BigDataStatMeth::hdf5Dataset* ds4, 
                                       BigDataStatMeth::hdf5Dataset* ds5,
                                       BigDataStatMeth::hdf5Dataset* ds6,
                                       BigDataStatMeth::hdf5Dataset* ds7,
                                       BigDataStatMeth::hdf5Dataset* ds8) 
    {
        try{
            
            if(ds1->getDatasetptr() !=nullptr){
                ds1->close_file();
            } else if(ds2->getDatasetptr() !=nullptr){
                ds2->close_file();
            } else if(ds3->getDatasetptr() !=nullptr){
                ds3->close_file();
            } else if(ds4->getDatasetptr() !=nullptr){
                ds4->close_file();
            } else if(ds5->getDatasetptr() !=nullptr){
                ds5->close_file();
            } else if(ds6->getDatasetptr() !=nullptr){
                ds6->close_file();
            } else if(ds7->getDatasetptr() !=nullptr){
                ds7->close_file();
            } else if(ds8->getDatasetptr() !=nullptr){
                ds8->close_file();
            }
            
        } catch( H5::FileIException& error ) { // catch failure caused by the H5File operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (File IException)\n";
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the DataSet operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSet IException)\n";
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            Rcpp::Rcerr<<"\nc++ exception checkClose_file (DataSpace IException)\n";
            return void();
        } catch(std::exception &ex) {
            Rcpp::Rcout<<"\nc++ exception checkClose_file \n"<< ex.what();
            return void();
        } catch (...) {
            Rcpp::Rcout<<"\nC++ exception checkClose_file (unknown reason)";
            return void();
        }
        
        return void();
    }

}

#endif // BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP