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
    
    
    // 
    // // Join multiple datasets in one dataset in the same group
    // 
    // template< typename T>
    // extern inline int join_datasets ( T* dsJoined, std::string strsubgroup, Rcpp::StringVector strinput, bool bremoveJoined )
    // {
    //     static_assert(std::is_same<T, BigDataStatMeth::hdf5Dataset >::value || 
    //                   std::is_same<T, BigDataStatMeth::hdf5DatasetInternal >::value,
    //                   "Error - type not allowed");
    // 
    //     try{
    //         
    //         H5::Exception::dontPrint();
    //         
    //         std::vector<hsize_t> stride = {1, 1},
    //                              block = {1, 1},
    //                              offset = {0, 0},
    //                              count = {0, 0};
    //         
    //         std::string stroutdataset = dsJoined->getDatasetName();
    //         
    //         
    //         
    //         BigDataStatMeth::hdf5Dataset* dstoJoin = new hdf5Dataset(dsJoined->getFileName(), strsubgroup, strinput[0], false);
    //         dstoJoin->openDataset();
    //         
    //         hsize_t* dims_out = dstoJoin->dim();
    //         
    //         //..ORIGINAL..// DataSet dataset = file->openDataSet(strsubgroup + "/" + strinput[0]);
    //         //..ORIGINAL..// IntegerVector dims_out = get_HDF5_dataset_size(dataset);
    //         
    //         //..ORIGINAL..// // Create unlimited dataset and add first dataset
    //         //..ORIGINAL..// if( exists_HDF5_element_ptr(file,stroutdataset))
    //         //..ORIGINAL..// remove_HDF5_element_ptr(file,stroutdataset);
    //         
    //         //..ORIGINAL..// create_HDF5_unlimited_matrix_dataset_ptr( file, stroutdataset, (unsigned long long)dims_out[0], (unsigned long long)dims_out[1], "real");
    //         //..ORIGINAL..// NumericMatrix readeddata(dims_out[0], dims_out[1]);
    //         //..ORIGINAL..// // read first dataset
    //         //..ORIGINAL..// read_HDF5_matrix_subset( file, &dataset, offset, dims_out, stride, block, REAL(readeddata) );
    //         //..ORIGINAL..// dataset.close();  
    //         
    //         //..ORIGINAL..// DataSet* unlimDataset = new DataSet(file->openDataSet(stroutdataset));
    //         
    //         // Add rows and needed cols to add the merged data in the new dataset
    //         Rcpp::Rcout<<"Extenem: "<<dims_out[0]<<" x "<<dims_out[1]<<"\n";
    //         dsJoined->extendUnlimitedDataset( (unsigned long long)dims_out[0], (unsigned long long)dims_out[1] );
    // 
    //         // Read data to merge
    //         std::vector<double> vreadeddata( dims_out[0] * dims_out[1] ); 
    //         dstoJoin->readDatasetBlock( {offset[0], offset[1]}, {dims_out[0], dims_out[1]}, stride, block, vreadeddata.data() );
    //         delete dstoJoin;
    //         
    //         // Write data to the new dataset
    //         dsJoined->writeDatasetBlock( Rcpp::wrap(vreadeddata), offset, dims_out, stride, block, false);
    // 
    //         //..ORIGINAL..// // write first dataset to new dataset
    //         //..ORIGINAL..// write_HDF5_matrix_subset_v2(file, unlimDataset, offset, dims_out, stride, block, transpose(readeddata));  
    //         
    //         // Update offset to new position
    //         offset[1] = offset[1] + dims_out[1];
    //         
    //         for( int i=1; i<strinput.size(); i++)
    //         {
    //             //..ORIGINAL..// dataset = file->openDataSet(strsubgroup + "/" + strinput[i]);
    //             //..ORIGINAL..// dims_out = get_HDF5_dataset_size(dataset);
    //             
    //             
    //             dstoJoin = new hdf5Dataset(dsJoined->getFileName(), strsubgroup, strinput[i], false);
    //             dstoJoin->openDataset();
    //             
    //             dims_out = dstoJoin->dim();
    //             
    //             // Update columns size
    //             //..ORIGINAL..// count = dims_out;
    //             
    //             // Extend dataset before put data
    //             //..ORIGINAL..// extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, dims_out[1]);
    //             dsJoined->extendUnlimitedDataset( (unsigned long long)0, (unsigned long long)dims_out[1] );
    //             
    //             //..ORIGINAL..// // Read data
    //             //..ORIGINAL..// NumericMatrix newdata(dims_out[0], dims_out[1]);
    //             //..ORIGINAL..// IntegerVector offsetRead = IntegerVector::create(0, 0);
    //             
    //             //..ORIGINAL..// // read dataset and close it.
    //             //..ORIGINAL..// read_HDF5_matrix_subset( file, &dataset, offsetRead, dims_out, stride, block, REAL(newdata) );
    //             //..ORIGINAL..// dataset.close();
    //             
    //             //..ORIGINAL..// write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, transpose(newdata)  );  
    //             
    //             // Read data to merge
    //             std::vector<double> vreadeddata( dims_out[0] * dims_out[1] ); 
    //             dstoJoin->readDatasetBlock( {offset[0], offset[1]}, {dims_out[0], dims_out[1]}, stride, block, vreadeddata.data() );
    //             delete dstoJoin;
    //             
    //             // Write data to the new dataset
    //             dsJoined->writeDatasetBlock( Rcpp::wrap(vreadeddata), offset, dims_out, stride, block, false);
    //             
    //             // Update offset
    //             offset[1] = offset[1] + dims_out[1];
    //         }
    //         
    //         if(bremoveJoined == true) {
    //             // Remove joined elements
    //             remove_elements(dsJoined->getDatasetptr(), strsubgroup, strinput);    
    //         }
    //         
    //         
    //     } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
    //         delete dsJoined;
    //         ::Rf_error( "c++ exception join_datasets (File IException)" );
    //         return -1;
    //     } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
    //         delete dsJoined;
    //         ::Rf_error( "c++ exception join_datasets (DataSet IException)" );
    //         return -1;
    //     } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
    //         delete dsJoined;
    //         ::Rf_error( "c++ exception join_datasets (Group IException)" );
    //         return -1;
    //     } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
    //         delete dsJoined;
    //         ::Rf_error( "c++ exception join_datasets (DataSpace IException)" );
    //         return -1;
    //     } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
    //         delete dsJoined;
    //         ::Rf_error( "c++ exception join_datasets (Data TypeIException)" );
    //         return -1;
    //     }
    //     return(0);
    // }


}

#endif // BIGDATASTATMETH_HDF5_UTILITIES_HPP

