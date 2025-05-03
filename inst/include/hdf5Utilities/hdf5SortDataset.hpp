#ifndef BIGDATASTATMETH_UTIL_SORT_DATASETS_HPP
#define BIGDATASTATMETH_UTIL_SORT_DATASETS_HPP

namespace BigDataStatMeth {


    // Return index where value val is found in array
    template<class _InputIterator, class T>
    std::vector<_InputIterator>
    find_all(_InputIterator begin, _InputIterator end, const T& val)
    {
        std::vector<_InputIterator> matches;
        while(begin != end)
        {
            if((*begin) == val)
                matches.push_back(begin);
            ++begin;
        }
        
        return matches;
    }



    // Internal call 
    extern inline void RcppSort_dataset_hdf5 ( BigDataStatMeth::hdf5Dataset* dsIn, 
                                 BigDataStatMeth::hdf5Dataset* dsOut,
                                 Rcpp::List blockedSortlist, std::string func)
    {
        
        try {
            
            Rcpp::NumericVector oper = {0, 1};
            oper.names() = Rcpp::CharacterVector({ "sortRows", "sortCols"});
            
            Rcpp::StringVector rownames, colnames;
            std::vector<hsize_t> stride = {1, 1},
                                 block = {1, 1},
                                 offset = {0, 0},
                                 count = {0, 0};
            
            hsize_t* dims_out = dsIn->dim();
            
            for( int i = 0; i < blockedSortlist.length(); i++) {
                
                Rcpp::DataFrame df(blockedSortlist[i]);
                std::vector<double> order = df[1];
                std::vector<double> neworder = df[3];
                std::vector<double> diagonal = df[2];
                
                auto indices_0 = find_all(diagonal.begin(), diagonal.end(), 0);
                
                if( indices_0.size() > 0) {
                    
                    // for(int t=0; t<indices_0.size(); t++){
                    //     Rcpp::Rcout<<"Indices val : " <<&indices_0[t]<<"\n";    
                    // }
                    
                } else {
                    
                    if( oper.findName( func ) == 0 ) {
                        offset[0] = order[0] - 1;
                        count[0] = order.size();
                        count[1] = dims_out[1]; 
                        
                    } else if( oper.findName( func ) == 1 ) {
                        offset[1] = order[0] - 1;
                        count[1] = dims_out[1]; 
                        count[0] = order[order.size() - order[0]];
                    } 
                    
                    std::vector<double> vdIn( count[0] * count[1] ); 
                    dsIn->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdIn.data() );
                    
                    if( oper.findName( func ) == 0 ) {
                        offset[0] = neworder[0]-1;
                    } else if( oper.findName( func ) == 1 ) {
                        offset[1] = neworder[0]-1;
                    }
                    
                    dsOut->writeDatasetBlock(vdIn, offset, count, stride, block);
                    
                }
            }
            
        } catch( H5::FileIException& error ) {
            ::Rf_error( "c++ exception RcppSort_dataset_hdf5 (File IException )" );
            return void();
        } catch( H5::DataSetIException& error ) { // catch failure caused by the dstosplit operations
            ::Rf_error( "c++ exception RcppSort_dataset_hdf5 (dstosplit IException )" );
            return void();
        } catch( H5::DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
            ::Rf_error( "c++ exception RcppSort_dataset_hdf5 (DataSpace IException )" );
            return void();
        } 
        
        return void();
        
    }

}

#endif // BIGDATASTATMETH_UTIL_SORT_DATASETS_HPP