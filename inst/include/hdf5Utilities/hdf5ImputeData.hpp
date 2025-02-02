#ifndef BIGDATASTATMETH_UTIL_IMPUTE_DATA_HPP
#define BIGDATASTATMETH_UTIL_IMPUTE_DATA_HPP

#include "Utilities/openme-utils.hpp"
#include<random>


namespace BigDataStatMeth {


    // Get value for imputation
    int get_value_to_impute_discrete(std::map<double, double> probMap)
    {
        std::vector <double> probs;
        
        // Get values and counts for each map element
        for( auto it = probMap.begin(); it != probMap.end(); ++it )
            probs.push_back( it->second );
        
        // remove last element (corresponds to 3=<NA>)
        probs.erase(probs.end() - 1);
        
        // Get total count
        double totalSNPS = std::accumulate(probs.begin(), probs.end(), decltype(probs)::value_type(0));
        
        // Get probabilities without <NA>
        for (std::vector<double>::iterator it = probs.begin() ; it != probs.end(); ++it)
            *it = *it/totalSNPS;
        
        // Generate value with given probabilities
        std::random_device rd;
        std::mt19937 gen(rd());
        
        std::discrete_distribution<> d(probs.begin(), probs.end());
        
        return (d(gen));
        
    }
    
    
    // Convert NumericVector to map (key:vlues - value: frequency value in vector)
    std::map<double, double> VectortoOrderedMap_SNP_counts( Eigen::VectorXd  vdata)
    {
        std::map<double, double> mapv;
        
        try 
        {
            int position = 0;
            std::vector<double> v(vdata.data(), vdata.data()+vdata.size());
            
            std::sort(v.begin(), v.end() ); // Sort vector to optimize search and count
            
            for (size_t i = 0; i <=  *std::max_element(v.begin(), v.end()) ; ++i)  
            {
                double mycount = std::count(v.begin() + position, v.end(), i);
                mapv[i] = mycount;
                position = position + mycount;
            }
            
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...) { 
            ::Rf_error("c++ exception (unknown reason)"); 
        }
        
        return mapv;
    }
    
    
    // Pedestrian dataset imputation .... 
    // TODO : 
    //    - perform better imputation
    void Rcpp_Impute_snps_hdf5(BigDataStatMeth::hdf5Dataset* dsIn, BigDataStatMeth::hdf5DatasetInternal* dsOut,
                         bool bycols, std::string stroutdataset, Rcpp::Nullable<int> threads  = R_NilValue)
    {
        
        // DataSet* outdataset = nullptr;
        
        try{
            
            std::vector<hsize_t> offset = {0,0},
                count = {0,0},
                stride = {1,1},
                block = {1,1};
            
            int ilimit,
                blocksize = 1000;
            
            hsize_t* dims_out = dsIn->dim();
            
            // id bycols == true : read all rows by group of columns ; else : all columns by group of rows
            if (bycols == true) {
                ilimit = dims_out[0];
                count[1] = dims_out[1];
                offset[1] = 0;
            } else {
                ilimit = dims_out[1];
                count[0] = dims_out[0];
                offset[0] = 0;
            };
            
            
            if( stroutdataset.compare("")!=0) {
                dsOut->createDataset(dims_out[0], dims_out[1], "real");
            } 
            
            dsOut->openDataset();
            
            
            int ithreads = get_number_threads(threads, R_NilValue);
            int chunks = (ilimit/blocksize)/ithreads;
            
            #pragma omp parallel num_threads(ithreads) shared(dsIn, dsOut, chunks)
            {
                #pragma omp for schleude
                for( int i=0; i<=(ilimit/blocksize); i++) 
                {
                    int iread;
                    
                    if( (i+1)*blocksize < ilimit) iread = blocksize;
                    else iread = ilimit - (i*blocksize);
                    
                    if(bycols == true) {
                        count[0] = iread; 
                        offset[0] = i*blocksize;
                    } else {
                        count[1] = iread; 
                        offset[1] = i*blocksize;
                    }
                    
                    // read block
                    std::vector<double> vdIn( count[0] * count[1] ); 
                    #pragma omp critical(accessFile)
                    {
                        dsIn->readDatasetBlock( { offset[0], offset[1]}, { count[0], count[1]}, stride, block, vdIn.data() );
                    }
                    Eigen::MatrixXd data = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vdIn.data(), count[0], count[1] );
                    
                    if(bycols == true) // We have to do it by rows
                    {
                        for( int row = 0; row<data.rows(); row++)  // COMPLETE EXECUTION
                        {
                            std::map<double, double> myMap;
                            myMap = VectortoOrderedMap_SNP_counts(data.row(row));
                            
                            Eigen::VectorXd ev = data.row(row);
                            std::vector<double> v(ev.data(), ev.data() + ev.size());
                            
                            auto it = std::find_if(std::begin(v), std::end(v), [](int i){return i == 3;});
                            while (it != std::end(v)) {
                                
                                if(*it==3) *it = get_value_to_impute_discrete(myMap);
                                it = std::find_if(std::next(it), std::end(v), [](int i){return i == 3;});
                            }
                            
                            Eigen::VectorXd X = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
                            data.row(row) = X;
                            
                        }
                        
                    } else {
                        for( int col = 0; col<data.cols(); col++) 
                        {
                            std::map<double, double> myMap;
                            myMap = VectortoOrderedMap_SNP_counts(data.col(col));
                            
                            Eigen::VectorXd ev = data.col(col);
                            std::vector<double> v(ev.data(), ev.data() + ev.size());
                            
                            auto it = std::find_if(std::begin(v), std::end(v), [](int i){return i == 3;});
                            while (it != std::end(v)) {
                                
                                if(*it==3) *it = get_value_to_impute_discrete(myMap);
                                it = std::find_if(std::next(it), std::end(v), [](int i){return i == 3;});
                            }
                            
                            Eigen::VectorXd X = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
                            data.col(col) = X;
                        }
                    }
                    
                    #pragma omp critical(accessFile)
                    {
                        dsOut->writeDatasetBlock( Rcpp::wrap(data), offset, count, stride, block, false);
                    }
                }
                
            }
            
            
            
            // outdataset->close();
            
        } catch( H5::FileIException& error) { // catch failure caused by the H5File operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Impute_snps_hdf5 (File IException)" );
        } catch( H5::DataSetIException& error) { // catch failure caused by the DataSet operations
            delete dsIn;
            delete dsOut;;
            ::Rf_error( "c++ exception Rcpp_Impute_snps_hdf5 (DataSet IException)" );
        } catch( H5::GroupIException& error) { // catch failure caused by the Group operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Impute_snps_hdf5 (Group IException)" );
        } catch( H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Impute_snps_hdf5 (DataSpace IException)" );
        } catch( H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
            delete dsIn;
            delete dsOut;
            ::Rf_error( "c++ exception Rcpp_Impute_snps_hdf5 (Data TypeIException)" );
        }
        
        return void();
    }
    


}


#endif // BIGDATASTATMETH_UTIL_IMPUTE_DATA_HPP