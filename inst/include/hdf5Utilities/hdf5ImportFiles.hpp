#ifndef BIGDATASTATMETH_HDF5_IMPORTFILES_HPP
#define BIGDATASTATMETH_HDF5_IMPORTFILES_HPP


// #include <fstream>
// #include <boost/algorithm/string.hpp>

namespace BigDataStatMeth {

    using matrix = std::vector< std::vector<double> >;


    extern inline matrix transpose( const matrix &M )
    {
        int rows = M.size();
        int cols = M[0].size();
        
        matrix T( cols, std::vector<double>( rows ) );
        
        for ( int i = 0; i < rows; i++ )
        {
            for ( int j = 0; j < cols; j++ ) T[j][i] = M[i][j];
        }
        
        return T;
    }

    // Split the line in several fields and store results in a Vector String.
    extern inline std::vector<std::string> get_SplitData_in_vectorString(std::string line, std::regex reg_expres)
    {
        std::vector<std::string> strValues;

        // Split line in columns by delim (reg_xpres)
        std::sregex_token_iterator
        begin(line.begin(), line.end(), reg_expres),
        end;

        // write all the words to strValues
        std::copy(begin, end, std::back_inserter(strValues));
        
        // for (auto i: strValues)
        //     std::cout << i << ' ';

        return(strValues);
    }


    // Test if read data is numeric or not
    bool is_number(const std::string& s)
    {
        char* end = nullptr;
        double val = strtod(s.c_str(), &end);

        return end != s.c_str() && *end == '\0' && val != HUGE_VAL;
    }


    // Convert readed block from string to double

    extern inline std::vector<double> get_data_as_Matrix(std::vector<std::string> strBlockValues)
    {
        
        std::vector<double> doubleVector(strBlockValues.size());
        std::transform(strBlockValues.begin(), strBlockValues.end(), doubleVector.begin(), [](const std::string& s) mutable
        {
            if( is_number(s) ==1 ){
                return (std::stod(s));
            } else {
                // Takes in to account the possible eol for the last column
                if( is_number(s.substr(0,s.length()-1) ) ==1 ){
                    return (std::stod(s));
                } else {
                    Rcpp::stop("Error: Column is not numeric. Only numeric data is allowed. Maybe there is a blank row at the end of the file or any field is empty");
                }
            }
        });

        return(doubleVector);

    }




    extern inline void Rcpp_Import_File_to_hdf5( Rcpp::CharacterVector filename,
                                   BigDataStatMeth::hdf5Dataset* dsOut,
                                   Rcpp::Nullable<std::string> sep = R_NilValue,
                                   Rcpp::Nullable<bool> header = false,
                                   Rcpp::Nullable<bool> rownames = false,
                                   Rcpp::Nullable<bool> bparal = R_NilValue,
                                   Rcpp::Nullable<int> threads = R_NilValue)
    {

        try {

            std::string path = Rcpp::as<std::string>(filename);
            std::string stdsep;

            // Colnames and rownames
            Rcpp::CharacterVector svrcolnames;

            // Blocks control
            double counter = 0;
            double blockCounter = 1000;


            if(sep.isNull()){
                stdsep = "\t";
            }else {
                stdsep = Rcpp::as<std::string>(sep);
            }

            std::string delim = "[^" + stdsep  + "]+";
            std::regex reg_expres(delim);


            std::string line;

            std::ifstream inFile(path.c_str()); //Opens the file. c_str is mandatory here so that ifstream accepts the string path
            std::getline(inFile,line,'\n'); //skip the first line (col names in our case). Remove those lines if note necessary

            // Number of columns
            std::ptrdiff_t const icols(std::distance(
                    std::sregex_iterator(line.begin(), line.end(), reg_expres),
                    std::sregex_iterator()));

            hsize_t incols = icols;

            if(Rcpp::as<bool>(rownames) == true) {
                // Read next line and count number of columns again depending on how file is created we can have
                // one empty space for rownames or not, then colnames will be different (-1 difference)
                std::getline(inFile,line,'\n'); //skip the first line (col names in our case). Remove those lines if not necessary

                // Number of columns
                std::ptrdiff_t const icols2(std::distance(
                        std::sregex_iterator(line.begin(), line.end(), reg_expres),
                        std::sregex_iterator()));

                if(icols2 == icols){
                    incols = icols-1; // Reduce in one the number of columns
                } else if ( icols == icols2 -1){
                    incols = icols;
                } else {
                    Rcpp::stop("Number of columns and headers are different, please review data, note that fields without values are not allowed");
                    // Rcpp::warning("Number of columns and headers are different, review data");
                    
                }
            }
            
            // Re-adjust block size
            if(incols < 100 ){
                blockCounter = 10000;
            }

            // Get number of rows (+1 to take in to account the last line without \n)
            int irows = std::count(std::istreambuf_iterator<char>(inFile),
                                   std::istreambuf_iterator<char>(), '\n') + 1 ;

            // Restore counter after read first line to get number of cols
            if( Rcpp::as<bool>(header)==false ){
                irows = irows + 1;
            }

            Rcpp::CharacterVector svrownames(irows);

            // Reset iterator to beginning
            inFile.clear();
            inFile.seekg(0);

            // Read again the first line if header = true
            line.clear();

            // If data contains header :  Store first row as a colnames and reads next line
            if(Rcpp::as<bool>(header) == true) {
                std::getline(inFile,line,'\n');
                svrcolnames = Rcpp::wrap( get_SplitData_in_vectorString(line, reg_expres));
                // If rownames then remove first column from header (belonging to the rownames)
                if(Rcpp::as<bool>(rownames) == true) {
                    if( incols ==  svrcolnames.size() || incols ==  (svrcolnames.size()-1)){
                        svrcolnames.erase(0);}
                }
                // Read next line
                line.clear();
                std::getline(inFile,line,'\n');
            }

            dsOut->createDataset( (hsize_t)irows, (hsize_t)incols, "real");
            
            std::vector<std::string> strBlockValues;
            std::vector<hsize_t> stride = {1,1},
                                 block = {1,1},
                                 count = { (hsize_t)incols, (hsize_t)irows},
                                 offset = {0,0};

            bool btowrite;
            std::vector<std::string> strValues;
            
            while( !inFile.eof()  )
            {

                std::stringstream is(line); // take the line into a stringstream

                btowrite = true;
                
                // Get splitted values
                boost::split(strValues, line, boost::is_any_of(delim), boost::token_compress_on);
                
                if( Rcpp::as<bool>(rownames) == true ) {
                    
                    svrownames[counter] =  strValues.front();
                    strValues.erase(strValues.begin());
                }

                // Concatenate Valutes to get a block with several rows
                std::move(strValues.begin(), strValues.end(), std::back_inserter(strBlockValues));

                // Empty vector
                strValues.clear();

                // Write block
                if( counter>0 && (int)counter % (int)blockCounter == 0)
                {
                    count[1] = strBlockValues.size() / incols;

                    std::vector<double> doubleVector = get_data_as_Matrix(strBlockValues);
                    
                    double *p = doubleVector.data();
                    
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> resMat (p, incols, strBlockValues.size() / incols );
                    dsOut-> writeDatasetBlock( Rcpp::wrap(resMat.transpose()), offset, count, stride, block, false);

                    offset[1] = offset[1] + (strBlockValues.size() / incols);
                    
                    // Clear Vector
                    strBlockValues.clear();
                    btowrite = false;

                }

                // Clear Buffer and Read next line
                line.clear();
                std::getline(inFile,line,'\n');

                // Increment counter
                counter++;

            }
            
            count[1] = strBlockValues.size() / incols;

            if((irows - (floor(irows/blockCounter)*blockCounter)>0 && strBlockValues.size()>0) || btowrite == true)
            {
                std::vector<double> doubleVector = get_data_as_Matrix(strBlockValues);
                
                double *p = doubleVector.data();
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> resMat (p, incols, count[1] );
                
                dsOut-> writeDatasetBlock( Rcpp::wrap(resMat.transpose()), offset, count, stride, block, false);
                    
            }

            
            BigDataStatMeth::hdf5Dims* dsdims;
            dsdims = new BigDataStatMeth::hdf5Dims(dsOut);
            
            if( Rcpp::as<bool>(rownames) == true || Rcpp::as<bool>(header) == true ) {
                if( Rcpp::as<bool>(rownames) == false){
                    Rcpp::StringVector svrownames(1);
                    dsdims->writeDimnames( Rcpp::wrap(svrownames), Rcpp::wrap(svrcolnames));
                } else if(Rcpp::as<bool>(header) == false){
                    Rcpp::StringVector svrcolnames(1);
                    dsdims->writeDimnames( svrownames, svrcolnames);
                } else {
                    // Write rownames and colnames
                    dsdims->writeDimnames( svrownames, svrcolnames);
                }
            }
            
            delete dsdims;

        } catch( H5::FileIException& error ) {
            ::Rf_error( "c++ exception Convert_text_to_HDF5 (File IException)" );
            return void();
        } catch( H5::GroupIException& error ) {
            ::Rf_error( "c++ exception Convert_text_to_HDF5 (Group IException)" );
            return void();
        } catch( H5::DataSetIException& error ) {
            ::Rf_error( "c++ exception Convert_text_to_HDF5 (DataSet IException)" );
            return void();
        } catch(const std::runtime_error& re) {
            Rcpp::Rcerr << "Runtime error: " << re.what() << std::endl;
            return void();
        } catch(const std::exception& ex) {
            Rcpp::Rcerr << "Error occurred: " << ex.what() << std::endl;
            return void();
        } catch(...) {
            Rcpp::Rcerr << "Unknown failure occurred. Possible memory corruption" << std::endl;
            return void();
        }

        return void();

    }



}

#endif // BIGDATASTATMETH_HDF5_IMPORTFILES_HPP

