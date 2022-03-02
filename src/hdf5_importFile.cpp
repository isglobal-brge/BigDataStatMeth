#include "include/hdf5_importFile.h"

using namespace Rcpp;

// Split the line in several fields and store results in a Vector String.
std::vector<std::string> get_SplitData_in_vectorString(std::string line, std::regex reg_expres)
{
   std::vector<std::string> strValues;
   
   // Split line in columns by delim (reg_xpres)
   std::sregex_token_iterator 
   begin(line.begin(), line.end(), reg_expres),
   end;
   
   // write all the words to strValues
   std::copy(begin, end, std::back_inserter(strValues));
   
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
 
std::vector<double> get_data_as_Matrix(std::vector<std::string> strBlockValues)
{
   
   std::vector<double> doubleVector(strBlockValues.size());
   
   std::transform(strBlockValues.begin(), strBlockValues.end(), doubleVector.begin(), [](const std::string& s)
   {
      if( is_number(s) ==1 ){
         return (std::stod(s));
      } else {
         // Takes in to account the possible eol for the last column
         if( is_number(s.substr(0,s.length()-1) ) ==1 ){
            return (std::stod(s));
         } else {
            Rcpp::Rcout<<"\nValue : "<<s<<"\n";
            stop("Error: Column is not numeric. Only numeric data is allowed");
         }
      }
   });
   
   return(doubleVector);

}


// If not exists, create group and dataset in hdf5 data file
bool manage_Dataset( H5File* file, std::string outGroup, std::string outDataset, bool overwrite, int irows, int icols  )
{
   int res;
   
   bool allOk = false;
   
   try {
      
      // Test group
      if( !exists_HDF5_element_ptr(file, outGroup)) {
         res = create_HDF5_group_ptr(file, outGroup );
      }
      
      // Test dataset
      if( !exists_HDF5_element_ptr(file, outGroup + "/" + outDataset)) {
         create_HDF5_dataset_ptr(file, outGroup + "/" + outDataset, icols, irows, "double" );
         
      } else {
         if( overwrite == true ) {
            
            remove_HDF5_element_ptr(file, outGroup + "/" + outDataset);
            create_HDF5_dataset_ptr(file, outGroup + "/" + outDataset, icols, irows, "double" );
         } else {
            return(false);
         }
      }

      // Test dimnames and colnames for dataset
      std::string strGroupDimnames = "." + outDataset + "_dimnames";
      StringVector dimdatasets;
      dimdatasets.push_back("1");
      dimdatasets.push_back("2");
      
      if( exists_HDF5_element_ptr(file, outGroup + "/" + strGroupDimnames)) {
         remove_HDF5_multiple_elements_ptr(file, outGroup + "/" + strGroupDimnames, dimdatasets );
      }
      
   } catch( FileIException& error ) { // catch failure caused by the H5File operations
      ::Rf_error( "c++ exception Convert_text_to_HDF5 (File IException)" );
      return(allOk);
   } catch( GroupIException& error ) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception Convert_text_to_HDF5 (Group IException)" );
      return(allOk);
   } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
      ::Rf_error( "c++ exception Convert_text_to_HDF5 (DataSet IException)" );
      return(allOk);
   } catch(const std::runtime_error& re) {
      // speciffic handling for runtime_error
      ::Rcerr << "Runtime error: " << re.what() << std::endl;
      return(allOk);
   } catch(const std::exception& ex) {
      // speciffic handling for all exceptions extending std::exception, except
      // std::runtime_error which is handled explicitly
      ::Rcerr << "Error occurred: " << ex.what() << std::endl;
      return(allOk);
   } catch(...) {
      // catch any other errors (that we have no information about)
      ::Rcerr << "Unknown failure occurred. Possible memory corruption" << std::endl;
      return(allOk);
   }

   return(true);
   
}


//' Converts text file to hdf5 data file
//' 
//' Converts text file to hdf5 data file
//' 
//' @param filename string file name with data to be imported
//' @param outputfile file name and path to store imported data
//' @param outGroup group name to store the dataset
//' @param outDataset dataset name to store the input file in hdf5
//' @param sep (optional), by default = "\\t". The field separator string. Values within each row of x are separated by this string.
//' @param header (optional) either a logical value indicating whether the column names of x are to be written along with x, or a character vector of column names to be written. See the section on ‘CSV files’ for the meaning of col.names = NA.
//' @param rownames (optional) either a logical value indicating whether the row names of x are to be written along with x, or a character vector of row names to be written.
//' @param overwrite (optional) either a logical value indicating whether the output file can be overwritten or not.
//' 
//' @return none value returned, data are stored in a dataset inside an hdf5 data file.
//' @export
// [[Rcpp::export]]
void bdImport_text_to_hdf5( Rcpp::CharacterVector filename, 
                           std::string outputfile, std::string outGroup, std::string outDataset, 
                           Rcpp::Nullable<std::string> sep = R_NilValue, 
                           Rcpp::Nullable<bool> header = false,
                           Rcpp::Nullable<bool> rownames = false,
                           Rcpp::Nullable<bool> overwrite = false)
{
   
   std::string path = as<std::string>(filename);
   std::string stdsep;
   int res;
   
   
   // Colnames and rownames
   //..// CharacterVector svrownames, svrcolnames;
   CharacterVector svrcolnames;

   // Blocks control
   double counter = 0;
   double blockCounter = 1000;
   
   // hdf5 variables
   H5File* file;
   DataSet* datasetOut;
   
   try{
      
      if(sep.isNull()){
          stdsep = "\t";
      }else {
          stdsep = as<std::string>(sep);
      }
      
      std::string delim = "[^" + stdsep  + "]+";
      std::regex reg_expres(delim);
      
      if(ResFileExist_filestream(path))
      {
         std::ifstream inFile(path.c_str()); //Opens the file. c_str is mandatory here so that ifstream accepts the string path
         
         std::string line;
         
         std::getline(inFile,line,'\n'); //skip the first line (col names in our case). Remove those lines if note necessary

         
         // Number of columns
         std::ptrdiff_t const icols(std::distance(
               std::sregex_iterator(line.begin(), line.end(), reg_expres),
               std::sregex_iterator()));

         int incols = icols;
         if(as<bool>(rownames) == true) {
            // Read next line and count number of columns again depending on how file is created we can have
            // one empty space for rownames or not, then colnames will be different (-1 difference)
            std::getline(inFile,line,'\n'); //skip the first line (col names in our case). Remove those lines if note necessary
             
            // Number of columns
            std::ptrdiff_t const icols2(std::distance(
                  std::sregex_iterator(line.begin(), line.end(), reg_expres),
                  std::sregex_iterator()));
            
            if(icols2 == icols){
               incols = icols-1; // Reduce in one the number of columns
            }else if ( icols == icols2 -1){
               incols = icols; 
            }else{
               warning("Number of columns and headers are different, review data");
            }
                  
         } 
         
         // Re-adjust block size
         if(incols < 100 ){
            blockCounter = 10000;
         } /*else if (incols > 1000 ){
            blockCounter = 1000;
         }*/
         
         // Get number of rows (+1 to take in to account the last line without \n)
         int irows = std::count(std::istreambuf_iterator<char>(inFile),
                                std::istreambuf_iterator<char>(), '\n') ;
         
         
         // Restore counter after read first line to get number of cols
         if( as<bool>(header)==false ){
            irows = irows + 1;
         }

         CharacterVector svrownames(irows);

         // Prepare hdf5 data file
         if( ! ResFileExist_filestream(outputfile)) {
            res = Create_hdf5_file(outputfile); 
         }
         
         file = new H5File( outputfile, H5F_ACC_RDWR );
         
         if( !manage_Dataset( file, outGroup, outDataset, as<bool>(overwrite), irows, incols ) ){
            stop("Dataset exists - please set overwrite = true if you want to overwrite the dataset");
         }
         
         datasetOut =  new DataSet(file->openDataSet(outGroup + "/" + outDataset));

         // Reset iterator to beginning
         inFile.clear(); 
         inFile.seekg(0);
         
         // Read again the first line if header = true
         line.clear();
         
         // If data contains header :  Store first row as a colnames and reads next line
         if(as<bool>(header) == true) {
            std::getline(inFile,line,'\n');
            svrcolnames = wrap( get_SplitData_in_vectorString(line, reg_expres));
            // Read next line
            line.clear();
            std::getline(inFile,line,'\n');
         }

         std::vector<std::string> strBlockValues;
         IntegerVector stride = {1,1};
         IntegerVector block = {1,1};
         IntegerVector count = {incols, irows};
         IntegerVector offset = {0,0};
         bool btowrite;
         std::vector<std::string> strValues;

         while( !inFile.eof()  )
         {
            
            //..// std::vector<double> numbers;
            std::stringstream is(line); // take the line into a stringstream
            
            btowrite = true;
            
            
            // Get splitted values
            boost::split(strValues, line, boost::is_any_of(delim), boost::token_compress_on);

            if(as<bool>(rownames) == true) {
               svrownames[counter] =  strValues.front(); 
               strValues.erase(strValues.begin());
            }
            
            
            // Concatenate Valutes to get a block with several rows
            std::move(strValues.begin(), strValues.end(), std::back_inserter(strBlockValues));
            
            // Empty vector
            //..// strValues.erase (strValues.begin(),strValues.end());
            strValues.clear();
            
            // Write block
            
            if( counter>0 && (int)counter % (int)blockCounter == 0)
            {
               
               // offset[1] = counter - blockCounter;
               // count[1] = blockCounter;
               
                
                count[1] = strBlockValues.size() / incols;
               
               
               std::vector<double> doubleVector = get_data_as_Matrix(strBlockValues);
               
               double *p = doubleVector.data();
               Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> resMat (p, incols, strBlockValues.size() / incols );
               write_HDF5_matrix_subset_v2(file, datasetOut, offset, count, stride, block, wrap(resMat));
               
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
         
         // if(counter - blockCounter <0){
         //    offset[1] = 0;
         // }else {
         //    offset[1] = (floor((counter-1) / blockCounter) * blockCounter);
         // }
         

         // if(irows - (floor(irows/blockCounter)*blockCounter) == 0 && btowrite == true){
         //    count[1] = blockCounter;
         // } else {
            // count[1] = irows - (floor(irows/blockCounter)*blockCounter);
            count[1] = strBlockValues.size() / incols;
         // }
         
         if((irows - (floor(irows/blockCounter)*blockCounter)>0 && strBlockValues.size()>0) || btowrite == true)
         {
            std::vector<double> doubleVector = get_data_as_Matrix(strBlockValues);
            double *p = doubleVector.data();

            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> resMat (p, incols, count[1] );
            write_HDF5_matrix_subset_v2(file, datasetOut, offset, count, stride, block, wrap(resMat));

         }
         
         if(as<bool>(rownames) == true || as<bool>(header) == true ) {
             if(as<bool>(rownames) == false){
                Rcpp::StringVector svrownames(1);
                write_hdf5_matrix_dimnames(file, outGroup, outDataset, svrownames, svrcolnames );
             } else if(as<bool>(header) == false){
                 Rcpp::StringVector svrcolnames(1);
                 write_hdf5_matrix_dimnames(file, outGroup, outDataset, svrownames, svrcolnames );
             } else {
                 // Write rownames and colnames
                 write_hdf5_matrix_dimnames(file, outGroup, outDataset, svrownames, svrcolnames );     
             }
         }
          
         datasetOut->close();
         file->close();
         
      }else {
         warning("File doesn't exists, please, review location");
      }
      
   } catch( FileIException& error ) { // catch failure caused by the H5File operations
      datasetOut->close();
      file->close();
      ::Rf_error( "c++ exception Import_text_to_hdf5 (File IException)" );
      return void();
   } catch( GroupIException& error ) { // catch failure caused by the DataSet operations
      datasetOut->close();
       file->close();
      ::Rf_error( "c++ exception Import_text_to_hdf5 (Group IException)" );
      return void();
   } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
      datasetOut->close();
       file->close();
      ::Rf_error( "c++ exception Import_text_to_hdf5 (DataSet IException)" );
      return void();
   } catch(const std::runtime_error& re) {
      // speciffic handling for runtime_error
      datasetOut->close();
       file->close();
      ::Rcerr << "Runtime error: " << re.what() << std::endl;
      return void();
   } catch(const std::exception& ex) {
      // datasetOut->close();
      file->close();
      ::Rcerr << "Error occurred: " << ex.what() << std::endl;
      return void();
   } catch(...) {
      // catch any other errors (that we have no information about)
      datasetOut->close();
       file->close();
      ::Rcerr << "Unknown failure occurred. Possible memory corruption" << std::endl;
      return void();
   }
   
   Rcpp::Rcout<<"\n"<<outDataset<<" dataset has been imported successfully.\n";
   return void();
   //..// return(0);
   
}




/*** R
*/
