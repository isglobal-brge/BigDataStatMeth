#ifndef hdf5_importFile
#define hdf5_importFile

   #include <RcppEigen.h>
   #include "rhdf5Utils.h"
   #include <regex>
   #include <iostream>
   #include <fstream>
   #include <vector>
   #include <algorithm>
   
   
   // Cpp functions
   std::vector<std::string> get_SplitData_in_vectorString(std::string line, std::regex reg_expres);
   bool is_number(const std::string& s);
   std::vector<double> get_data_as_Matrix(std::vector<std::string> strBlockValues);
   bool manage_Dataset( H5File* file, std::string outGroup, std::string outDataset, bool overwrite, int irows, int icols  );
   
   
   // R functions
   int Convert_text_to_HDF5( Rcpp::CharacterVector filename, 
                             std::string outputfile, std::string outGroup, std::string outDataset, 
                             Rcpp::Nullable<std::string> sep, 
                             Rcpp::Nullable<bool> header ,
                             Rcpp::Nullable<bool> rownames,
                             Rcpp::Nullable<bool> overwrite);
   
#endif