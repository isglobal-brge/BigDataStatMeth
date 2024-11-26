#include <BigDataStatMeth.hpp>
#include "hdf5Utilities/hdf5ImportFiles.hpp"


//' Converts text file to hdf5 data file
//'
//' Converts text file to hdf5 data file
//'
//' @param filename string file name with data to be imported
//' @param outputfile file name and path to store imported data
//' @param outGroup group name to store the dataset
//' @param outDataset dataset name to store the input file in hdf5
//' @param sep (optional), by default = "\\t". The field separator string. 
//' Values within each row of x are separated by this string.
//' @param header (optional) either a logical value indicating whether the 
//' column names of x are to be written along with x, or a character vector of 
//' column names to be written. See the section on ‘CSV files’ for the meaning 
//' of col.names = NA.
//' @param rownames (optional) either a logical value indicating whether the 
//' row names of x are to be written along with x, or a character vector of 
//' row names to be written.
//' @param overwrite (optional) either a logical value indicating whether the 
//' output file can be overwritten or not.
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel 
//' computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent 
//' threads in parallelization if threads is null then threads =  maximum 
//' number of threads available
//'
//' @return none value returned, data are stored in a dataset inside an hdf5 data file.
//' @export
// [[Rcpp::export]]
void bdImportTextFile_hdf5( std::string filename,
                            std::string outputfile, std::string outGroup, std::string outDataset,
                            Rcpp::Nullable<std::string> sep = R_NilValue,
                            Rcpp::Nullable<bool> header = false,
                            Rcpp::Nullable<bool> rownames = false,
                            Rcpp::Nullable<bool> overwrite = false,
                            Rcpp::Nullable<bool> paral = R_NilValue,
                            Rcpp::Nullable<int> threads = R_NilValue)
{

    BigDataStatMeth::hdf5File* objFile;
    BigDataStatMeth::hdf5Dataset* datasetOut;

    try{
        
        bool boverwrite;
        
        if( overwrite.isNull()) { boverwrite = false; 
        } else { boverwrite = Rcpp::as<bool> (overwrite); }
        
        // Check if exists file to import
        if( BigDataStatMeth::Rcpp_FileExist(filename) ) {

            // Create file if does not exists
            objFile = new BigDataStatMeth::hdf5File(outputfile, false);
            int iRes = objFile->createFile();

            if(iRes == EXEC_WARNING) {
                objFile->openFile("rw");
            }
            
            // Create dataset
            datasetOut = new BigDataStatMeth::hdf5Dataset(objFile, outGroup, outDataset, boverwrite);
            // datasetOut->openDataset();

            Rcpp_Import_File_to_hdf5( filename, datasetOut, sep, header, rownames, paral, threads) ;

        } else {
            Rcpp::Rcerr << "File "<<filename<<" doesn't exists, please, review location" << std::endl;
        }

        
        
    }catch(const std::runtime_error& re) {
        // speciffic handling for runtime_error
        delete datasetOut;
        Rcpp::Rcerr << "c++ exception bdImportTextFile_hdf5 - Runtime error: " << re.what() << std::endl;
        return void();
    } catch(const std::exception& ex) {
        delete datasetOut;
        Rcpp::Rcerr << "c++ exception bdImportTextFile_hdf5 - Error occurred: " << ex.what() << std::endl;
        return void();
    } catch(...) {
        // catch any other errors (that we have no information about)
        delete datasetOut;
        Rcpp::Rcerr << "c++ exception bdImportTextFile_hdf5 - Unknown failure occurred. Possible memory corruption" << std::endl;
        return void();
    }
    
    Rcpp::message(Rcpp::wrap("The file has been imported"));
    
    delete datasetOut;
    delete objFile;
    return void();

}