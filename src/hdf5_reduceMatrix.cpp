#include "include/hdf5_splitMatrix.h"


// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;


// Internal call 
int RcppReduce_matrix_hdf5 ( H5File* file,  std::string strgroup, std::string stroutgroup, std::string stroutdataset, std::string strreducefunction, bool bremove )
{
    DataSet* ploaddataset = nullptr;
    Eigen::MatrixXd fullReduced;
    IntegerVector dims_out;
    int ndatasets;
    std::string strdatasetout = stroutgroup + "/" + stroutdataset;
    
    try {
    
        // Get dataset names without prefix, all datasets inside the group
        StringVector joindata =  get_dataset_names_from_group(file, strgroup, "");
    
        ndatasets = joindata.size();
        
        for ( int i=0; i< ndatasets; i++)
        {
            ploaddataset = new DataSet(file->openDataSet(strgroup + "/" + joindata[i]));
            dims_out = get_HDF5_dataset_size(*ploaddataset);
            
            if( i == 0 ) {
                fullReduced = GetCurrentBlock_hdf5_Original( file, ploaddataset, 0, 0, dims_out[0],dims_out[1]);
                
            } else {
                
                // If readed block is smaller than full matrix adds rows or columns
                Eigen::MatrixXd newRead = GetCurrentBlock_hdf5_Original( file, ploaddataset, 0, 0, dims_out[0],dims_out[1]);
                if( newRead.rows() != fullReduced.rows()){
                  int difference = fullReduced.rows() - newRead.rows();
                  newRead.conservativeResize(newRead.rows() + difference, newRead.cols()); 
                }
                
                if( newRead.cols() != fullReduced.cols()){
                  int difference = fullReduced.cols() - newRead.cols();
                  newRead.conservativeResize(newRead.rows(), newRead.cols() + difference); 
                }
                
                // Reduce matrix
                if( strreducefunction.compare("+")==0) {
                    fullReduced = fullReduced + newRead;
                } else if (strreducefunction.compare("-")==0) {
                    fullReduced = fullReduced - newRead;
                } 
            }
            
            ploaddataset->close();
            
            if( bremove == true){
                // Remove used dataset
                remove_HDF5_element_ptr(file, strgroup + "/" + joindata[i]);
            }
            
        }
        
      // Transform to rowmajor
        Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > mapBlock(fullReduced.data(), 
                                                                                                  fullReduced.cols() , fullReduced.rows());
        
        write_HDF5_matrix_from_R_ptr(file, strdatasetout, Rcpp::wrap(mapBlock.transpose()), false);
        
    }catch( FileIException& error ) {
        ploaddataset->close();
        file->close();
        ::Rf_error( "c++ exception (File IException )" );
        return -1;
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        ploaddataset->close();
        file->close();
        ::Rf_error( "c++ exception (DataSet IException )" );
        return -1;
    } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        ploaddataset->close();
        file->close();
        ::Rf_error( "c++ exception (DataSpace IException )" );
        return -1;
    } 
    
    return(0);

  
}





//' Reduce hdf5 dataset
//'
//' Reduce hdf5 datasets inside a group by rows or columns and store complete matrix inside hdf5 data file.
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data sets are stored 
//' @param reducefunction, single character with function to apply, can be '+' or  '-' array indicating the input dataset to be imputed
//' @param outgroup, optional character array indicating group where the data set will be saved after imputation if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, optional character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
//' @param force, boolean if true, previous results in same location inside hdf5 will be overwritten.
//' @param remove, boolean if true, removes original matrices, by default bremove = false.
//' @return Full matrix with results from reduction
//' @export
// [[Rcpp::export]]
Rcpp::RObject bdReduce_matrix_hdf5( std::string filename, std::string group, 
                                    std::string reducefunction, 
                                    Rcpp::Nullable<std::string> outgroup = R_NilValue, Rcpp::Nullable<std::string> outdataset = R_NilValue,
                                    Rcpp::Nullable<bool> force = false , Rcpp::Nullable<bool> remove = false )
{
  
  H5File* file;

  try
  {
    std::string stroutgroup,  
                strdatasetout,
                stroutdata;
                
    bool bforce, bremove;
    

    if (reducefunction.compare("+") != 0 && reducefunction.compare("-") != 0 ) {
      throw std::range_error( "Function to apply must be \"+\" or \"-\" other values are not allowed" );
      return(Rcpp::wrap(-1));
    }
    
    if(outgroup.isNull()) {  stroutgroup = group ;
    } else {   stroutgroup = Rcpp::as<std::string>(outgroup);}
    
    if(remove.isNull()) { bremove = false ;
    } else {   bremove = Rcpp::as<bool>(remove);}
    
    if(force.isNull()) { bforce = false; } 
    else {   bforce = Rcpp::as<bool>(force); }

    if(outdataset.isNull()){  strdatasetout = group ;
    } else {   strdatasetout = Rcpp::as<std::string>(outdataset);}

    if (exist_FileGroupDataset (filename, group, "")!= 0 ) {
      file = new H5File( filename, H5F_ACC_RDWR );
    }
    
    prepare_outGroup(file, stroutgroup, bforce);

    prepare_outDataset(file, stroutgroup + "/" + strdatasetout, bforce);

    RcppReduce_matrix_hdf5( file, group, stroutgroup, strdatasetout, reducefunction, bremove);
    
    
  }
  catch( FileIException& error ) { // catch failure caused by the H5File operations
    file->close();
    ::Rf_error( "c++ exception (File IException)" );
    return(wrap(-1));
  }
  
  file->close();
  
  return(wrap(0));
}



/***R

*/
