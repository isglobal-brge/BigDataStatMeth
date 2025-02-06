#include <BigDataStatMeth.hpp>
#include "hdf5Omics/hdf5RemoveMAF.hpp"





//' Remove SNPs in hdf5 omic dataset with low data
//'
//' Remove SNPs in hdf5 omic dataset with low data
//' 
//' @param filename, character array indicating the name of the file to create
//' @param group, character array indicating the input group where the data set to be imputed is. 
//' @param dataset, character array indicating the input dataset to be imputed
//' @param outgroup, character array indicating group where the data set will be saved after remove data with if `outgroup` is NULL, output dataset is stored in the same input group. 
//' @param outdataset, character array indicating dataset to store the resulting data after imputation if `outdataset` is NULL, input dataset will be overwritten. 
//' @param maf, by default maf = 0.05. Numeric indicating the percentage to be considered to remove SNPs, SNPS with higest MAF will be removed from data
//' @param bycols, boolean by default = true, if true, indicates that SNPs are in cols, if SNPincols = false indicates that SNPs are in rows.
//' @param blocksize, integer, block size dataset to read/write and calculate MAF, by default this operations is made in with 100 rows if byrows = true or 100 cols if byrows = false.
//' @param overwrite, optional boolean if true, previous results in same location 
//' inside hdf5 will be overwritten, by default overwrite = false, data was not overwritten.
//' @return Original hdf5 data file with imputed data
//' @export
// [[Rcpp::export]]
void bdRemoveMAF_hdf5( std::string filename, std::string group, std::string dataset, 
                       std::string outgroup, std::string outdataset, 
                       Rcpp::Nullable<double> maf, Rcpp::Nullable<bool> bycols, 
                       Rcpp::Nullable<int> blocksize, Rcpp::Nullable<bool> overwrite = R_NilValue )
{
    
    // H5File* file = nullptr;
    int iremoved = 0;
    
    BigDataStatMeth::hdf5Dataset* dsIn = nullptr;
    BigDataStatMeth::hdf5DatasetInternal* dsOut = nullptr;
    
    try
    {
        
        bool bcols, bforce;
        double dpcent;
        int iblocksize = 100;
        
        std::string stroutdata = outgroup +"/" + outdataset;
        std::string strdataset = group +"/" + dataset;
        
        if(bycols.isNull()){  
            bcols = false ;
        }else{    
            bcols = Rcpp::as<bool>(bycols);
        }
        
        if(maf.isNull()){ dpcent = 0.05 ; } 
        else { dpcent = Rcpp::as<double>(maf); }
        
        if(!blocksize.isNull()){  
            iblocksize = Rcpp::as<int>(blocksize);
        }
        
        if(overwrite.isNull()) { bforce = false ; }
        else { bforce = Rcpp::as<bool>(overwrite); }
        
            
        if( strdataset.compare(stroutdata)!= 0)
        {
            dsIn = new BigDataStatMeth::hdf5Dataset(filename, group, dataset, false);
            dsIn->openDataset();
            
            dsOut = new BigDataStatMeth::hdf5DatasetInternal(filename, outgroup, outdataset, bforce);
            
        } else {
            throw std::range_error("Input and output dataset must be different");  
            return void();
        }
        
        iremoved = Rcpp_Remove_MAF_hdf5( dsIn, dsOut, bcols, dpcent, iblocksize);
        
        delete dsIn; dsIn = nullptr;
        delete dsOut; dsOut = nullptr;
        
        Rcpp::Function warning("warning");
        if (!bcols )
            warning( std::to_string(iremoved) + " Rows have been removed");
        else
            warning( std::to_string(iremoved) + " Columns have been removed");
        
    } catch( H5::FileIException& error ){
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdRemoveMAF_hdf5 (File IException)\n";
        return void();
    } catch( H5::DataSetIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdRemoveMAF_hdf5 (DataSet IException)\n";
        return void();
    } catch( H5::DataSpaceIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdRemoveMAF_hdf5 (DataSpace IException)\n";
        return void();
    } catch( H5::DataTypeIException& error ) { 
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdRemoveMAF_hdf5 (DataType IException)\n";
        return void();
    } catch(std::exception &ex) {
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nc++ c++ exception bdRemoveMAF_hdf5: "<< ex.what()<<"\n";
        return void();
    }  catch (...) {
        checkClose_file(dsIn, dsOut);
        Rcpp::Rcerr<<"\nC++ exception bdRemoveMAF_hdf5 (unknown reason)";
        return void();
    }
    
    return void();
    
}

