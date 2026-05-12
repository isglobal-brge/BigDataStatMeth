/**
 * @file hdf5Dims.hpp
 * @brief Dimension management utilities for HDF5 datasets
 * 
 * This file provides functionality for managing dimension names (row and column names)
 * in HDF5 datasets. It implements efficient storage and retrieval of dimension
 * information, with support for large datasets through block-wise processing.
 * 
 * Key features:
 * - Row and column name management
 * - Efficient string storage
 * - Block-wise processing for large datasets
 * - Automatic dimension validation
 * - Comprehensive error handling
 * 
 * @note This module is part of the BigDataStatMeth library
 */

#ifndef BIGDATASTATMETH_HDF5_DIMS_HPP
#define BIGDATASTATMETH_HDF5_DIMS_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

namespace BigDataStatMeth
{
    /**
     * @brief Class for managing dimension names in HDF5 datasets
     * 
     * This class extends hdf5Group to provide specialized functionality for
     * storing and managing dimension names (row and column names) associated
     * with HDF5 datasets.
     */
    class hdf5Dims : public hdf5Group
    {
    public:
        /**
         * @brief Constructs a dimension manager for an HDF5 dataset
         *
         * A single constructor serves both write mode (default) and read mode.
         * In write mode the existing dimnames group is removed and a fresh one
         * is created, matching the previous behaviour.  In read mode the group
         * is left untouched so that @c readDimnames() can access it.
         *
         * @param pdataset Pointer to the open hdf5Dataset to manage dimensions for.
         * @param bWrite   When @c true (default) the dimnames group is rebuilt for
         *                 writing.  Pass @c false to open in read-only mode.
         *
         * @note In write mode any pre-existing dimnames group is deleted first.
         */
        hdf5Dims(BigDataStatMeth::hdf5Dataset* pdataset, bool bWrite = true) :
        hdf5Group(pdataset->getFullPath(), pdataset->getGroup() + "/." + pdataset->getDatasetName() + "_dimnames")
        {
            pmaindataset = pdataset;

            if (bWrite) {
                if( exists_HDF5_element(pfile, groupname) ) {
                    remove_elements(pmaindataset->getFileptr(), pmaindataset->getGroup(),
                                    "." + pmaindataset->getDatasetName() + "_dimnames" );
                }
                create_HDF5_groups(groupname);
            }
        }

        /**
         * @brief Writes dimension names to the HDF5 file
         * 
         * @param rownames Vector of row names
         * @param colnames Vector of column names
         * 
         * @throws H5::FileIException on file operation errors
         * @throws H5::DataSetIException on dataset operation errors
         * @throws H5::GroupIException on group operation errors
         * @throws H5::DataSpaceIException on dataspace operation errors
         * @throws H5::DataTypeIException on datatype operation errors
         * 
         * @note Validates dimensions against the main dataset
         * @note Overwrites existing dimension names if they exist
         * 
         * Performance considerations:
         * - Uses block-wise processing for large string vectors
         * - Implements efficient string storage with fixed-length buffers
         * - Handles memory cleanup automatically
         */
        void writeDimnames( Rcpp::StringVector rownames, Rcpp::StringVector colnames)
        {
            try {
                // static_assert(std::is_same<T, Rcpp::RObject >::value,
                //               "Error - type not allowed");
                int nrows, ncols;

                if ( rownames.length() > 0 ) {
                    nrows = rownames.length();
                } else  {
                    nrows = 0;
                }

                if (  colnames.length() > 0 ) {
                    ncols = colnames.length();
                } else  {
                    ncols = 0;
                }
                
                if( nrows<1 && ncols<1) {
                    Rcpp::warning("Data not provided to write dimensions");
                } else {

                    // Write rownames
                    //. 20260222 .// if( (unsigned)nrows == pmaindataset->nrows_file()) {
                    if( (unsigned)nrows == pmaindataset->ncols()) {
                        
                        std::string fullDatasetPath = groupname + "/" + strrows;

                        dimrownames[0] = nrows; dimrownames[1] = 1;
                        H5::DataSpace dataspace( RANK2, dimrownames );
                        bool bexists = exists_HDF5_element(pfile, fullDatasetPath);
                        if( bexists == true ) {
                            Rcpp::warning ("Rownames already exits and will be overwritten");
                            remove_elements(pmaindataset->getFileptr(), groupname, {strrows} );
                        }
                        
                        pmaindataset->setRownamesDatasetPath( fullDatasetPath );
                        writeStringVector( pdsrownames, fullDatasetPath, rownames );
                    }

                    // Write colnames
                    //. 20260222 .// if( (unsigned)ncols == pmaindataset->ncols_file()) {
                    if( (unsigned)ncols == pmaindataset->nrows()) {
                        
                        std::string fullDatasetPath = groupname + "/" + strcols;

                        dimcolnames[0] = ncols; dimcolnames[1] = 1;

                        // H5::DataSpace dataspace( RANK2, dimcolnames );
                        bool bexists = exists_HDF5_element(pfile, fullDatasetPath);
                        if( bexists == true ) {
                            Rcpp::warning ("Rownames already exits and will be overwritten");
                            remove_elements(pmaindataset->getFileptr(), groupname, {strcols} );
                        }
                        pmaindataset->setColnamesDatasetPath( fullDatasetPath );
                        writeStringVector( pdscolnames, fullDatasetPath, colnames );

                    }
                }

            } catch(H5::FileIException& error) { // catch failure caused by the H5File operations
                close_datasets();
                throw std::runtime_error("c++ exception writeDimnames (File IException)");
            } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
                close_datasets();
                throw std::runtime_error("c++ exception writeDimnames (DataSet IException)");
            } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
                close_datasets();
                throw std::runtime_error("c++ exception writeDimnames (Group IException)");
            } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
                close_datasets();
                throw std::runtime_error("c++ exception writeDimnames (DataSpace IException)");
            } catch(H5::DataTypeIException& error) { // catch failure caused by the DataSpace operations
                close_datasets();
                throw std::runtime_error("c++ exception writeDimnames (Data TypeIException)");
            }

          return void();

        }


        // /**
        //  * @brief Reads dimension names stored alongside an HDF5 dataset
        //  *
        //  * Opens the hidden @c _dimnames group written by @c writeDimnames() and
        //  * reads back the row-names vector (@c /1) and the column-names vector
        //  * (@c /2).  If a vector does not exist an empty @c Rcpp::CharacterVector
        //  * is returned for that component.
        //  *
        //  * @return Named @c Rcpp::List with two elements:
        //  *   - @c "rownames" : @c CharacterVector of row names (length 0 if absent)
        //  *   - @c "colnames" : @c CharacterVector of col names (length 0 if absent)
        //  *
        //  * @pre The @c hdf5Dims object must have been constructed with
        //  *      @c bWrite = false (read mode).
        //  * @pre The underlying hdf5Dataset must be open.
        //  *
        //  * @throws std::runtime_error on HDF5 I/O failures
        //  */
        // Rcpp::List readDimnames()
        // {
        //     Rcpp::CharacterVector rn(0), cn(0);
        // 
        //     try {
        //         H5::Exception::dontPrint();
        // 
        //         const std::string path_rows = groupname + "/" + strrows;  // .ds_dimnames/1
        //         const std::string path_cols = groupname + "/" + strcols;  // .ds_dimnames/2
        // 
        //         rn = readStringDataset(path_rows);
        //         cn = readStringDataset(path_cols);
        // 
        //     } catch (H5::FileIException& error) {
        //         throw std::runtime_error("c++ exception readDimnames (File IException)");
        //     } catch (H5::DataSetIException& error) {
        //         throw std::runtime_error("c++ exception readDimnames (DataSet IException)");
        //     } catch (H5::GroupIException& error) {
        //         throw std::runtime_error("c++ exception readDimnames (Group IException)");
        //     } catch (std::exception& ex) {
        //         throw std::runtime_error(std::string("c++ exception readDimnames: ") + ex.what());
        //     }
        // 
        //     return Rcpp::List::create(Rcpp::Named("rownames") = rn,
        //                               Rcpp::Named("colnames") = cn);
        // }


        /**
         * @brief Virtual destructor
         */
        virtual ~hdf5Dims(){
        }

    protected:
        /**
         * @brief Structure for storing name strings
         */
        typedef struct names {
            char chr[MAXSTRING];  ///< Fixed-length character array for name storage
        } names;

        BigDataStatMeth::hdf5Dataset* pmaindataset;  ///< Pointer to main dataset
        H5::DataSet* pdsrownames = nullptr;          ///< Dataset for row names
        H5::DataSet* pdscolnames = nullptr;          ///< Dataset for column names
        std::string strrows = "1";                   ///< Identifier for row names
        std::string strcols = "2";                   ///< Identifier for column names
        hsize_t dimcolnames[2];                      ///< Dimensions for column names
        hsize_t dimrownames[2];                      ///< Dimensions for row names


        // /**
        //  * @brief Reads one string-CompType dataset from the dimnames group
        //  *
        //  * Opens the dataset at @p path inside the current file and reads it as
        //  * a @c CharacterVector using the same @c CompType layout used by
        //  * @c writeStringVector().  Returns an empty vector when @p path does
        //  * not exist.
        //  *
        //  * @param path  Full HDF5 path of the string dataset (e.g. group + "/1")
        //  * @return      @c Rcpp::CharacterVector with the stored strings, or
        //  *              @c CharacterVector(0) when the dataset is absent.
        //  */
        // Rcpp::CharacterVector readStringDataset(const std::string& path)
        // {
        //     if (!exists_HDF5_element(pfile, path))
        //         return Rcpp::CharacterVector(0);
        // 
        //     // Mirror the CompType used in writeStringVector()
        //     typedef struct name_t { char chr[MAXSTRING]; } name_t;
        // 
        //     H5::DataSet   hds  = pfile->openDataSet(path);
        //     H5::DataSpace spc  = hds.getSpace();
        // 
        //     hsize_t dims[1] = {0};
        //     spc.getSimpleExtentDims(dims);
        //     const hsize_t n = dims[0];
        // 
        //     if (n == 0) {
        //         hds.close();
        //         return Rcpp::CharacterVector(0);
        //     }
        // 
        //     H5::CompType mtype(sizeof(name_t));
        //     mtype.insertMember("chr", HOFFSET(name_t, chr),
        //                        H5::StrType(H5::PredType::C_S1, MAXSTRING));
        // 
        //     std::vector<name_t> buf(n);
        //     hds.read(buf.data(), mtype);
        //     hds.close();
        // 
        //     Rcpp::CharacterVector cv(n);
        //     for (hsize_t i = 0; i < n; ++i)
        //         cv[static_cast<int>(i)] = std::string(buf[i].chr);
        //     return cv;
        // }


        /**
         * @brief Closes all open datasets
         */
        void close_datasets()
        {
            pdsrownames->close();
            pdscolnames->close();

        }

        /**
         * @brief Converts R data frame to HDF5-compatible range list
         * 
         * @param DatasetValues R object containing values to convert
         * @param rowscols Identifier for row or column conversion
         * @param bFullDataset Flag indicating if full dataset conversion is needed
         * @return names* Pointer to array of converted names
         * 
         * @note Caller is responsible for freeing returned memory
         */
        names* convert_DataFrame_to_RangeList(Rcpp::RObject DatasetValues, std::string rowscols, bool bFullDataset)
        {

            int datarows = Rcpp::as<Rcpp::StringVector>(DatasetValues).size();
            int isizetoWrite = datarows;

            if(bFullDataset) {
                if( rowscols == "rows") {
                    isizetoWrite = dimrownames[0] * dimrownames[1];
                } else {
                    isizetoWrite = dimcolnames[0] * dimcolnames[1];
                }
            }

            names *names_list = new names[isizetoWrite];  // Convert to range list

            // Write data to dataset, if data to write is smaller than dataset then empty positions='\0'
            for( int i = 0; i < isizetoWrite; i++ ) {
                int j = 0;
                if(i< datarows) {
                    Rcpp::String wchrom = Rcpp::as<Rcpp::StringVector>(DatasetValues)(i);
                    std::string word = wchrom.get_cstring();

                    for( j = 0; (unsigned)j < word.size() && j < (MAXSTRING-1); j++ ) {
                        names_list[i].chr[j] = word[j];
                    }
                }
                names_list[i].chr[j] = '\0'; // insert hdf5 end of string
            }
            return(names_list);
        }

    private:
        /**
         * @brief Writes a string vector to an HDF5 dataset
         * 
         * @param dataset Pointer to HDF5 dataset
         * @param datasetname Name of the dataset
         * @param DatasetValues Vector of strings to write
         * 
         * Implementation details:
         * 1. Creates appropriate HDF5 datatype for strings
         * 2. Processes data in blocks for large vectors
         * 3. Handles string truncation and null termination
         * 4. Manages memory efficiently
         * 
         * @note Uses MAXSTRBLOCK for block size in processing
         * @note Automatically handles quotation mark removal
         */
        void writeStringVector( H5::DataSet* dataset, std::string datasetname, Rcpp::StringVector DatasetValues)
        {
            
            try
            {
                typedef struct name {
                    char chr[MAXSTRING];
                } name;
                
                // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
                H5::Exception::dontPrint();
                
                // Create the data space for the dataset.
                hsize_t vectorsize;
                
                if (Rcpp::is<Rcpp::StringVector>(DatasetValues))
                {
                    vectorsize = DatasetValues.length();
                    
                    // Define hdf5 dataspace size
                    hsize_t dims[] = {vectorsize};
                    H5::DataSpace dataspace(RANK1, dims);
                    
                    // Create the memory datatype.
                    H5::CompType mtype(sizeof(name));
                    mtype.insertMember("chr", HOFFSET(name, chr), H5::StrType(H5::PredType::C_S1, MAXSTRING ));
                    
                    // Create the dataset.
                    dataset = new H5::DataSet(pfile->createDataSet(datasetname, mtype, dataspace));
                    
                    //  Get dataset dataspace.
                    dataspace = dataset->getSpace();
                    
                    if(vectorsize > MAXSTRBLOCK) {
                        
                        // Number of blocks to process
                        int iblocsks = vectorsize/MAXSTRBLOCK;
                        
                        for(int i=0; i<=iblocsks; i++)
                        {
                            // Gets block size to read
                            hsize_t ilength = MAXSTRBLOCK;
                            if(i == iblocsks){
                                ilength = vectorsize - (i * MAXSTRBLOCK);
                            }
                            
                            // Convert Dataframe to range list
                            name *names_list = new name[ilength];
                            
                            for(int row=0; (unsigned)row< ilength; row++ )
                            {
                                Rcpp::String wchrom = Rcpp::as<Rcpp::StringVector>(DatasetValues)((i*MAXSTRBLOCK) + row);
                                std::string word = wchrom.get_cstring();
                                
                                //.. 20260423..//// boost::algorithm::erase_all(word, "\"");
                                word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
                                
                                int j=0;
                                for( j = 0; (unsigned)j < word.size() && j < (MAXSTRING-1); j++ ){
                                    names_list[row].chr[j] = word[j]; }
                                
                                names_list[row].chr[j] = '\0'; // insert hdf5 end of string
                            }
                            
                            // HyperSlab position and length
                            hsize_t start[1];
                            start[0] = (i*MAXSTRBLOCK);
                            hsize_t count[] = {ilength};
                            
                            H5::DataSpace memspace(RANK1, count, NULL);
                            
                            // Get position and write data in dataset
                            dataspace.selectHyperslab(H5S_SELECT_SET, count, start); 
                            dataset->write(names_list, mtype, memspace, dataspace);
                            
                            // Release resources
                            delete[] names_list;
                            memspace.close();
                        }
                        
                    } else {
                        
                        int datarows = Rcpp::as<Rcpp::StringVector>(DatasetValues).size();
                        
                        // Convert Dataframe to range list
                        name *names_list = new name[datarows];
                        
                        for(int i=0; i< datarows; i++ )
                        {
                            //..// name n;
                            Rcpp::String wchrom = Rcpp::as<Rcpp::StringVector>(DatasetValues)(i);
                            std::string word = wchrom.get_cstring();
                            
                            //.. 20260423..////  boost::erase_all(word, "\"");
                            word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
                            
                            int j=0;
                            for( j=0; (unsigned)j < word.size() && j < (MAXSTRING-1); j++ ) {
                                        names_list[i].chr[j] = word[j];
                            }
                            
                            names_list[i].chr[j] = '\0'; // insert hdf5 end of string
                        }
                        
                        dataset->write(names_list, mtype);
                        delete[] names_list;
                    }
                    
                    // Release resources
                    dataspace.close();
                }
            } 
            catch(H5::FileIException& error) { // catch failure caused by the H5File operations
                throw std::runtime_error("c++ exception write_hdf5_string_vector (File IException)");
            } catch(H5::DataSetIException& error) { // catch failure caused by the DataSet operations
                throw std::runtime_error("c++ exception write_hdf5_string_vector (DataSet IException)");
            } catch(H5::GroupIException& error) { // catch failure caused by the Group operations
                throw std::runtime_error("c++ exception write_hdf5_string_vector (Group IException)");
            } catch(H5::DataSpaceIException& error) { // catch failure caused by the DataSpace operations
                throw std::runtime_error("c++ exception write_hdf5_string_vector (DataSpace IException)");
            } catch(std::exception &ex) {
                throw std::runtime_error(std::string("c++ exception write_hdf5_string_vector: ") + ex.what());
            } catch (...) {
                throw std::runtime_error("C++ exception write_hdf5_string_vector (unknown reason)");
            }
            
            dataset->close();
            return void();
        }
    };
}

#endif // BIGDATASTATMETH_HDF5_DIMS_HPP
