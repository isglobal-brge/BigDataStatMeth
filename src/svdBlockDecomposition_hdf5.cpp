#include "include/svdBlockDecomposition_hdf5.h"


/***
 * Que haurem de fer per fer la descomposició llegint des del fitxer ??
 * 
 * nota : al escriure les dades al fitxer implementar alguna cosa que guardi la mida de la matriu sense haver-la de llegir
 *        només accedint als atributs i prou.
 *        
 * Llegir blocs de tamany x (quin?, fer proves per decidir-ho.... matriu/4 ; matriu/8; ... )
 * 
 * Desprès  implementar : 
 * https://scicomp.stackexchange.com/questions/26562/block-matrix-svd-and-rank-bounds
 * 
 *            A     C
 *      M = 
 *            B     D
 *  
 *      Fem descomposició (SVD) per blocs i tenim : 
 *  
 *            Ua  0   Uc  0       Sa    0     0     0           Va'   0
 *        M =                 x   0     Sb    0     0     x     Vb'   0
 *            0   Ub  0  Ud       0     0     Sc    0           0     Vc'
 *                                0     0     0     Sd          0     Vd'
 *                  X                       S                       Y'
 * 
 *        Calculem QR factors : 
 *            X = Qx * Rx            Y = Qy * Ry
 *    
 *        A partir de la matriu W = Rx * S * Ry'    y obtenint SVD    W = U_w * S * V_w'
 *        
 *            M = Qx * W * Qy' = Qx*Uw * S * Vw'*Qy'  = USV'
 *
 *        Finalment obtenim U i V :
 *        
 *            U = Qx*Uw  ;  S = S   ;   Vt = Vw'*Qy'
 */


// FALTA IMPLEMENTAR-HO --> 1r VULL IMPLEMENTAR EL SVD DIRECTAMENT DES DE FITXER PERÒ
//    UTILITZANT LA FUNCIÓ QUE JA TINC FETA.... ES A DIR, LA IMPLEMENTACIÓ SENSE BLOCS.

svdeig bdSVD_hdf5(Rcpp::RObject X)
{
  
}


// Reads big matrix from hdf5 file in blocks and perform a svd descomposition from
// each block,results are saved in hdf5 datasets under temporal group to be processed
// if necessary
int First_level_SvdBlock_decomposition_hdf5(H5File* file, DataSet* dataset, int k, int q, int nev, bool bcenter, bool bscale, 
                                            int irows, int icols, Rcpp::Nullable<int> threads = R_NilValue)
{
  
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  Eigen::MatrixXd nX;
  svdeig retsvd;
  int  M, p, n;
  int maxsizetoread;
  bool transp = false;
  std::string strGroupName  = "tmpgroup";
  Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2,icols);
  

  if( exists_HDF5_element_ptr(file,strGroupName))
    remove_HDF5_element_ptr(file,strGroupName);

  int ret = create_HDF5_group_ptr(file, strGroupName);
  
  try{
    
    if(irows > icols) {
      // Work with transposed matrix
      n = icols;
      p = irows;
      int i2 = 2;
      transp = true;
      
      // Get data to normalize matrix
      get_HDF5_mean_sd_by_column_ptr( file, dataset, datanormal);

    } else {
      n = irows;
      p = icols;
    }
    
    M = pow(k, q);
    if(M>p)
      throw std::runtime_error("k^q must not be greater than the number of columns of the matrix");
    
    double block_size = std::ceil((double)p/(double)M); // prèviament ja em controlat si p son files o columnes tenint en compte les 


    // Get data from M blocks in initial matrix
    for( int i = 0; i< M ; i++) 
    {
      
      // 1.- Get SVD from all blocks
      //    a) Get all blocks from initial matrix 
      
      IntegerVector offset = getInitialPosition( transp, (unsigned long long)(i*block_size) ); // Posició inicial lectura
      IntegerVector count; // Tamany de block
      maxsizetoread = block_size;

      // Get max block size to read - for blocks smaller than default block size 
      if(transp == true)
      {
        if( ((i+1)*block_size) > irows)
          maxsizetoread = irows - (i*block_size);
      } else {
        if( ((i+1)*block_size) > icols)
          maxsizetoread = icols - (i*block_size);
      }

      count = getSizetoRead(transp, (unsigned long long)(maxsizetoread), icols, irows );
      Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, dataset, offset[0], offset[1], count[0], count[1]);
      
      if(transp==false)
        X.transposeInPlace();

      
      // Normalize data
      if (bcenter==true || bscale==true)
        X = RcppNormalize_Data_hdf5(X, bcenter, bscale, transp, datanormal);


      //    b) SVD for each block
      retsvd = RcppbdSVD_lapack(X, false, false);

      
      //    c)  U*d
      // Create diagonal matrix from svd decomposition d
      int isize = (retsvd.d).size();
      Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
      d.diagonal() = retsvd.d;
      

      std::string strDatasetName = strGroupName + "/A" + std::to_string(i/(M/k));
      
      Eigen::MatrixXd restmp = Bblock_matrix_mul_parallel(retsvd.u, d, 128, threads);
      

      //    d) Escriure els resultats provisionals en alguna part tmp del fitxer hdf5
      offset[0] = 0; offset[1] = 0;
      if(transp == 1 )
      {
        count[0] = restmp.rows();
        count[1] = restmp.cols();
      }else {
        count[0] = restmp.cols();
        count[1] = restmp.rows();
      }

      if(i%(M/k) == 0) {
        // If dataset exists --> remove dataset
        if( exists_HDF5_element_ptr(file,strDatasetName))
          remove_HDF5_element_ptr(file,strDatasetName);
        
        
        // Create unlimited dataset in hdf5 file
        if( transp == true )
          create_HDF5_unlimited_dataset_ptr(file, strDatasetName, count[0], count[1], "numeric");
        else
          create_HDF5_unlimited_dataset_ptr(file, strDatasetName, restmp.cols(), restmp.rows(), "numeric");
        
      } else {
        // Get write position
        if(maxsizetoread == block_size)
        {
          if(transp == 1)
            offset[1] = (i%(M/k))*block_size;
          else
            offset[0] = (i%(M/k))*block_size;
        }
        else{
          if(transp==1)
            offset[1] = ( (i%(M/k))-1 )*block_size + maxsizetoread ;
          else
            offset[0] = ( (i%(M/k))-1 )*block_size + maxsizetoread ;
        }
      }
      
      DataSet* unlimDataset = new DataSet(file->openDataSet(strDatasetName));
      
      // Extend dataset before put data
      if((i%(M/k)) != 0)
      {
        if(transp == true)
          extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
        else
          extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[0]);
        
      }
        
      
      if(transp == true)
        write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp)  );  
      else
        write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp.transpose())  );  
      unlimDataset->close();
      
    }
    
    
  }catch( FileIException error ){ // catch failure caused by the H5File operations
    error.printErrorStack();
    file->close();
    return -1;
  } catch( DataSetIException error ) { // catch failure caused by the DataSet operations
    error.printErrorStack();
    file->close();
    return -1;
  } catch( DataSpaceIException error ) { // catch failure caused by the DataSpace operations
    error.printErrorStack();
    file->close();
    return -1;
  } catch( DataTypeIException error ) { // catch failure caused by the DataSpace operations
    error.printErrorStack();
    file->close();
    return -1;
  }catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    file->close();
    return -1;
  }
  
  return 0;
  
  
}



// Reads small datasets from hdf5 and perform a svd descomposition from each block,
// results are saved in hdf5 datasets under temporal group to be processed if necessary
int Next_level_SvdBlock_decomposition_hdf5(H5File* file, std::string strGroupName, int k, int q, 
                                           bool bcenter, bool bscale, Rcpp::Nullable<int> threads = R_NilValue)
{
  
  IntegerVector stride = IntegerVector::create(1, 1);
  IntegerVector block = IntegerVector::create(1, 1);
  IntegerVector count = IntegerVector::create(0, 0);
  IntegerVector offset = IntegerVector::create(0, 0);
  Eigen::MatrixXd nX;
  svdeig retsvd;
  int M;
  // int nconv, M, p, n;
  // int maxsizetoread;
  
  CharacterVector strvmatnames = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
  
  try{
    
    // Get dataset names
    StringVector joindata =  get_dataset_names_from_group(file, strGroupName, (std::string)strvmatnames[q-1]);
  
    M = joindata.size();
    
    // Get data from M blocks in initial matrix
    for( int i = 0; i< M ; i++) 
    {
      
      //    a) Get dataset
      DataSet currentdataset = file->openDataSet(strGroupName + "/" + joindata[i]);
      IntegerVector dims_out = get_HDF5_dataset_size(currentdataset);
      
      Eigen::MatrixXd X = GetCurrentBlock_hdf5( file, &currentdataset, 0, 0, dims_out[0], dims_out[1]);
      
      //    b) Get dataset svd
      retsvd = RcppbdSVD_lapack(X, false, false);
      
      //    c) U*d
      int isize = (retsvd.d).size();
      Eigen::MatrixXd d = Eigen::MatrixXd::Zero(isize, isize);
      d.diagonal() = retsvd.d;
      
      std::string strDatasetName = strGroupName + "/" + strvmatnames[q] + std::to_string(i/(M/k));
      
      Eigen::MatrixXd restmp = Bblock_matrix_mul_parallel(retsvd.u, d, 128, threads);
      
      //    d) Write results to dataset
      int cbefore = restmp.cols();
      offset[0] = 0; offset[1] = 0;
      count[0] = restmp.rows();
      count[1] = restmp.cols();
      
      if(i%(M/k) == 0) {
        // If dataset exists --> remove dataset
        if( exists_HDF5_element_ptr(file,strDatasetName))
          remove_HDF5_element_ptr(file,strDatasetName);
        // Create unlimited dataset in hdf5 file
        create_HDF5_unlimited_dataset_ptr(file, strDatasetName, restmp.rows(), restmp.cols(), "numeric");
        
      } else {
        // Get initial write position
        offset[1] = offset[1] + cbefore;
      }
      
      DataSet* unlimDataset = new DataSet(file->openDataSet(strDatasetName));
      // Extend dataset before put data
      if((i%(M/k)) != 0)
        extend_HDF5_matrix_subset_ptr(file, unlimDataset, 0, count[1]);
      
      write_HDF5_matrix_subset_v2(file, unlimDataset, offset, count, stride, block, Rcpp::wrap(restmp)  );  
      unlimDataset->close();
      
    }
    
    remove_HDF5_multiple_elements_ptr(file, strGroupName, joindata);
    
  }catch(std::exception &ex) {
    Rcpp::Rcout<< ex.what();
    return -1;
  }
  
  return 0;
  
  
}

