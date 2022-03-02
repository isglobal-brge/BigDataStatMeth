#include "include/utils_sparse.h"


// Read data from hdf5 file and Test if loaded block is sparse
bool is_hdf5_readed_block_sparse(H5File* file, DataSet* dataset,
                                 hsize_t offsetx, hsize_t offsety, 
                                 hsize_t countx, hsize_t county)
{
   
   bool bsparse = true;
   
   Eigen::MatrixXd mat = GetCurrentBlock_hdf5_Original(file, dataset, offsetx, offsety, countx, county );
   
   Eigen::SparseMatrix<double> mat_sp = mat.sparseView();
   
   int nnz = mat_sp.nonZeros();
   double numel = mat_sp.rows() * mat_sp.cols();
   double sparsity = (numel - nnz)/numel;
   
   if(sparsity < 0.75)
      bsparse = false;
   
   return(bsparse);
   
}


// Test if eigen matrix is sparse
bool is_sparse( Eigen::MatrixXd mat)
{
   bool bsparse = true;
   
   Eigen::SparseMatrix<double> mat_sp = mat.sparseView();
   
   int nnz = mat_sp.nonZeros();
   double numel = mat_sp.rows() * mat_sp.cols();
   double sparsity = (numel - nnz)/numel;
   
   if(sparsity < 0.75)
      bsparse = false;
   
   return(bsparse);
   
}
