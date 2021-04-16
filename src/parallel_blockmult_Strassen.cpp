#include "include/parallel_blockmult_Strassen.h"

using namespace Rcpp;


/*******************************************************************************************
 This function is an adaptation of the function MPI56.cpp available at 
   https://github.com/tjnowak/MatrixMultiplicationParallelization/blob/master/MPI56.cpp 
 to Rcpp
 
 This algorithm is a parallel implementation of Strassen's matrix multiplication algorithm
 using MPI. This program MUST be called with 56 nodes in the command line.
 The implementation will recursively use up to 56 nodes for matrix multiplication if 
 needed. If more than one recursive call to the matrix multiplication algorithm is required,
 a single node will be used for subsequent calls.
 To compile/run in linux:
 mpicxx -std=c++11 -g -Wall -o mpi56 MPI56.cpp
 mpirun -np 56 ./mpi56 [dimension] [max integer]
 *******************************************************************************************/



// Multiply two matrices using Strassen's algorithm
void StrassenMult(double matrix1[], double matrix2[], double matrix3[], int dim)
{
   
   // Check for matrices with 1 element
   if (dim == 1)
   {
      matrix3[0] = matrix1[0] * matrix2[0];          // only int multipication needed
      return;
   }
   
   if (dim == 2)
   {
      double p1, p2, p3, p4, p5, p6, p7;         // hold results of Strassen's 7 equations
      
      // Find 7 equation results for 2 x 2 matrices
      p1 = matrix1[0] * (matrix2[1] - matrix2[3]);   // a(f-h)
      p2 = (matrix1[0] + matrix1[1]) * matrix2[3];   // (a+b)h
      p3 = (matrix1[2] + matrix1[3]) * matrix2[0];   // (c+d)e
      p4 = matrix1[3] * (matrix2[2] - matrix2[0]);   // d(g-e)
      p5 = (matrix1[0] + matrix1[3]) *
         (matrix2[0] + matrix2[3]);                // (a+d)(e+h)
      p6 = (matrix1[1] - matrix1[3]) *
         (matrix2[2] + matrix2[3]);                // (b-d)(g+h)
      p7 = (matrix1[0] - matrix1[2]) *
         (matrix2[0] + matrix2[1]);                // (a-c)(e+f)
      
      // Fill result matrix3 based on p1-p7
      matrix3[0] = p5 + p4 - p2 + p6;
      matrix3[1] = p1 + p2;
      matrix3[2] = p3 + p4;
      matrix3[3] = p1 + p5 - p3 - p7;
   }
   else
   {
      int subDim = dim / 2;                        // dim for each quadrant of the matrices
      int numElements = subDim * subDim;           // number of elements in sub-matrices
      
      // Initialize sub-matrices
      double* a = new double[numElements];         // top, left quadrant of matrix1
      double* b = new double[numElements];         // top, right quadrant of matrix1
      double* c = new double[numElements];         // bottom, left quadrant of matrix1
      double* d = new double[numElements];         // bottom, right quadrant of matrix1
      double* e = new double[numElements];         // top, left quadrant of matrix2
      double* f = new double[numElements];         // top, right quadrant of matrix2
      double* g = new double[numElements];         // bottom, left quadrant of matrix2
      double* h = new double[numElements];         // bottom, right quadrant of matrix2
      double* result1 = new double[numElements];   // the result of a sub-matrix operation
      double* result2 = new double[numElements];   // the result of a sub-matrix operation
      double* m1 = new double[numElements];        // matrix with result of Strassen's eq. 1
      double* m2 = new double[numElements];        // matrix with result of Strassen's eq. 2
      double* m3 = new double[numElements];        // matrix with result of Strassen's eq. 3
      double* m4 = new double[numElements];        // matrix with result of Strassen's eq. 4
      double* m5 = new double[numElements];        // matrix with result of Strassen's eq. 5
      double* m6 = new double[numElements];        // matrix with result of Strassen's eq. 6
      double* m7 = new double[numElements];        // matrix with result of Strassen's eq. 7
      double* quad1 = new double[numElements];     // top, left quadrant of matrix3
      double* quad2 = new double[numElements];     // top, right quadrant of matrix3
      double* quad3 = new double[numElements];     // bottom, left quadrant of matrix3
      double* quad4 = new double[numElements];     // bottom, right quadrant of matrix3
      
      // Fill sub-matrices a-h from matrix1 and matrix2
      FillSubmatrices(matrix1, dim, a, b, c, d, subDim);
      FillSubmatrices(matrix2, dim, e, f, g, h, subDim);
      
      // Find matrices m1-m7 with results for equations 1-7
      SubtractMatrices(f, h, result1, subDim);            // f-h
      StrassenMult(a, result1, m1, subDim);               // a(f-h)
      
      AddMatrices(a, b, result1, subDim);                 // a+b
      StrassenMult(result1, h, m2, subDim);               // (a+b)h
      
      AddMatrices(c, d, result1, subDim);                  // c+d
      StrassenMult(result1, e, m3, subDim);                // (c+d)e
      
      SubtractMatrices(g, e, result1, subDim);             // g-e
      StrassenMult(d, result1, m4, subDim);                // d(g-e)
      
      AddMatrices(a, d, result1, subDim);                  // a+d
      AddMatrices(e, h, result2, subDim);                  // e+h
      StrassenMult(result1, result2, m5, subDim);          // (a+d)(e+h)
      
      SubtractMatrices(b, d, result1, subDim);             // b-d
      AddMatrices(g, h, result2, subDim);                  // g+h
      StrassenMult(result1, result2, m6, subDim);          // (b-d)(g+h)
      
      SubtractMatrices(a, c, result1, subDim);             // a-c
      AddMatrices(e, f, result2, subDim);                  // e+f
      StrassenMult(result1, result2, m7, subDim);          // (a-c)(e+f)
      
      // Determine quadrants of matrix3 based on m1-m7
      AddMatrices(m5, m4, result1, subDim);                // m5+m4
      SubtractMatrices(result1, m2, result2, subDim);      // m5+m4-m2
      AddMatrices(result2, m6, quad1, subDim);             // m5+m4-m2+m6
      
      AddMatrices(m1, m2, quad2, subDim);                  // m1+m2
      
      AddMatrices(m3, m4, quad3, subDim);                  // m3+m4
      
      AddMatrices(m1, m5, result1, subDim);                // m1+m5
      SubtractMatrices(result1, m3, result2, subDim);      // m1+m5-m3
      SubtractMatrices(result2, m7, quad4, subDim);        // m1+m5-m3-m7
      
      // Fill matrix3 from quadrants
      FillWithQuads(quad1, quad2, quad3, quad4, subDim, matrix3, dim);
      
      // Deallocate sub-matrices
      delete [] a;
      delete [] b;
      delete [] c;
      delete [] d;
      delete [] e;
      delete [] f;
      delete [] g;
      delete [] h;
      delete [] result1;
      delete [] result2;
      delete [] m1;
      delete [] m2;
      delete [] m3;
      delete [] m4;
      delete [] m5;
      delete [] m6;
      delete [] m7;
      delete [] quad1;
      delete [] quad2;
      delete [] quad3;
      delete [] quad4;
   }
}

// Multiply two matrices using Strassen's algorithm and OpenMP
void StrassenMultOpenMP(double matrix1[], double matrix2[], double matrix3[],
                        int dim, int recursions)
{
   
   // Check for matrices with 1 element
   if (dim == 1)
   {
      matrix3[0] = matrix1[0] * matrix2[0];      // only int multipication needed
   }
   
   if (dim == 2)
   {
      double p1, p2, p3, p4, p5, p6, p7;         // hold results of Strassen's 7 equations
      
      // Find 7 equation results for 2 x 2 matrices
      p1 =  matrix1[0] * (matrix2[1]  - matrix2[3]);  // a(f-h)
      p2 = (matrix1[0] +  matrix1[1]) * matrix2[3];   // (a+b)h
      p3 = (matrix1[2] +  matrix1[3]) * matrix2[0];   // (c+d)e
      p4 =  matrix1[3] * (matrix2[2]  - matrix2[0]);  // d(g-e)
      p5 = (matrix1[0] +  matrix1[3]) *
         (matrix2[0] + matrix2[3]);                 // (a+d)(e+h)
      p6 = (matrix1[1] -  matrix1[3]) *
         (matrix2[2] + matrix2[3]);                 // (b-d)(g+h)
      p7 = (matrix1[0] -  matrix1[2]) *
         (matrix2[0] + matrix2[1]);                 // (a-c)(e+f)
      
      // Fill result matrix3 based on p1-p7
      matrix3[0] = p5 + p4 - p2 + p6;
      matrix3[1] = p1 + p2;
      matrix3[2] = p3 + p4;
      matrix3[3] = p1 + p5 - p3 - p7;
   }
   if (dim > 2)
   {
      int subDim = dim / 2;                        // dim for each quadrant of the matrices
      int numElements = subDim * subDim;           // number of elements in sub-matrices
      
      // Initialize sub-matrices
      double* a = new double[numElements];         // top, left quadrant of matrix1
      double* b = new double[numElements];         // top, right quadrant of matrix1
      double* c = new double[numElements];         // bottom, left quadrant of matrix1
      double* d = new double[numElements];         // bottom, right quadrant of matrix1
      double* e = new double[numElements];         // top, left quadrant of matrix2
      double* f = new double[numElements];         // top, right quadrant of matrix2
      double* g = new double[numElements];         // bottom, left quadrant of matrix2
      double* h = new double[numElements];         // bottom, right quadrant of matrix2
      double* m1 = new double[numElements];        // matrix with result of Strassen's eq. 1
      double* m2 = new double[numElements];        // matrix with result of Strassen's eq. 2
      double* m3 = new double[numElements];        // matrix with result of Strassen's eq. 3
      double* m4 = new double[numElements];        // matrix with result of Strassen's eq. 4
      double* m5 = new double[numElements];        // matrix with result of Strassen's eq. 5
      double* m6 = new double[numElements];        // matrix with result of Strassen's eq. 6
      double* m7 = new double[numElements];        // matrix with result of Strassen's eq. 7
      double* temp1 = new double[numElements];     // matrix with intermediate results
      double* temp2 = new double[numElements];     // matrix with intermediate results
      double* quad1 = new double[numElements];     // top, left quadrant of matrix3
      double* quad2 = new double[numElements];     // top, right quadrant of matrix3
      double* quad3 = new double[numElements];     // bottom, left quadrant of matrix3
      double* quad4 = new double[numElements];     // bottom, right quadrant of matrix3
      
      // Fill sub-matrices a-d from matrix1
      FillSubmatrices(matrix1, dim, a, b, c, d, subDim);
      
      // Fill sub-matrices e-f from matrix2
      FillSubmatrices(matrix2, dim, e, f, g, h, subDim);
      
      // Find matrices m1-m7 with results for equations 1-7
#       pragma omp parallel num_threads(7)   // fork 7 threads
{
   if (recursions > 0)
   {
#               pragma omp sections          // allow only one thread to execute a section
{
#                   pragma omp section
{
   double* result1 = new double[numElements];  
   SubtractMatrices(f, h, result1, subDim);            // f-h
   StrassenMultOpenMP(a, result1, m1, subDim,          // a(f-h)
                      (recursions - 1));
   delete [] result1;
}
   
#                   pragma omp section
{
   double* result1 = new double[numElements];   
   AddMatrices(a, b, result1, subDim);                 // a+b
   StrassenMultOpenMP(result1, h, m2, subDim,          // (a+b)h
                      (recursions - 1));    
   delete [] result1;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];  
   AddMatrices(c, d, result1, subDim);                  // c+d
   StrassenMultOpenMP(result1, e, m3, subDim,           // (c+d)e
                      (recursions - 1)); 
   delete [] result1;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];  
   SubtractMatrices(g, e, result1, subDim);             // g-e
   StrassenMultOpenMP(d, result1, m4, subDim,           // d(g-e)
                      (recursions - 1));       
   delete [] result1;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];
   double* result2 = new double[numElements];  
   AddMatrices(a, d, result1, subDim);                  // a+d
   AddMatrices(e, h, result2, subDim);                  // e+h
   StrassenMultOpenMP(result1, result2, m5, subDim,     // (a+d)(e+h)
                      (recursions - 1));    
   delete [] result1;
   delete [] result2;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];  
   double* result2 = new double[numElements];  
   SubtractMatrices(b, d, result1, subDim);             // b-d
   AddMatrices(g, h, result2, subDim);                  // g+h
   StrassenMultOpenMP(result1, result2, m6, subDim,     // (b-d)(g+h)
                      (recursions - 1));    
   delete [] result1;
   delete [] result2;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];   
   double* result2 = new double[numElements];  
   SubtractMatrices(a, c, result1, subDim);             // a-c
   AddMatrices(e, f, result2, subDim);                  // e+f
   StrassenMultOpenMP(result1, result2, m7, subDim,     // (a-c)(e+f)
                      (recursions - 1));       
   delete [] result1;
   delete [] result2;
}
}
   }
   else   // if (recursions == 0)
   {
#               pragma omp sections
{
#                   pragma omp section
{
   double* result1 = new double[numElements];   
   SubtractMatrices(f, h, result1, subDim);            // f-h
   StrassenMult(a, result1, m1, subDim);               // a(f-h)
   delete [] result1;
}
   
#                   pragma omp section
{
   double* result1 = new double[numElements];  
   AddMatrices(a, b, result1, subDim);                 // a+b
   StrassenMult(result1, h, m2, subDim);               // (a+b)h
   delete [] result1;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];   
   AddMatrices(c, d, result1, subDim);                  // c+d
   StrassenMult(result1, e, m3, subDim);                // (c+d)e
   delete [] result1;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];   
   SubtractMatrices(g, e, result1, subDim);             // g-e
   StrassenMult(d, result1, m4, subDim);                // d(g-e)
   delete [] result1;
}

#                   pragma omp section
{
   double* result1 = new double[numElements]; 
   double* result2 = new double[numElements];   
   AddMatrices(a, d, result1, subDim);                  // a+d
   AddMatrices(e, h, result2, subDim);                  // e+h
   StrassenMult(result1, result2, m5, subDim);          // (a+d)(e+h)
   delete [] result1;
   delete [] result2;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];  
   double* result2 = new double[numElements];   
   SubtractMatrices(b, d, result1, subDim);             // b-d
   AddMatrices(g, h, result2, subDim);                  // g+h
   StrassenMult(result1, result2, m6, subDim);          // (b-d)(g+h)
   delete [] result1;
   delete [] result2;
}

#                   pragma omp section
{
   double* result1 = new double[numElements];
   double* result2 = new double[numElements];  
   SubtractMatrices(a, c, result1, subDim);             // a-c
   AddMatrices(e, f, result2, subDim);                  // e+f
   StrassenMult(result1, result2, m7, subDim);          // (a-c)(e+f)
   delete [] result1;
   delete [] result2;
}
}
   }
}

// Determine quadrants of matrix3 based on m1-m7
AddMatrices(m5, m4, temp1, subDim);                // m5+m4
SubtractMatrices(temp1, m2, temp2, subDim);        // m5+m4-m2
AddMatrices(temp2, m6, quad1, subDim);             // m5+m4-m2+m6

AddMatrices(m1, m2, quad2, subDim);                // m1+m2

AddMatrices(m3, m4, quad3, subDim);                // m3+m4

AddMatrices(m1, m5, temp1, subDim);                // m1+m5
SubtractMatrices(temp1, m3, temp2, subDim);        // m1+m5-m3
SubtractMatrices(temp2, m7, quad4, subDim);        // m1+m5-m3-m7

// Fill matrix3 from quadrants
FillWithQuads(quad1, quad2, quad3, quad4, subDim, matrix3, dim);

// Deallocate sub-matrices
delete [] a;
delete [] b;
delete [] c;
delete [] d;
delete [] e;
delete [] f;
delete [] g;
delete [] h;
delete [] m1;
delete [] m2;
delete [] m3;
delete [] m4;
delete [] m5;
delete [] m6;
delete [] m7;
delete [] temp1;
delete [] temp2;
delete [] quad1;
delete [] quad2;
delete [] quad3;
delete [] quad4;
   }
}

// Divide a matrix up by quadrant into 4 sub-matrices.  They are
// sub1, sub2, sub3, and sub4 starting in the top, left quadrant
// of the matrix, going left to right, and top to bottom.
void FillSubmatrices(double matrix[], int dim, double sub1[],
                     double sub2[], double sub3[], double sub4[],
                                                              int subDim)
{
   int index1 = 0;                                   // index of a sub-matrix 1 element
   int index2 = 0;                                   // index of a sub-matrix 2 element
   int index3 = 0;                                   // index of a sub-matrix 3 element
   int index4 = 0;                                   // index of a sub-matrix 4 element
   int matrixElement = 0;                            // index of a matrix element
   
   for (int row=0; row<subDim; row++)
      for (int col=0; col<dim; col++)
      {
         if (col < subDim)
         {
            // Set sub1 element from matrix quadrant 1
            sub1[index1] = matrix[matrixElement];
            index1++;
         }
         else
         {
            // Set sub2 element from matrix quadrant 2
            sub2[index2] = matrix[matrixElement];
            index2++;
         }
         matrixElement++;
      }
      
      for (int row=subDim; row<dim; row++)
         for (int col=0; col<dim; col++)
         {
            if (col < subDim)
            {
               // Set sub3 element from matrix quadrant 3
               sub3[index3] = matrix[matrixElement];
               index3++;
            }
            else
            {
               // Set sub4 element from matrix quadrant 4
               sub4[index4] = matrix[matrixElement];
               index4++;
            }
            matrixElement++;
         }
}

// Add two matrices together into a result matrix
void AddMatrices(double matrix1[], double matrix2[], double matrix3[], int dim)
{
   for (int i=0; i<(dim*dim); i++)
      matrix3[i] = matrix1[i] + matrix2[i];
}

// Subtract two matrices.  Subtract matrix2 from matrix1 and put the difference
// into a result matrix.
void SubtractMatrices(double matrix1[], double matrix2[], double matrix3[], int dim)
{
   for (int i=0; i<(dim*dim); i++)
      matrix3[i] = matrix1[i] - matrix2[i];
}

// // Display a matrix
// void DisplayMatrix(double matrix[], int dim)
// {
//    for (int i=0; i<(dim*dim); i++)
//    {
//       if (i % dim == 0)                             // return at end of matrix line
//          cout << endl;
//       cout << " " << matrix[i];
//    }
//    cout << endl << endl;
// }

// Fill a matrix with elements from four separate quadrant sub-matrices
void FillWithQuads(double quad1[], double quad2[], double quad3[],
                   double quad4[], int subDim, double matrix[], int dim)
{
   int index1 = 0;                                   // index of a quad1 element
   int index2 = 0;                                   // index of a quad2 element
   int index3 = 0;                                   // index of a quad3 element
   int index4 = 0;                                   // index of a quad4 element
   int matrixElement = 0;                            // index of a matrix element
   
   for (int row=0; row<subDim; row++)
      for (int col=0; col<dim; col++)
      {
         if (col < subDim)
         {
            // Set matrix element from quad1
            matrix[matrixElement] = quad1[index1];
            index1++;
         }
         else
         {
            // Set matrix element from quad2
            matrix[matrixElement] = quad2[index2];
            index2++;
         }
         matrixElement++;
      }
      
      for (int row=subDim; row<dim; row++)
         for (int col=0; col<dim; col++)
         {
            if (col < subDim)
            {
               // Set matrix element from quad3
               matrix[matrixElement] = quad3[index3];
               index3++;
            }
            else
            {
               // Set matrix element from quad4
               matrix[matrixElement] = quad4[index4];
               index4++;
            }
            matrixElement++;
         }
}





void eigen_matrixXd_to_double_array(const Eigen::MatrixXd& evector, double* destination, int rows, int cols)
{
   typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajMat;
   Eigen::Map<RowMajMat>(destination, rows, cols) = evector;
}








//' Block matrix multiplication with Delayed Array Object
//' 
//' This function performs a block matrix-matrix multiplication with numeric matrix or Delayed Arrays
//' 
//' @param a a double matrix.
//' @param b a double matrix.
//' @param block_size (optional, defalut = 128) block size to make matrix multiplication, if `block_size = 1` no block size is applied (size 1 = 1 element per block)
//' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel computation else performs seria computation
//' @param threads (optional) only if bparal = true, number of concurrent threads in parallelization if threads is null then threads =  maximum number of threads available
//' @param bigmatrix (optiona, default = 5000) maximum number of rows or columns to consider as big matrix and work with
//' hdf5 files, by default a matrix with more than 5000 rows or files is considered big matrix and computation is made in disk 
//' @param mixblock_size (optiona, default = 128), only if we are working with big matrix and parallel computation = true. 
//' Block size for mixed computation in big matrix parallel. Size of the block to be used to perform parallelized memory 
//' memory of the block read from the disk being processed.
//' @param outfile (optional) file name to work with hdf5 if we are working with big matrix in disk.
//' @return A List with : 
//' \itemize{
//'   \item{"matrix"}{ Result matrix if execution has been performed in memory}
//'   \item{"filename"}{ HDF5 filename if execution has been performed in disk, HDF5 file contains : 
//'     \itemize{
//'       \item{"INPUT"}{hdf5 group with input matrix A and B}
//'       \item{"OUTPUT"}{hdf5 group with output matrix C}
//'     }with input and output matrix.
//'   }
//' }
//' @examples
//' 
//' library(DelayedArray)
//' 
//' # with numeric matrix
//' m <- 500
//' k <- 1500
//' n <- 400
//' A <- matrix(rnorm(n*k), nrow=n, ncol=k)
//' B <- matrix(rnorm(n*k), nrow=k, ncol=n)
//' 
//' blockmult(A,B,128, TRUE)
//' 
//' # with Delaeyd Array
//' AD <- DelayedArray(A)
//' BD <- DelayedArray(B)
//' 
//' blockmult_Strassen(AD,BD,128, TRUE)
//' @export
// [[Rcpp::export]]
Rcpp::List blockmult_Strassen( Rcpp::RObject a, Rcpp::RObject b, 
                               Rcpp::Nullable<int> block_size = R_NilValue, 
                               Rcpp::Nullable<bool> paral = R_NilValue,
                               Rcpp::Nullable<int> threads = R_NilValue
                               )
{
   
   //..// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajMat;
   
   
   const int MAX_RECURSIONS = 1;                   // max StrassenMultOpenMP recursions
   
   int dim;                                        // row/col dim of input matrices
   //..// int maxInt;                                     // number of possible element values
   double* firstMatrix = NULL;                     // first input matrix
   double* secondMatrix = NULL;                    // second input matrix
   double* resultMatrix = NULL;                    // matrix mult result
   double startTime;                               // start of matrix mult
   double endTime;                                 // end of matrix mult
   
   // Allow nested OpenMP parallelism
   omp_set_nested(1);
   
//    dim = Rcpp::as<Rcpp::NumericMatrix>(a).nrow();
   dim = Rcpp::as<Rcpp::NumericMatrix>(a).nrow();
 

   // Initialize matrices
   // firstMatrix = new double[dim * dim];
   // secondMatrix = new double[dim * dim];
   resultMatrix = new double[dim * dim];

/***   GEN
 Block per generar les matrius --> Aquí no serveix per a res 
 
   // Initialize random number generator
   srand(time(NULL));
   
   // Fill two matrices
   FillMatrix(firstMatrix, dim, maxInt);
   FillMatrix(secondMatrix, dim, maxInt);
 ***/ 
   
   // Initialize
   Eigen::Map<Eigen::MatrixXd> Xf(Rcpp::as< Eigen::Map<Eigen::MatrixXd> >(a)); 
   Xf.transposeInPlace();
   firstMatrix = Xf.data();
   
   
   Eigen::Map<Eigen::MatrixXd> Xs(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(b));
   secondMatrix = Xs.transpose().data();
   
// 
//    
//    Rcpp::Rcout<<"\n First Matrix : "<<firstMatrix <<" \n ";
//    
//    Rcpp::Rcout<<"\n\n First Matrix :  \n"<< *(firstMatrix+0) <<"\t"<<*(firstMatrix+8)<<"\t"<<*(firstMatrix+16)<<"\t"<<*(firstMatrix+24) <<" \n ";
//    Rcpp::Rcout<<"\n\n\t\tAcabem de mostrar firstmatrix elements\n\n";
//    
//    Rcpp::Rcout<<"\n\n Second Matrix :  \n"<< *(secondMatrix+0) <<"\t"<<*(secondMatrix+8)<<"\t"<<*(secondMatrix+16)<<"\t"<<*(secondMatrix+24) <<" \n ";
//    
//    Rcpp::Rcout<<"\n Dim ... : "<<dim <<" \n\n ";
//    
//    
//    //Rcpp::Rcout<<"\n Second Matrix : "<<secondMatrix << " \n "<< Xs;
//       
//    // vector<vector<double> >::pointer firstMatrix_ptr = Rcpp::as<NumericMatrix>(b).dataptr()[0];
//    
//    // firstMatrix= malloc(sizeof(VECTOR_PTR(Rcpp::as<NumericVector>(a))));
//    // 
//    // //firstMatrix = VECTOR_PTR(Rcpp::as<NumericVector>(a));
//    // secondMatrix = Rcpp::as<NumericMatrix>(b).dataptr();
   
   // Multiply the two matrices using Strassen's algorithm
   startTime = omp_get_wtime();
   StrassenMultOpenMP(firstMatrix, secondMatrix, resultMatrix,
                      dim, MAX_RECURSIONS);
   endTime = omp_get_wtime();
   
   // Display results
   printf("\nMultiplication of dimension %d matrices took %f"
             " seconds\n\n", dim, endTime - startTime);




   //..// printf("first matrix:\n");
   //..// DisplayMatrix(firstMatrix, dim);
   //..// printf("second matrix:\n");
   //..// DisplayMatrix(secondMatrix, dim);
   //..// printf("result:\n");
   //..// DisplayMatrix(resultMatrix, dim);
   
   // Deallocate matrices
   // delete [] firstMatrix;
   // delete [] secondMatrix;
   // delete [] resultMatrix;
   
   return 0;   
   
   
   
   
}




/*** R
library(devtools) # Only to reload package
library(Matrix)
library(BigDataStatMeth)


random_matrix <- function(k, seed=NULL, compute.dense=FALSE){
   set.seed(seed)
   ## Sparse singular matrix ##
   D <- diff(Diagonal(k),differences=1)
   S.sparse <- t(D)%*%D
   ## Sparse non-singular matrix ##
   l <- runif(1)
   NS.sparse <- l*S.sparse+(1-l)*Diagonal(k) 
   ## Dense singular/non-singular matrices ##
   if(compute.dense){
      M <- matrix(rnorm(k^2), nrow=k)
      S.dense <- S.sparse%*%M
      NS.dense <- l*S.dense+(1-l)*Diagonal(k)
   }else{
      S.dense <- NULL
      NS.dense <- NULL
   }
   return(list(S.sparse=S.sparse, NS.sparse=NS.sparse, S.dense=S.dense, NS.dense=NS.dense))
}



k <- 4096

M1 <- random_matrix(k, seed=123, compute.dense=TRUE)
M2 <- random_matrix(k, seed=321, compute.dense=TRUE)

reload(pkgload::inst("BigDataStatMeth"))

res <- microbenchmark::microbenchmark( blockmult_Strassen( as.matrix(M1$S.dense), as.matrix(M2$S.dense)),
                                       BigDataStatMeth::blockmult(as.matrix(M1$S.dense), as.matrix(M2$S.dense)),
                                       as.matrix(M1$S.dense)%*% as.matrix(M2$S.dense),
                                       times = 3)

print(res)



library(devtools) # Only to reload package
library(Matrix)


library(BigDataStatMeth)


a <- matrix(c(4.08871, 6.24779, 4.24913, 1.68276, 6.71266, 5.53369, 3.6224, 5.86116,
                3.3737, 6.17535, 6.66233, 4.91872, 5.30112, 3.92453, 6.60796, 3.42434,
                4.39426, 5.4531, 6.05339, 3.25049, 5.40845, 7.83502, 4.8802, 5.78476,
                6.07524, 6.75907, 3.22968, 8.29333, 6.54057, 2.33085, 6.26718, 4.82074,
                6.90316, 6.93738, 3.78272, 4.65016, 6.34638, 6.14092, 4.43826, 5.00616,
                7.36931, 6.77012, 1.19809, 3.78208, 6.11539, 3.32425, 2.53041, 5.9477,
                5.67361, 7.02902, 4.02655, 5.35905, 7.14915, 2.98546, 6.49869, 4.25027,
                7.14703, 4.17957, 5.16678, 4.241, 4.87737, 1.04653, 3.77343, 6.70771), nrow=8, byrow=TRUE )



b <- matrix(c( 1.77284, 8.12854, 6.01416, 4.90193, 8.60554, 6.73623, 8.16294, 5.11955,
               3.13505, 2.11654, 1.32864, 3.11144, 5.27269, 4.51073, 3.3419, 6.07763,
               5.36682, 5.55259, 1.30046, 3.87732, 4.89205, 5.61549, 4.53863, 8.1588,
               5.90918, 6.2177, 2.26226, 7.65156, 3.8853, 5.53137, 0.340948, 4.65508,
               6.20073, 2.4694, 5.60392, 5.124, 4.25624, 6.33544, 1.32355, 4.94445,
               7.5545, 2.97776, 3.64464, 3.37875, 6.21592, 4.94152, 4.43024, 7.30229,
               4.7526, 4.47249, 7.61619, 5.41288, 7.85611, 1.12197, 6.29337, 5.54141,
               8.28032, 5.21989, 3.26278, 5.92957, 6.71789, 4.78047, 5.19883, 5.3405), nrow=8, byrow=TRUE )

#library(devtools) # Only to reload package
# reload(pkgload::inst("BigDataStatMeth"))
blockmult_Strassen( a, b)


*/
