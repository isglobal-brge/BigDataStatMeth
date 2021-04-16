#ifndef parallel_blockmult_Strassen
#define parallel_blockmult_Strassen

   #include <RcppEigen.h>
   #include <iostream>
   #include <time.h>
   #include <iomanip>
   #include <string>
   #include <math.h>
   #include <omp.h>
   using namespace std;
   
   //..// void FillMatrix(double matrix[], int dimension, int maxInt);
   void StrassenMult(double matrix1[], double matrix2[], double matrix3[], int dim);
   void StrassenMultOpenMP(double matrix1[], double matrix2[], double matrix3[], 
                           int dim, int recursions);
   void FillSubmatrices(double matrix[], int dim, double sub1[],
                        double sub2[], double sub3[], double sub4[],
                                                                 int subDim);
   void AddMatrices(double matrix1[], double matrix2[], double matrix3[], int dim);
   void SubtractMatrices(double matrix1[], double matrix2[], double matrix3[], int dim);
   //..// void DisplayMatrix(double matrix[], int dim);
   void FillWithQuads(double quad1[], double quad2[], double quad3[],
                      double quad4[], int subDim, double matrix[], int dim);

#endif