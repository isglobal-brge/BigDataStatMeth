## ########################## ##
##    Test BDSM functions     ##
## ########################## ##

# library definition

require(DelayedArray)
require(BigDataStatMeth)
require(microbenchmark)

# Working directory
#..# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/BigDataStatMeth") # Project directory
setwd( "/Users/mailos/BitDataStatMeth_test/BasicMatVect" ) # Working directory

# Define object sizes

n <- 1000
m <- 1500

p <- 1500
q <- 1200

# R Object
 
matA <- matrix(runif(n*m), nrow = n, ncol = m)
matB <- matrix(runif(p*q), nrow = p, ncol = q)

# DelayedArray Object

matDA <- DelayedArray(matA)
matDB <- DelayedArray(matB)
 
 
# HDF5 File

Create_HDF5_matrix_file("BasicMatVect.hdf5", matA, "INPUT", "matA")
Create_HDF5_matrix(matB, "BasicMatVect.hdf5", "INPUT", "matB")


## -------------------------------------------
##  Basic functions with vectors and matrices 
## -------------------------------------------



# 1.- Matrix product :
 
#  a) R Objects
  
#   a.1) Serial

    res.s <- BigDataStatMeth::blockmult(matA, matB, paral = FALSE)

#   a.2) Parallel

    res.p <- BigDataStatMeth::blockmult(matA, matB, paral = TRUE)

#   a.3) R function

    res.r <- matA%*%matB

#   a.4) Time and equallity

    all.equal(res.s, res.p)
    all.equal(res.s, res.r)
    
     restime <- microbenchmark(
       BigDataStatMeth::blockmult(matA, matB, paral = FALSE),
       BigDataStatMeth::blockmult(matA, matB, paral = TRUE),
       matA%*%matB,
       times = 10)
    
     restime 
     

     # Best time : parallel implementation
    
        
#  b) DelayedArray Objects
   
#   b.1) Serial
     
     res.s <- BigDataStatMeth::blockmult(matDA, matDB, paral = FALSE)
     
#   b.2) Parallel
     
     res.p <- BigDataStatMeth::blockmult(matDA, matDB, paral = TRUE)
     
#   b.3) R function
     
     res.r <- matA%*%matB
     
#   b.4) All equal and Time
     
     all.equal(res.s, res.p)
     all.equal(res.s, res.r)
     
     restime <- microbenchmark(
       BigDataStatMeth::blockmult(matDA, matDB, paral = FALSE),
       BigDataStatMeth::blockmult(matDA, matDB, paral = TRUE),
       matDA%*%matDB,
       times = 10)
     
     restime 
     
     # Best time : parallel implementation

  
#  c) HDF5 files
  
#   c.1) hdf5 datasets 
     
     res <- blockmult_hdf5("BasicMatVect.hdf5", "INPUT","matA","matB")
     r <- rhdf5::H5Fopen(res$file)
     res.hdf5 <- r$OUTPUT$matA_x_matB
     rhdf5::h5closeAll()
       
#   c.2) R function
     
     res.r <- matA%*%matB
     
#   c.3) Time and equallity
     
     all.equal(res.hdf5, res.r)
     
     restime <- microbenchmark(
       blockmult_hdf5("BasicMatVect.hdf5", "INPUT","matA","matB"),
       matA%*%matB,
       times = 1)
     
     restime 
     
     # Best time : HDF5 with big mattrix --> Relation : more big, moredifference, better hdf5 performance in relation to R Objects.
#   