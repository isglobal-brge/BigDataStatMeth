## ----setup, include = FALSE---------------------------------------------------
library(knitr)
library(BiocStyle)

# knitr::opts_chunk$set(collapse = TRUE, comment = "", cache = FALSE, message = FALSE, width = 180, crop = NULL)

## ----cleanup, echo=FALSE, include=FALSE---------------------------------------
hdf5_close_all <- function(...) invisible(NULL)  # safe no-op before package loads
for (f in c('delayed.hdf5', 'robject.hdf5', 'colesterol_file.hdf5',
            'subtraction_example.hdf5', 'addition_example.hdf5')) {
  if (isTRUE(file.exists(f))) unlink(f)
}

## ----install_required, eval=FALSE---------------------------------------------
# # Install BiocManager (if not previously installed)
# install.packages("BiocManager")
# 
# # Install required packages
# BiocManager::install(c("Matrix", "RcppEigen", "RSpectra", "HDF5Array"))

## ----install, eval=FALSE------------------------------------------------------
# # Install devtools and load library (if not previously installed)
# install.packages("devtools")
# library(devtools)
# 
# # Install BigDataStatMeth
# install_github("isglobal-brge/BigDataStatMeth")

## ----load, cache=FALSE--------------------------------------------------------
library(BigDataStatMeth)

## ----hdf5Img, out.width = '100%', fig.align = 'center', fig.cap = "HDF5 hierarchical structure", echo=FALSE----
# knitr::include_graphics("imgs/hdf5_squema.jpg")
# 
## ----hdf5Create---------------------------------------------------------------
set.seed(5234)
n <- 500
m <- 600
A <- matrix(rnorm(n*m,mean=0,sd=1), n,m)

# Create a dataset from R matrix object; returns an HDF5Matrix
A_hdf5 <- hdf5_create_matrix("robject.hdf5", "INPUT/A", data = A)
A_hdf5

## ----ls-----------------------------------------------------------------------
list.files(pattern = "*.hdf5")

## ----hdf5AddDataset-----------------------------------------------------------
set.seed(5234)
n <- 500
m <- 1000
A <- matrix(rnorm(n*m,mean=0,sd=1), n, m)

set.seed(5234)
n <- 1000
m <- 12000
B <- matrix(rnorm(n*m,mean=3,sd=0.5), n, m)

# Path to HDF5 file
example_fn <- "delayed.hdf5"

# Create file with matrix A in group INPUT
A_hdf5 <- hdf5_create_matrix(example_fn, "INPUT/A", data = A)

# Add matrix B to the same file and group
B_hdf5 <- hdf5_create_matrix(example_fn, "INPUT/B", data = B)

## ----hdf5Show-----------------------------------------------------------------
# List datasets in the file
list_datasets(example_fn, group = "/", recursive = TRUE)

## ----hdf5Open, cache=FALSE----------------------------------------------------
# Open an existing dataset from file — returns an HDF5Matrix
A_hdf5 <- hdf5_matrix("robject.hdf5", "INPUT/A")
A_hdf5

## ----hdf5Dataset--------------------------------------------------------------
B_hdf5[1:3, 1:5]

## ----hdf5DatasetClose---------------------------------------------------------
hdf5_close_all()

## ----convert_HDF5, cache=FALSE------------------------------------------------
import_hdf5 <- hdf5_import(
  source   = "/Users/mailos/PhD/TREBALLANT/BigDataStatMeth_20250822/vignettes/colesterol.csv",
  filename = "colesterol_file.hdf5",
  dataset  = "COLESTEROL/COLESTEROLDATA",
  sep      = ",",
  header   = TRUE,
  overwrite = TRUE
)

## ----read_data_col_HDF5-------------------------------------------------------
# Show content
list_datasets(import_hdf5)

# Show first rows and columns
import_hdf5[1:5, 1:6]

# Show colnames
head(colnames(import_hdf5))

## ----blockmult_hdf5_exec------------------------------------------------------
# Reopen datasets (handles were closed above)
A_hdf5 <- hdf5_matrix(example_fn, "INPUT/A")
B_hdf5 <- hdf5_matrix(example_fn, "INPUT/B")

# Perform blockwise matrix multiplication
res <- A_hdf5 %*% B_hdf5

# Show the content of the HDF5 file
list_datasets(res)

## ----blockmult_hdf5_res-------------------------------------------------------
# Extract a subset of the result from HDF5
result_hdf5 <- res[1:3, 1:5]
result_hdf5

# Compute the same multiplication in R
result_r <- (A %*% B)[1:3, 1:5]
result_r

# Compare both results
all.equal(A %*% B, as.matrix(res))

## ----crossprod_sing-----------------------------------------------------------
# Create example matrices
set.seed(123)
A <- matrix(rnorm(1000 * 200), nrow = 1000, ncol = 200)
B <- matrix(rnorm(1000 * 150), nrow = 1000, ncol = 150)
C <- matrix(rnorm(800 * 200),  nrow = 800,  ncol = 200) 

# Save matrices to HDF5 file
hdf5_close_all()
unlink(example_fn)

A_hdf5 <- hdf5_create_matrix(example_fn, "INPUT/A", data = A)
B_hdf5 <- hdf5_create_matrix(example_fn, "INPUT/B", data = B)
C_hdf5 <- hdf5_create_matrix(example_fn, "INPUT/C", data = C)

# Compute t(A) %*% A
res_cross <- crossprod(A_hdf5)

# Show where the result is stored
list_datasets(res_cross)

# Compare with R's crossprod
res_hdf5 <- as.matrix(res_cross)
res_r <- crossprod(A)

all.equal(res_r, res_hdf5)

## ----crossprod_dbl------------------------------------------------------------
# Compute t(A) %*% B
res_cross2 <- crossprod(A_hdf5, B_hdf5)

# Compare with R
res_hdf5 <- as.matrix(res_cross2)
res_r <- crossprod(A, B)

all.equal(res_r, res_hdf5)

## ----tcrossprod_sing----------------------------------------------------------
# Compute A %*% t(A)
res_tcross <- tcrossprod(A_hdf5)

list_datasets(res_tcross)

res_hdf5_t <- as.matrix(res_tcross)
res_r_t <- tcrossprod(A)

all.equal(res_r_t, res_hdf5_t)

## ----tcrossprod_dbl-----------------------------------------------------------
# Compute A %*% t(C)
res_tcross2 <- tcrossprod(A_hdf5, C_hdf5)

res_hdf5_t <- as.matrix(res_tcross2)
res_r_t <- tcrossprod(A, C)

all.equal(res_r_t, res_hdf5_t)

## ----substract_init-----------------------------------------------------------
# Create two matrices of the same dimensions
set.seed(42)
A_sub <- matrix(rnorm(1000 * 300), nrow = 1000, ncol = 300)
B_sub <- matrix(rnorm(1000 * 300), nrow = 1000, ncol = 300)

# Save them to HDF5
fn_sub <- "subtraction_example.hdf5"

A_sub_hdf5 <- hdf5_create_matrix(fn_sub, "INPUT/A_sub", data = A_sub)
B_sub_hdf5 <- hdf5_create_matrix(fn_sub, "INPUT/B_sub", data = B_sub)

# Perform subtraction: A - B
res_sub <- A_sub_hdf5 - B_sub_hdf5

# Compare a subset with R
result_hdf5 <- as.matrix(res_sub)
result_r <- A_sub - B_sub

all.equal(result_r, result_hdf5)

## ----add----------------------------------------------------------------------
# Create two compatible matrices
set.seed(99)
A_add <- matrix(rnorm(800 * 250), nrow = 800, ncol = 250)
B_add <- matrix(rnorm(800 * 250), nrow = 800, ncol = 250)

# Save them to HDF5
fn_add <- "addition_example.hdf5"

A_add_hdf5 <- hdf5_create_matrix(fn_add, "INPUT/A", data = A_add)
B_add_hdf5 <- hdf5_create_matrix(fn_add, "INPUT/B", data = B_add)

# Perform addition: A + B
res_add <- A_add_hdf5 + B_add_hdf5

# Compare result with R
result_hdf5 <- as.matrix(res_add)
result_r <- A_add + B_add

all.equal(result_r, result_hdf5)

## ----BSVDImg, out.width = '100%', fig.align = 'center', fig.cap = "Flowchart for a two-level hierarchical Block SVD algorithm", echo=FALSE----
# knitr::include_graphics("imgs/blocksvd.png")

## ----BlockSVDNorm-------------------------------------------------------------
# Create dataframe data with 'omicdata' matrix in delayed hdf5 file at OMICS group
set.seed(5234)
n <- 100
m <- 15000
omicdata <- matrix(rnorm(n*m, mean=0, sd=1), n, m)

omic_hdf5 <- hdf5_create_matrix(example_fn, "OMICS/data", data = omicdata, overwrite = TRUE)

# SVD directly on HDF5Matrix — d is returned as a numeric vector
svdh5 <- svd(omic_hdf5)

# Results: d is already a numeric vector, no need to read from file
svdh5$d[1:7]

svd_r <- svd(scale(omicdata))
svd_r$d[1:7]

## ----BlockSVDNotNorm----------------------------------------------------------
# SVD without normalization
svdh5 <- svd(omic_hdf5, center = FALSE, scale = FALSE, overwrite = TRUE)

## ----BlockSVDNotNormResults---------------------------------------------------
# SVD (d) from file - data not normalized
svdh5$d[1:7]

# with R implementation from data in memory
svd(omicdata)$d[1:7]

## ----BlockSVDk4---------------------------------------------------------------
# Block decomposition with 1 level and 4 local SVDs at each level using 
# two threads (as maximum)
svdh5 <- svd(omic_hdf5, k = 4, q = 1, threads = 2, overwrite = TRUE)

# SVD (d) — data centered and scaled by default
svdh5$d[1:7]

# with R implementation from data in memory
svd(scale(omicdata))$d[1:7]

## ----cholDesc-----------------------------------------------------------------
N <- 100
set.seed(5234)
Y <- matrix(rnorm(N*N), N, N)
Ycp <- crossprod(Y)

hdf5_close_all()
unlink(example_fn)

Ycp_hdf5 <- hdf5_create_matrix(example_fn, "chol/data", data = Ycp)

cholh5 <- chol(Ycp_hdf5)
choldesc_hdf5 <- as.matrix(cholh5)
choldesc_hdf5[1:3, 1:5]

choldesc_r <- chol(Ycp)
choldesc_r[1:3, 1:5]

all.equal(choldesc_hdf5, choldesc_r)

## ----cholDesc_50--------------------------------------------------------------
cholh5 <- chol(Ycp_hdf5,
               block_size = 50,
               overwrite  = TRUE)

choldesc_hdf5 <- as.matrix(cholh5)
choldesc_hdf5[1:3, 1:5]

choldesc_r <- chol(Ycp)
choldesc_r[1:3, 1:5]

all.equal(choldesc_hdf5, choldesc_r)

## ----QRdec--------------------------------------------------------------------
Ycp_hdf5 <- hdf5_matrix(example_fn, "chol/data")

QRh5 <- qr(Ycp_hdf5, thin = TRUE)

QR_Q_hdf5 <- as.matrix(QRh5$Q)
QR_R_hdf5 <- as.matrix(QRh5$R)

# Q matrix
QR_Q_hdf5[1:3, 1:5]

# R matrix
QR_R_hdf5[1:3, 1:5]

# Q matrix in R
QR_Q_r <- qr.Q(qr(Ycp))
QR_Q_r[1:3, 1:5]

all.equal(QR_Q_hdf5, QR_Q_r)

## ----QRdec_blocksize----------------------------------------------------------
QRh5 <- qr(Ycp_hdf5,
            block_size = 256,
            overwrite  = TRUE)

QR_Q_hdf5 <- as.matrix(QRh5$Q)
QR_R_hdf5 <- as.matrix(QRh5$R)

# Q matrix
QR_Q_hdf5[1:3, 1:5]

# R matrix
QR_R_hdf5[1:3, 1:5]

# Q matrix in R
QR_Q_r <- qr.Q(qr(Ycp))
QR_Q_r[1:3, 1:5]

all.equal(QR_Q_hdf5, QR_Q_r)

## ----CholInv------------------------------------------------------------------
invCholh5 <- solve(Ycp_hdf5)

invChol_hdf5 <- as.matrix(invCholh5)
invChol_hdf5[1:5, 1:5]

## ----CholInv_full-------------------------------------------------------------
invCholh5 <- solve(Ycp_hdf5, full_matrix = TRUE, overwrite = TRUE)

invChol_hdf5 <- as.matrix(invCholh5)
invChol_hdf5[1:5, 1:5]

# Inverse matrix in R
invChol_r <- solve(Ycp)
invChol_r[1:5, 1:5]

all.equal(invChol_hdf5, invChol_r)

## ----sesinfo------------------------------------------------------------------
sessionInfo()

