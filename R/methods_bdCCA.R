
#' QR by blocks
#' 
#' This function is an application of the BigDataStatMeth functions to generate new methods. This function perform a QR
#' from two matrices stored in hdf5 data file. This function applies matrix partitioning, merge bloks to create a full matrix, apply a function to different blocks...
#' 
#' @export
#' 
#' @param strdataset string, dataset path within the hdf5 data file from which we want to calculate the QR
#' @param file string file name where dataset to normalize is stored
#' @param mblocks number of blocks in which we want to partition the matrix to perform the calculations
#' @param center, boolean, if true, dataset is centered to perform calculus
#' @param bcols boolean if bcols = TRUE matrix it´s splitted by columns if bcols = FALSE, then matrix or dataset is splitted by rows.
#' @param scale, boolean, if true, dataset is centered to perform calculus
#' @param overwrt, boolean, if true, datasets existing inside a file must be overwritten if we are using the same names
#' @return hdf5 data file with CCA results, 
#' @examples
#' 
#'    print ("Example in vignette")
#' 
getQRbyBlocks <- function(strdataset, file, mblocks, center, scale, bcols, overwrt)
{
    
    strgroup <- gsub("/.?$", "", strdataset)
    strdataset <- gsub("^.*/", "", strdataset)
    
    
    # Prepare data - Normalize data (only Center)
    bdNormalize_hdf5(filename = file, 
                     group = strgroup, dataset = strdataset, 
                     bcenter = center, bscale = scale, force = overwrt) 

    # Review m size block must be equal than number of samples
    # if( bdgetDim_hdf5(file, paste0(strgroup, "/",strdataset))[1] / mblocks <  bdgetDim_hdf5(file, paste0(strgroup, "/",strdataset))[2] ) {
    #     mblocks <- floor( bdgetDim_hdf5(file, paste0(strgroup, "/",strdataset))[1] / bdgetDim_hdf5(file, paste0(strgroup, "/",strdataset))[2]) - 1
    #     if(mblocks <= 0) mblocks = 1
    #     message("m set to ", mblocks)
    # }

    # Step 1
    # Split datasets X abd Y by rows and store data to data file
    bdSplit_matrix_hdf5( filename = file, 
                         group = paste0("NORMALIZED/",strgroup), 
                         dataset = strdataset, 
                         outgroup = paste0( "Step1/", strdataset, "rows"), 
                         nblocks = mblocks, bycols = bcols, force = overwrt)
    
    # Step 2
    # Get splitted dataset names
    blocks <- BigDataStatMeth::bdgetDatasetsList_hdf5(file, paste0( "Step1/", strdataset, "rows"))
    bdapply_Function_hdf5( filename = file, group = paste0( "Step1/", strdataset, "rows"), 
                           datasets = blocks, 
                           outgroup = paste0( "Step2/", strdataset, "rows"), 
                           func = "QR", 
                           force = overwrt )
    
    # stop("Reviewing")
    
    # Step 3
    blocks.qr <- bdgetDatasetsList_hdf5(file, paste0( "Step2/", strdataset, "rows"))
    bdBind_hdf5(filename = file, group =  paste0( "Step2/", strdataset, "rows"), 
                datasets = blocks.qr[which(blocks.qr %like% ".R")],
                outgroup = "Step3/merged", outdataset =  paste0( strdataset, "Rt"), 
                func = "bindRows", force = overwrt )
    
    # stop("Reviewing")
    bdapply_Function_hdf5( file, "Step3/merged", paste0( strdataset, "Rt"), "Step3/Final_QR", "QR", force = overwrt )
    
    # Step 4
    bdSplit_matrix_hdf5(file, "Step3/Final_QR", paste0( strdataset, "Rt.Q"), 
                        outgroup = "Step4/splitted", 
                        nblocks = mblocks, 
                        bycols = bcols, force = overwrt )
    
    
    # Step 5
    # Get splitted matrices names
    tmp <- bdgetDatasetsList_hdf5(file, "Step4/splitted")
    Rt.Q.divide <- tmp[which(tmp %like% paste0( strdataset, "Rt.Q"))]
    # multiply previous splitted matrices with Q descomposed matrices from model (X)
    bdapply_Function_hdf5(  filename = file, group = paste0( "Step2/", strdataset, "rows"), 
                            datasets = blocks.qr[which(blocks.qr %like% ".Q")], 
                            outgroup = "Step5", func = "blockmult",
                            b_group = "Step4/splitted", b_datasets = Rt.Q.divide,
                            force = TRUE )
    
    # Step 6
    blocks.Q <- bdgetDatasetsList_hdf5(file, "Step5")
    bdBind_hdf5(filename = file, group = "Step5", datasets = blocks.Q[which(blocks.Q %like% paste0(strdataset,"."))],
                outgroup = "Step6", outdataset = paste0(strdataset,"Q"), 
                func = "bindRows", force = TRUE )
    
    
}


writeCCAComponents_hdf5 <- function(filename, ncolsX, ncolsY)
{
    
    if(!file.exists(filename)){
        message("ERROR - File doesn't exists")
        return()
    }
    
    # Read data from file
    h5f = H5Fopen(filename)
    XQ <- h5f$Step6$XQ[1:ncolsX, 1:ncolsX]
    YQ <- h5f$Step6$YQ[1:ncolsY, 1:ncolsY]
    XR <- h5f$Step3$Final_QR$XRt.R
    YR <- h5f$Step3$Final_QR$YRt.R
    d <- h5f$SVD$CrossProd_XQxYQ$d
    u <- h5f$SVD$CrossProd_XQxYQ$u
    v <- h5f$SVD$CrossProd_XQxYQ$v
    xcenter <- h5f$NORMALIZED$data$X.mean
    ycenter <- h5f$NORMALIZED$data$Y.mean
    x.names <- h5f$data$.X_dimnames$`2`
    y.names <- h5f$data$.Y_dimnames$`2`
    h5closeAll()
    
    # Get qr compact (more or less)
    XR[lower.tri(XR, diag = F)] <- 0
    XQ[upper.tri(XQ, diag = T)] <- 0
    XQR <- XR + XQ
    
    YR[lower.tri(YR, diag = F)] <- 0
    YQ[upper.tri(YQ, diag = T)] <- 0
    YQR <- YR + YQ
    
    xcoef <- bdSolve(XQR, u)
    ycoef <- bdSolve(YQR, v)
    
    rownames(xcoef) <- as.matrix(x.names) 
    rownames(ycoef) <- as.matrix(y.names)
    
    # Store results to hdf5 data file under Results group/folder
    #     cor, xcoef, ycoef, xcenter, ycenter, xscores, yscores
    #     corr.X.xscores, corr.Y.xscores, corr.X.yscores, corr.Y.yscores

    print("Getting and writting xcoef and ycoef")
    bdAdd_hdf5_matrix( xcoef, filename,  "Results", "xcoef", force = TRUE)
    bdAdd_hdf5_matrix( ycoef, filename,  "Results", "ycoef", force = TRUE)
    
    print("Getting and writting cor matrix")
    bdAdd_hdf5_matrix( diag(d), filename,  "Results", "cor", force = TRUE)
    
    print("Getting and writting xcenter and ycenter vectors")
    bdAdd_hdf5_matrix( xcenter, filename,  "Results", "xcenter", force = TRUE)
    bdAdd_hdf5_matrix( ycenter, filename,  "Results", "ycenter", force = TRUE)
    
    print("Getting and writting xscores and yscores matrix")
    bdblockmult_hdf5(filename, group = "data", a = "X", b = "xcoef", groupB = "Results", outgroup = "Results", outdataset = "xscores")
    bdblockmult_hdf5(filename, group = "data", a = "Y", b = "ycoef", groupB = "Results", outgroup = "Results", outdataset = "yscores")
    
}


#' Canonical Correlation Analysis
#' 
#' This function is an application of the BigDataStatMeth functions to generate new methods. This function perform a Canonical Correlation Analysis
#' from two matrices stored in hdf5 data file. This function applies matrix partitioning, merge bloks to create a full matrix, apply a function to different blocks, etc.
#' 
#' @export
#' 
#' @param filename string file name where dataset to normalize is stored.
#' @param X Dataset, path inside the hdf5 data file.
#' @param Y Dataset, path inside the hdf5 data file.
#' @param m Integer, number of blocks in which we want to partition the matrix to perform the calculations.
#' @param bcenter, Boolean, if true, dataset is centered to perform calculus.
#' @param bscale, Boolean, if true, dataset is centered to perform calculus.
#' @param bycols, Boolean by default = true, true indicates that the imputation will be done by columns, otherwise, the imputation will be done by rows.
#' @param overwriteResults, Boolean, if true, datasets existing inside a file must be overwritten if we are using the same names.
#' @param keepInteResults, Boolean, if false, intermediate results will be removed.
#' @return hdf5 data file with CCA results, 
#' @examples
#'    print ("Example in vignette")
#' 
bdCCA_hdf5 <- function(filename, X, Y, m = 10, bcenter = TRUE, bscale = FALSE, bycols = FALSE, overwriteResults = FALSE, keepInteResults = FALSE)
{
    
    # if(file.exists(filename)){
    #     if(overwritefile == FALSE ){
    #         stop("File already exists, please set overwrite = TRUE if you want to remove previous file")
    #     } else {
    #         message("File will be overwritten")
    #     }
    # } 
    
    # if( bdgetDim_hdf5(filename, Y)[1] < (m * 2) |  bdgetDim_hdf5(filename, X)[2] < (m * 2) ) {
    #     stop("m is too big for data size. Please set m to a small value")
    # }

    matrices <- c(X, Y)
    sapply( matrices, getQRbyBlocks, file = filename, mblocks = m, center = bcenter, scale = bscale, bcols = bycols, overwrt = overwriteResults )

    # Step 7
    #   tQXQY <- crossprod(t(QX), QY)[1:ncol(x), ]
    res <- bdCrossprod_hdf5(filename = filename, 
                            group = "Step6", A = "XQ",
                            groupB = "Step6", B = "YQ", 
                            outgroup = "Step7")
    
    # Step 8 : 
    # z <- svd( tQXQY )
    res <- bdSVD_hdf5(file = filename, 
                      group = "Step7", dataset = "CrossProd_XQxYQ",
                      bcenter = FALSE, bscale = FALSE, k = 16, q = 2, threads = 3)
    
    res <- sapply( matrices, bdgetDim_hdf5, filename = filename )
    
    writeCCAComponents_hdf5( filename, res[2,X], res[2,Y])
    
    if( keepInteResults == FALSE){
        sapply(paste0 ("Step",1:7), function (x) {
            invisible(bdRemove_hdf5_element( filename, element = x))
            print(paste0 ("Step",x, "Removed"))
        })
    }
    
}

