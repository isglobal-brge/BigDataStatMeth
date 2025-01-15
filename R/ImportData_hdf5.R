#' Import data from url or a file
#' 
#' This function download data from an url and decompress data (if needed), 
#' then imports the file to hdf5 data file
#' 
#' @export
#' 
#' @param inFile string file name or url with data to import
#' @param destFile file name and path to store imported data
#' @param destGroup group name to store the dataset
#' @param destDataset dataset name to store the input file in hdf5
#' @param header (optional) either a logical value indicating whether the column 
#' names of x are to be written along with x, or a character vector of column 
#' names to be written. See the section on ‘CSV files’ for the meaning of 
#' col.names = NA.
#' @param rownames (optional) either a logical value indicating whether the row 
#' names of x are to be written along with x, or a character vector of row names 
#' to be written.
#' @param overwrite (optional) either a logical value indicating whether the 
#' output file can be overwritten or not.
#' @param overwriteFile logical (optional), CAUTION, if TRUE, file will be 
#' overwritten with imported dataset, by default `fileoverwrite = FALSE` to avoid
#' file overwritting.
#' @param sep (optional), by default = "\\t". The field separator string. Values 
#' within each row of x are separated by this string.
#' @param paral, (optional, default = TRUE) if paral = TRUE performs parallel 
#' computation else performs seria computation
#' @param threads (optional) only if bparal = true, number of concurrent 
#' threads in parallelization if threads is null then threads =  maximum 
#' number of threads available
#' @examples
#'    print ("Example in vignette")
#'    
#' @return none value returned, data are stored in a dataset inside an hdf5 data file.
bdImportData_hdf5 <- function( inFile, destFile, destGroup, destDataset, 
                               header = TRUE, rownames = FALSE, 
                               overwrite = FALSE, overwriteFile = FALSE, 
                               sep = NULL,  paral = NULL, threads = NULL)
{
    
    untarExtensions <- c("tar.gz", "gzip", "bzip2", "gz", "tgz")
    
    extension <- substr(inFile, regexpr("\\.[^\\.]*$", inFile)[1]+1, nchar(inFile))
    filename <- substr(inFile, regexpr("\\/[^\\/]*$", inFile)[1]+1, nchar(inFile)-nchar(extension)-1)
    importfile <- ""
    
    # if( url.exists(inFile))
    # {
    dfile <- try (download.file(url = inFile, destfile = paste0(getwd(), "/", filename, ".", extension)), TRUE)
    if (inherits(dfile, "try-error")) {
        if(!file.exists(inFile)) {
            stop("File does not exists, please review the route")
        } 
        
        importfile <- inFile
        
    } else {
        inFile <- paste0(filename,".", extension)
    }
    
    
    if(extension == "zip") {
        importfile <- unzip(inFile, list = TRUE )$Name[1]
        unzip(inFile)
    }else if(extension %in% untarExtensions) {
        importfile <- untar(inFile, list = TRUE )[1]
        untar(inFile)
    } else {
        importfile <- inFile
    }
    
    # Import files to hdf5
    bdImportTextFile_hdf5(filename = importfile,
                          outputfile = destFile,
                          outGroup = destGroup, 
                          outDataset = destDataset,
                          sep = sep,
                          header = header,
                          rownames = rownames,
                          overwrite = overwrite,
                          overwriteFile = overwriteFile,
                          paral = paral, 
                          threads = threads)
    # unlink(importfile)
}