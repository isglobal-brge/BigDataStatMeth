#' Close all HDF5Matrix objects
#' 
#' @description
#' Finds and closes all HDF5Matrix objects in the specified environment.
#' 
#' @param envir Environment to search (default: .GlobalEnv)
#' @param verbose Show details (default: TRUE)
#' 
#' @return Invisible vector of closed filenames
#' 
#' @details
#' This function:
#' \itemize{
#'   \item Searches for HDF5Matrix objects in the environment
#'   \item Calls \code{$close()} on each valid object
#'   \item Forces garbage collection
#'   \item Reports closed files
#' }
#' 
#' **Note:** Only finds objects in the specified environment.
#' Objects inside functions or other environments are not affected.
#' 
#' @examples
#' \donttest{
#' tmp1 <- tempfile(fileext = ".h5")
#' tmp2 <- tempfile(fileext = ".h5")
#' rhdf5::h5createFile(tmp1)
#' rhdf5::h5createFile(tmp2)
#' rhdf5::h5write(matrix(1:100, 10, 10), tmp1, "data/A")
#' rhdf5::h5write(matrix(1:100, 10, 10), tmp2, "data/B")
#' 
#' X <- hdf5_matrix(tmp1, "data/A")
#' Y <- hdf5_matrix(tmp2, "data/B")
#' 
#' # Both open
#' X$is_valid()  # TRUE
#' Y$is_valid()  # TRUE
#' 
#' # Close all at once
#' hdf5_close_all()
#' 
#' # Both closed
#' X$is_valid()  # FALSE
#' Y$is_valid()  # FALSE
#' 
#' # Cleanup
#' unlink(c(tmp1, tmp2))
#' }
#' 
#' @export
hdf5_close_all <- function(envir = .GlobalEnv, verbose = TRUE) {
    obj_names <- ls(envir = envir)
    closed_count <- 0
    closed_files <- character(0)
    
    for (obj_name in obj_names) {
        obj <- tryCatch({
            get(obj_name, envir = envir)
        }, error = function(e) NULL)
        
        if (is.null(obj)) next
        
        if (inherits(obj, "HDF5Matrix")) {
            if (obj$is_valid()) {
                info <- tryCatch({
                    obj$info()
                }, error = function(e) list(filename = "unknown"))
                
                closed_files <- c(closed_files, info$filename)
                obj$close()
                closed_count <- closed_count + 1
            }
        }
    }
    
    # Close any remaining HDF5 handles not tracked as HDF5Matrix objects
    tryCatch(
        rcpp_hdf5_close_all_registry(),
        error = function(e) invisible(NULL)
    )
    
    # Force garbage collection
    gc()
    
    # Report
    if (verbose) {
        if (closed_count > 0) {
            message("Closed ", closed_count, " HDF5Matrix object(s)")
            unique_files <- unique(closed_files)
            if (length(unique_files) <= 5) {
                message("Files: ", paste(basename(unique_files), collapse = ", "))
            } else {
                message("Files: ", length(unique_files), " unique file(s)")
            }
        } else {
            message("No open HDF5Matrix objects found")
        }
    }
    
    invisible(unique(closed_files))
}


#' List open HDF5 files
#' 
#' @description
#' List all currently open HDF5 files in the R session.
#' Useful for debugging file locking issues.
#' 
#' @return Character vector of open file paths
#' 
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' rhdf5::h5createFile(tmp)
#' rhdf5::h5write(matrix(1:100, 10, 10), tmp, "data/matrix")
#' 
#' X <- hdf5_matrix(tmp, "data/matrix")
#' 
#' # Check open files
#' hdf5_list_open()  # Shows tmp file
#' 
#' X$close()
#' hdf5_list_open()  # Empty
#' 
#' unlink(tmp)
#' }
#' 
#' @export
hdf5_list_open <- function() {
    # Use rhdf5 to list open files
    open_files <- tryCatch({
        rhdf5::h5listIdentifier()
    }, error = function(e) {
        list()
    })
    
    if (length(open_files) == 0) {
        message("No open HDF5 files")
        return(character(0))
    }
    
    # Extract filenames
    file_ids <- open_files[names(open_files) == "type" & open_files == "H5I_FILE"]
    
    message(length(file_ids), " HDF5 file(s) open")
    
    return(names(file_ids))
}