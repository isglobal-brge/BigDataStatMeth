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
#' 
#' X  <- hdf5_create_matrix(tmp1, "data/A", data = matrix(rnorm(100), 10, 10))
#' Y  <- hdf5_create_matrix(tmp2, "data/B", data = matrix(rnorm(100), 10, 10))
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
    
    # Safety net: close any handles that escaped the registry
    # Closes datasets/groups/attrs per file, then file handle — releases lock
    tryCatch(
        rcpp_hdf5_close_all_file_handles(),
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


#' Close all HDF5 handles for a specific file
#'
#' @description
#' Closes all open \code{HDF5Matrix} objects and HDF5 C library handles
#' associated with a single HDF5 file, without affecting other open files.
#'
#' @param x An \code{HDF5Matrix} object, or a character string with the
#'   path to the HDF5 file.
#'
#' @return Invisibly, the absolute path of the closed file.
#'
#' @examples
#' \donttest{
#' fn1 <- tempfile(fileext = ".h5")
#' fn2 <- tempfile(fileext = ".h5")
#' A <- hdf5_create_matrix(fn1, "data/A", data = matrix(1:9, 3, 3))
#' B <- hdf5_create_matrix(fn2, "data/B", data = matrix(1:9, 3, 3))
#'
#' # Close only fn1 — B remains open and usable
#' hdf5_close_file(fn1)
#' dim(B)   # still works
#'
#' hdf5_close_all()
#' unlink(c(fn1, fn2))
#' }
#'
#' @seealso \code{\link{hdf5_close_all}} to close all files at once.
#' @export
hdf5_close_file <- function(x) {
    if (inherits(x, "HDF5Matrix")) {
        filename <- normalizePath(x$get_filename(), mustWork = FALSE)
    } else if (is.character(x) && length(x) == 1) {
        filename <- normalizePath(x, mustWork = FALSE)
    } else {
        stop("x must be an HDF5Matrix or a single file path")
    }
    
    # Close all HDF5Matrix objects pointing to this file in global env
    for (obj_name in ls(envir = .GlobalEnv)) {
        obj <- tryCatch(get(obj_name, envir = .GlobalEnv), error = function(e) NULL)
        if (inherits(obj, "HDF5Matrix") && obj$is_valid()) {
            if (identical(normalizePath(obj$get_filename(), mustWork = FALSE),
                          filename)) {
                obj$close()
            }
        }
    }
    
    # Close remaining handles at C++ and HDF5 library level
    tryCatch(
        rcpp_hdf5_close_file_handles(filename),
        error = function(e) invisible(NULL)
    )
    
    gc()
    invisible(filename)
}


# #' List open HDF5 files
# #' 
# #' @description
# #' List all currently open HDF5 files in the R session.
# #' Useful for debugging file locking issues.
# #' 
# #' @return Character vector of open file paths
# #' 
# #' @examples
# #' \donttest{
# #' tmp <- tempfile(fileext = ".h5")
# #' 
# #' X <- hdf5_create_matrix(tmp, "data/matrix", data = matrix(1:100, 10, 10))
# #' X <- hdf5_matrix(tmp, "data/matrix")
# #' 
# #' # Check open files
# #' hdf5_list_open()  # Shows tmp file
# #' 
# #' X$close()
# #' hdf5_list_open()  # Empty
# #' 
# #' unlink(tmp)
# #' }
# #' 
# #' @export
# hdf5_list_open <- function() {
#     ids <- tryCatch({
#         rcpp_hdf5_list_open_files()
#     }, error = function(e) character(0))
#     
#     if (length(ids) == 0) {
#         message("No open HDF5 files")
#         return(character(0))
#     }
#     ids
# }


# hdf5_list_open <- function() {
#     # Use rhdf5 to list open files
#     open_files <- tryCatch({
#         rhdf5::h5listIdentifier()
#     }, error = function(e) {
#         list()
#     })
#     
#     if (length(open_files) == 0) {
#         message("No open HDF5 files")
#         return(character(0))
#     }
#     
#     # Extract filenames
#     file_ids <- open_files[names(open_files) == "type" & open_files == "H5I_FILE"]
#     
#     message(length(file_ids), " HDF5 file(s) open")
#     
#     return(names(file_ids))
# }