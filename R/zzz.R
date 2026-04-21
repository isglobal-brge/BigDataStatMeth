#' Package hook: cleanup on unload
#'
#' @description
#' Closes all open HDF5Dataset C++ objects and HDF5 file handles
#' when the package is unloaded. This prevents finalizers from running
#' after the C++ library is gone, which would cause crashes.
#'
#' @param libpath Library path (unused)
#'
#' @keywords internal
.onUnload <- function(libpath) {
    # First: close all C++ objects via the registry
    # This must happen BEFORE the shared library is unloaded
    tryCatch(
        rcpp_hdf5_close_all_registry(),
        error = function(e) invisible(NULL)
    )
    
    # tryCatch(
    #     rcpp_hdf5_close_all_handles_unsafe(), 
    #     error = function(e) invisible(NULL)
    # )
    # Force GC to trigger any remaining finalizers while C++ is still loaded
    gc()
}