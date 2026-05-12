# R/HDF5Matrix_import.R
# Modern wrapper for data importation

#' Import data from file or URL into HDF5 format
#'
#' @description
#' Modern wrapper for importing CSV, TSV, or other delimited text files into
#' HDF5 format. Returns an \code{HDF5Matrix} object ready for use.
#'
#' @param source Character. Path to local file or URL to import.
#'   Supports compressed files (.gz, .tar.gz, .zip, .bz2).
#' @param filename Character. Path to HDF5 output file (created if doesn't exist).
#' @param dataset Character. Full dataset path (e.g., "data/imported" or "group/dataset").
#' @param sep Character. Field separator. Default \code{NULL} (auto-detect from extension:
#'   "," for .csv, "\\t" for .tsv, "\\t" otherwise).
#' @param header Logical or character vector. If \code{TRUE}, first row contains column
#'   names. If character vector, use these as column names. Default \code{TRUE}.
#' @param rownames Logical or character vector. If \code{TRUE}, first column contains
#'   row names. If character vector, use these as row names. Default \code{FALSE}.
#' @param overwrite Logical. If \code{TRUE}, overwrite dataset if exists.
#'   Default \code{FALSE}.
#' @param parallel Logical. Use parallel processing for import. Default \code{TRUE}.
#' @param threads Integer. Number of threads for parallel processing. Default \code{NULL}
#'   (uses all available cores).
#'
#' @return \code{HDF5Matrix} object pointing to the imported data.
#'
#' @details
#' This function is a modern, user-friendly wrapper around \code{\link{bdImportData_hdf5}}
#' and \code{\link{bdImportTextFile_hdf5}}. It:
#' 
#' \itemize{
#'   \item Automatically detects file format from extension
#'   \item Handles compressed files (.gz, .tar.gz, .zip)
#'   \item Downloads from URLs automatically
#'   \item Returns ready-to-use \code{HDF5Matrix} object
#'   \item Uses sensible defaults for most use cases
#' }
#'
#' **Supported formats:**
#' \itemize{
#'   \item CSV files (.csv) - comma-separated
#'   \item TSV files (.tsv, .txt) - tab-separated
#'   \item Compressed files (.gz, .tar.gz, .zip, .bz2)
#'   \item Remote files (http://, https://, ftp://)
#' }
#'
#' **Memory efficiency:**
#' Import is done in a streaming fashion, so very large files can be imported
#' without loading them entirely into memory.
#'
#' @examples
#' \donttest{
#' csv_file  <- tempfile(fileext = ".csv")
#' hdf5_file <- tempfile(fileext = ".h5")
#' 
#' # Write sample numeric data
#' write.table(matrix(rnorm(50), nrow = 10, ncol = 5),
#'             csv_file, sep = ",", row.names = FALSE, col.names = TRUE)
#' 
#' # Import CSV to HDF5
#' mat <- hdf5_import(
#'   source   = csv_file,
#'   filename = hdf5_file,
#'   dataset  = "raw/data",
#'   sep      = ","
#' )
#' dim(mat)
#' 
#' hdf5_close_all()
#' unlink(c(csv_file, hdf5_file))
#' }
#'
#' @seealso
#' \code{\link{bdImportData_hdf5}} for the underlying implementation,
#' \code{\link{hdf5_create_matrix}} for creating matrices from R objects
#'
#' @export
hdf5_import <- function(source,
                        filename,
                        dataset,
                        sep = NULL,
                        header = TRUE,
                        rownames = FALSE,
                        overwrite = FALSE,
                        parallel = TRUE,
                        threads = NULL) {
  
  # Validate inputs
  if (!is.character(source) || length(source) != 1)
    stop("source must be a single character string")
  if (!is.character(filename) || length(filename) != 1)
    stop("filename must be a single character string")
  if (!is.character(dataset) || length(dataset) != 1)
    stop("dataset must be a single character string")
  if (!is.logical(overwrite) || length(overwrite) != 1)
    stop("overwrite must be TRUE or FALSE")
  
  # Parse dataset path: "group/subgroup/name" → (group_path, dataset_name)
  parts <- strsplit(dataset, "/")[[1]]
  if (length(parts) < 2) {
    stop("dataset must be in format 'group/dataset' or 'group/subgroup/dataset'")
  }
  
  dataset_name <- parts[length(parts)]
  group_path <- paste(parts[-length(parts)], collapse = "/")
  
  # Auto-detect separator from file extension if not provided
  if (is.null(sep)) {
    ext <- tools::file_ext(source)
    sep <- switch(tolower(ext),
      "csv" = ",",
      "tsv" = "\t",
      "\t"  # Default: tab-separated
    )
  }
  
  # Check if source is a URL
  is_url <- grepl("^(http://|https://|ftp://)", source, ignore.case = TRUE)
  
  # Check source file exists (only for local files)
  if (!is_url && !file.exists(source)) {
      stop("Source file not found: ", source)
  }
  
  # Use appropriate import function
  if (is_url || grepl("\\.(gz|tar\\.gz|zip|bz2|tgz)$", source, ignore.case = TRUE)) {
    # Use bdImportData_hdf5 for URLs and compressed files
    # (handles download and decompression)
    result <- bdImportData_hdf5(
      inFile = source,
      destFile = filename,
      destGroup = group_path,
      destDataset = dataset_name,
      header = header,
      rownames = rownames,
      overwrite = overwrite,
      overwriteFile = FALSE,  # Never overwrite entire file
      sep = sep,
      paral = parallel,
      threads = threads
    )
  } else {
      
      # # Normalize paths (symlinks solver)
      # source <- normalizePath(source, mustWork = TRUE)
      # filename <- normalizePath(filename, mustWork = FALSE)  # FALSE porque aún no existe
      
    # Use bdImportTextFile_hdf5 for local files
    # (more efficient, no intermediate download)
    result <- bdImportTextFile_hdf5(
      filename = source,
      outputfile = filename,
      outGroup = group_path,
      outDataset = dataset_name,
      sep = sep,
      header = header,
      rownames = rownames,
      overwrite = overwrite,
      overwriteFile = FALSE,  # Never overwrite entire file
      paral = parallel,
      threads = threads
    )
  }
  
  # Open and return HDF5Matrix object
  full_path <- paste0(group_path, "/", dataset_name)
  
  message("Import complete. Opening HDF5Matrix...")
  mat <- hdf5_matrix(filename, full_path)
  
  message("Dataset dimensions: ", dim(mat)[1], " x ", dim(mat)[2])
  
  return(mat)
}


#' Import multiple files into HDF5
#'
#' @description
#' Imports multiple files into the same HDF5 file, each as a separate dataset.
#' Useful for batch importing related datasets.
#'
#' @param sources Character vector. Paths to files or URLs to import.
#' @param filename Character. Path to HDF5 output file.
#' @param datasets Character vector. Dataset paths for each source file.
#'   Must be same length as \code{sources}.
#' @param ... Additional arguments passed to \code{\link{hdf5_import}}
#'
#' @return Named list of \code{HDF5Matrix} objects, one for each imported file.
#'
#' @examples
#' \donttest{
#' # Create temporary CSV files
#' f1 <- tempfile(fileext = ".csv")
#' f2 <- tempfile(fileext = ".csv")
#' f3 <- tempfile(fileext = ".csv")
#' hdf5_file <- tempfile(fileext = ".h5")
#' 
#' write.table(matrix(rnorm(20), 4, 5), f1, sep = ",",
#'             row.names = FALSE, col.names = TRUE)
#' write.table(matrix(rnorm(20), 4, 5), f2, sep = ",",
#'             row.names = FALSE, col.names = TRUE)
#' write.table(matrix(rnorm(20), 4, 5), f3, sep = ",",
#'             row.names = FALSE, col.names = TRUE)
#' 
#' mats <- hdf5_import_multiple(
#'   sources  = c(f1, f2, f3),
#'   filename = hdf5_file,
#'   datasets = c("data/exp1", "data/exp2", "data/exp3"),
#'   sep      = ","
#' )
#' 
#' dim(mats$exp1)
#' 
#' hdf5_close_all()
#' unlink(c(f1, f2, f3, hdf5_file))
#' }
#'
#' @export
hdf5_import_multiple <- function(sources, filename, datasets, ...) {
  
  if (length(sources) != length(datasets)) {
    stop("sources and datasets must have the same length")
  }
  
  # Import each file
  results <- list()
  for (i in seq_along(sources)) {
    message("\nImporting ", i, " of ", length(sources), ": ", sources[i])
    
    results[[i]] <- hdf5_import(
      source = sources[i],
      filename = filename,
      dataset = datasets[i],
      ...
    )
  }
  
  # Name the results based on dataset names
  names(results) <- sapply(datasets, function(x) {
    parts <- strsplit(x, "/")[[1]]
    parts[length(parts)]  # Use final part as name
  })
  
  message("\n", length(sources), " files imported successfully.")
  
  return(results)
}
