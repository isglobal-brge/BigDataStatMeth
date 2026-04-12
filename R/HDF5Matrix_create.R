#' Create an HDF5 dataset and return an HDF5Matrix object
#'
#' @description
#' Creates a new HDF5 dataset (optionally writing data) and returns an
#' \code{HDF5Matrix} object pointing to it.
#'
#' @param filename Character. Path to the HDF5 file (created if it does not exist).
#' @param dataset Character. Full path inside the HDF5 file in
#'   \code{"group/dataset"} or \code{"group/subgroup/dataset"} format.
#' @param nrow Integer or NULL. Number of rows. Required when \code{data = NULL}.
#' @param ncol Integer or NULL. Number of columns. Required when \code{data = NULL}.
#' @param data Numeric matrix, integer matrix, or numeric vector, or NULL.
#'   When non-NULL, the data are written to the new dataset. When NULL,
#'   an empty (zero-filled) dataset of size \code{nrow x ncol} is created.
#' @param dtype Character. Element type: \code{"double"} (default),
#'   \code{"integer"}, or \code{"logical"}.
#' @param overwrite Logical. If \code{TRUE}, an existing dataset is replaced.
#' @param compression Integer (0-9) or NULL. gzip compression level.
#'   \code{NULL} uses the global option set by \code{\link{hdf5matrix_options}}
#'   (default 6). Use \code{0} to disable compression.
#'
#' @return An \code{HDF5Matrix} object pointing to the created dataset.
#'
#' @details
#' Replaces the legacy \code{bdCreate_hdf5_matrix()} / \code{bdCreate_hdf5_emptyDataset()}
#' calls in the R6+S3 interface. The legacy functions remain available for
#' backward compatibility.
#'
#' Row and column names stored in the \code{dimnames} attribute of \code{data}
#' are written to the HDF5 file automatically.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#'
#' # Create from matrix data
#' mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' X <- hdf5_create_matrix(tmp, "data/X", data = mat)
#' dim(X)  # 20 x 10
#'
#' # Create empty dataset
#' Y <- hdf5_create_matrix(tmp, "data/Y", nrow = 1000, ncol = 500)
#' dim(Y)  # 1000 x 500
#'
#' # No compression (useful for benchmarks or intermediate results)
#' Z <- hdf5_create_matrix(tmp, "data/Z", data = mat, compression = 0)
#'
#' X$close(); Y$close(); Z$close()
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{hdf5matrix_options}} to set global compression default.
#' @export
hdf5_create_matrix <- function(filename,
                               dataset,
                               nrow = NULL,
                               ncol = NULL,
                               data = NULL,
                               dtype = c("double", "integer", "logical"),
                               overwrite = FALSE,
                               compression = NULL) {

    # --- input validation ---
    if (!is.character(filename) || length(filename) != 1)
        stop("filename must be a single character string")
    if (!is.character(dataset) || length(dataset) != 1)
        stop("dataset must be a single character string")
    if (!is.logical(overwrite) || length(overwrite) != 1)
        stop("overwrite must be TRUE or FALSE")

    dtype <- match.arg(dtype)

    # Parse "group/subgroup/name" -> (group_path, dataset_name)
    parts <- strsplit(dataset, "/")[[1]]
    if (length(parts) < 2)
        stop("dataset must be in format 'group/dataset' or 'group/subgroup/dataset'")
    dataset_name <- parts[length(parts)]
    group_path   <- paste(parts[-length(parts)], collapse = "/")

    # --- resolve dimensions ---
    if (!is.null(data)) {
        if (is.vector(data) && !is.list(data)) {
            if (is.null(nrow) && is.null(ncol)) {
                data <- matrix(data, ncol = 1)
            } else if (!is.null(nrow) && !is.null(ncol)) {
                nrow <- as.integer(nrow)[1]
                ncol <- as.integer(ncol)[1]
                if (length(data) < nrow * ncol) {
                    data <- rep_len(data, nrow * ncol)
                } else if (length(data) > nrow * ncol) {
                    warning("data length truncated to nrow*ncol")
                    data <- data[seq_len(nrow * ncol)]
                }
                data <- matrix(data, nrow = nrow, ncol = ncol)
            } else {
                stop("If data is a vector with nrow or ncol, both must be specified")
            }
        } else if (is.matrix(data)) {
            if (!is.null(nrow) && !is.null(ncol)) {
                nrow <- as.integer(nrow)[1]
                ncol <- as.integer(ncol)[1]
                if (nrow(data) != nrow || ncol(data) != ncol)
                    stop("data dimensions (", nrow(data), "x", ncol(data),
                         ") don't match nrow/ncol (", nrow, "x", ncol, ")")
            }
        } else {
            stop("data must be a matrix or vector")
        }
        nrow_final <- nrow(data)
        ncol_final <- ncol(data)
    } else {
        if (is.null(nrow) || is.null(ncol))
            stop("If data is NULL, nrow and ncol must be specified")
        if (!is.numeric(nrow) || !is.numeric(ncol))
            stop("nrow and ncol must be numeric")
        nrow_final <- as.integer(nrow)[1]
        ncol_final <- as.integer(ncol)[1]
        if (nrow_final < 1 || ncol_final < 1)
            stop("nrow and ncol must be >= 1")
    }

    # --- compression: override > global option > 6 ---
    compression_eff <- as.integer(.get_option("compression", default = 6L,
                                               override = compression))

    # --- dtype mapping ---
    dtype_bd <- switch(dtype,
                       "double"  = "real",
                       "integer" = "int",
                       "logical" = "logical",
                       "real")

    # --- delegate to C++ R6 wrapper ---
    rcpp_hdf5_create_matrix(
        filename          = filename,
        group             = group_path,
        dataset           = dataset_name,
        nrows             = nrow_final,
        ncols             = ncol_final,
        data              = data,
        dtype             = dtype_bd,
        overwrite_file    = FALSE,
        overwrite_dataset = overwrite,
        compression       = compression_eff
    )

    # Return HDF5Matrix object
    hdf5_matrix(filename, paste0(group_path, "/", dataset_name))
}
