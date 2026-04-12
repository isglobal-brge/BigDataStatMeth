# S3 methods - Elementwise arithmetic
# Ops.HDF5Matrix dispatches +, -, *, / to the R6 methods


#' Elementwise arithmetic operators for HDF5Matrix objects
#'
#' @description
#' Standard R arithmetic operators applied element-wise to \code{HDF5Matrix}
#' objects stored on disk. Both operands must be \code{HDF5Matrix} objects
#' with identical dimensions.
#'
#' @param e1 An \code{HDF5Matrix} object (left-hand side)
#' @param e2 An \code{HDF5Matrix} object (right-hand side)
#'
#' @return A new \code{HDF5Matrix} containing the result, stored in the same
#'   HDF5 file as \code{e1} under a temporary dataset name.
#'
#' @details
#' Supported operators:
#' \describe{
#'   \item{\code{+}}{Element-wise addition}
#'   \item{\code{-}}{Element-wise subtraction}
#'   \item{\code{*}}{Element-wise multiplication (Hadamard product)}
#'   \item{\code{/}}{Element-wise division. Division by zero produces
#'     \code{NaN} or \code{Inf}, matching base R behaviour.}
#' }
#'
#' All operations use block-wise processing and optional OpenMP parallelisation,
#' controlled via \code{\link{hdf5matrix_options}}.
#'
#' **Performance settings:**
#'
#' Global options set via \code{\link{hdf5matrix_options}} are applied.
#' For explicit control use the R6 methods directly:
#' \code{A$add(B, paral = TRUE, threads = 4)}.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' bdCreate_hdf5_matrix(tmp, matrix(1:12, 3, 4), "data", "A")
#' bdCreate_hdf5_matrix(tmp, matrix(2,    3, 4), "data", "B")
#'
#' A <- hdf5_matrix(tmp, "data/A")
#' B <- hdf5_matrix(tmp, "data/B")
#'
#' C <- A + B     # element-wise addition
#' D <- A - B     # element-wise subtraction
#' E <- A * B     # Hadamard product
#' F <- A / B     # element-wise division
#'
#' # Verify against R
#' all.equal(as.matrix(C), as.matrix(A) + as.matrix(B))
#'
#' A$close(); B$close(); C$close(); D$close(); E$close(); F$close()
#' unlink(tmp)
#' }
#'
#' @seealso
#' \code{\link{hdf5matrix_options}} for global performance settings,
#' \code{\link{HDF5Matrix}} for R6 methods with explicit parameters
#'
#' @export
Ops.HDF5Matrix <- function(e1, e2) {

    if (!inherits(e1, "HDF5Matrix") || !inherits(e2, "HDF5Matrix"))
        stop("Both operands must be HDF5Matrix objects")

    switch(.Generic,
        "+" = e1$add(e2,
                     paral      = .get_option("paral"),
                     block_size = .get_option("block_size"),
                     threads    = .get_option("threads"),
                     compression = .get_option("compression", default = NULL)),

        "-" = e1$subtract(e2,
                          paral      = .get_option("paral"),
                          block_size = .get_option("block_size"),
                          threads    = .get_option("threads"),
                          compression = .get_option("compression", default = NULL)),

        "*" = e1$multiply_ew(e2,
                             paral      = .get_option("paral"),
                             block_size = .get_option("block_size"),
                             threads    = .get_option("threads"),
                             compression = .get_option("compression", default = NULL)),

        "/" = e1$divide_ew(e2,
                           paral      = .get_option("paral"),
                           block_size = .get_option("block_size"),
                           threads    = .get_option("threads"),
                           compression = .get_option("compression", default = NULL)),

        stop(paste0("operator '", .Generic,
                    "' is not supported for HDF5Matrix objects"))
    )
}
