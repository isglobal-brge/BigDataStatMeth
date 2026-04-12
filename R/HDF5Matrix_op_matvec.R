# HDF5Matrix_op_matvec.R
#
# Adds five methods to the HDF5Matrix R6 class:
#   $sweep()       -- broadcast vector over matrix (like base R sweep)
#   $diag()        -- extract diagonal -> numeric vector in memory
#   $diag_assign() -- write diagonal in-place
#   $diag_op()     -- binary op on diagonals of two HDF5Matrix objects
#   $diag_scale()  -- apply scalar to diagonal elements


# -- $sweep() -----------------------------------------------------------------

HDF5Matrix$set("public", "sweep",
# @description
# Broadcasts a 1-row or 1-column HDF5Matrix (vector) across the rows or
# columns of this matrix.  Equivalent to base R sweep().
#
# @param v       An HDF5Matrix with exactly one row or one column.
#   Must be in the same HDF5 file.
#   1-row (1xN) and byrows = FALSE: N must equal ncol(self).
#   1-row (1xM) and byrows = TRUE:  M must equal nrow(self).
# @param FUN     Character. Operation: "+", "-", "*" (default), "/", "pow".
# @param byrows  Logical. FALSE (default): vector applied along each
#   row (scales columns). TRUE: vector applied along each column (scales rows).
# @param paral      Logical or NULL. OpenMP parallelisation.
# @param threads    Integer or NULL. Thread count.
# @param compression Integer (0-9) or NULL. gzip level for result.
# @return A new HDF5Matrix with the broadcast result.
function(v,
         FUN         = "*",
         byrows      = FALSE,
         paral       = NULL,
         threads     = NULL,
         compression = NULL) {

    if (!self$is_valid())     stop("HDF5Matrix is closed or invalid")
    if (!inherits(v, "HDF5Matrix")) stop("v must be an HDF5Matrix")
    if (!v$is_valid())        stop("v is closed or invalid")

    d <- dim(v)
    if (d[1] != 1L && d[2] != 1L)
        stop("sweep: v must be a 1-row or 1-column HDF5Matrix")

    if (private$filename != v$get_filename())
        stop("sweep: both HDF5Matrix objects must be in the same HDF5 file")

    FUN <- match.arg(FUN, c("+", "-", "*", "/", "pow"))

    compression_eff <- .get_option("compression", default = NULL, override = compression)
    paral_eff   <- .get_option("paral",   default = NULL, override = paral)
    threads_eff <- .get_option("threads", default = NULL, override = threads)

    v_ptr <- v$.__enclos_env__$private$ptr

    res <- rcpp_hdf5dataset_sweep(
        ptr_mat     = private$ptr,
        ptr_vec     = v_ptr,
        func        = FUN,
        byrows      = isTRUE(byrows),
        paral       = paral_eff,
        threads     = threads_eff,
        compression = compression_eff
    )

    hdf5_matrix(res$filename, res$path)
})


# -- $diag() ------------------------------------------------------------------

HDF5Matrix$set("public", "diag",
# @description
# Extracts the diagonal elements of this HDF5Matrix and returns them as an
# in-memory numeric vector (length = min(nrow, ncol)).
#
# @return Numeric vector of diagonal elements.
function() {
    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")
    rcpp_hdf5dataset_diag_get(private$ptr)
})


# -- $diag_assign() -----------------------------------------------------------

HDF5Matrix$set("public", "diag_assign",
# @description
# Replaces the diagonal elements of this HDF5Matrix in-place with values.
# The matrix must be square.  The HDF5 file is modified.
#
# @param values Numeric vector of length min(nrow, ncol).
# @return Invisibly self (for chaining).
function(values) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")
    if (!is.numeric(values)) stop("values must be a numeric vector")

    # rcpp_hdf5dataset_diag_set opens its own fresh HDF5 handle internally,
    # so we can pass the active ptr directly without closing it first.
    rcpp_hdf5dataset_diag_set(private$ptr, as.numeric(values))
    invisible(self)
})


# -- $diag_op() ---------------------------------------------------------------

HDF5Matrix$set("public", "diag_op",
# @description
# Element-wise binary operation on the diagonal elements of this HDF5Matrix
# and another HDF5Matrix (other).  Inputs may be square matrices
# (diagonal extracted) or diagonal vectors (1xN or Nx1, used directly).
# Result is a new HDF5Matrix (vector form).
#
# @param other       An HDF5Matrix. Must be in the same HDF5 file.
# @param op          Character. "+" (default), "-", "*", "/".
# @param paral       Logical or NULL.
# @param threads     Integer or NULL.
# @param compression Integer (0-9) or NULL.
# @return A new HDF5Matrix containing the diagonal result.
function(other,
         op          = "+",
         paral       = NULL,
         threads     = NULL,
         compression = NULL) {

    if (!self$is_valid())         stop("HDF5Matrix is closed or invalid")
    if (!inherits(other, "HDF5Matrix")) stop("other must be an HDF5Matrix")
    if (!other$is_valid())        stop("other is closed or invalid")
    if (private$filename != other$get_filename())
        stop("diag_op: both objects must be in the same HDF5 file")

    op <- match.arg(op, c("+", "-", "*", "/"))

    compression_eff <- .get_option("compression", default = NULL,
                                   override = compression)
    paral_eff   <- .get_option("paral",   default = NULL, override = paral)
    threads_eff <- .get_option("threads", default = NULL, override = threads)

    other_ptr <- other$.__enclos_env__$private$ptr

    res <- rcpp_hdf5dataset_diag_op(
        ptr_a       = private$ptr,
        ptr_b       = other_ptr,
        op          = op,
        paral       = paral_eff,
        threads     = threads_eff,
        compression = compression_eff
    )

    hdf5_matrix(res$filename, res$path)
})


# -- $diag_scale() ------------------------------------------------------------

HDF5Matrix$set("public", "diag_scale",
# @description
# Applies a scalar arithmetic operation to the diagonal elements of this
# HDF5Matrix.  Off-diagonal elements are unchanged.  The input must be
# square (or a diagonal vector 1xN / Nx1).  Result is a new HDF5Matrix.
#
# @param scalar      Numeric scalar.
# @param op          Character. "multiply" (default), "add", "subtract", "divide".
# @param paral       Logical or NULL.
# @param threads     Integer or NULL.
# @param compression Integer (0-9) or NULL.
# @return A new HDF5Matrix with scaled diagonal.
function(scalar,
         op          = "multiply",
         overwrite   = FALSE, 
         paral       = NULL,
         threads     = NULL,
         compression = NULL) {

    if (!self$is_valid()) stop("HDF5Matrix is closed or invalid")
    if (!is.numeric(scalar) || length(scalar) != 1L)
        stop("scalar must be a single numeric value")

    op <- match.arg(op, c("add", "subtract", "multiply", "divide"))
    op_code <- switch(op, add = 0L, subtract = 1L, multiply = 2L, divide = 3L)

    compression_eff <- .get_option("compression", default = NULL,
                                   override = compression)
    paral_eff   <- .get_option("paral",   default = NULL, override = paral)
    threads_eff <- .get_option("threads", default = NULL, override = threads)

    res <- rcpp_hdf5dataset_diag_scale(
        ptr_mat     = private$ptr,
        scalar      = as.double(scalar),
        op_code     = op_code,
        paral       = paral_eff,
        threads     = threads_eff,
        compression = compression_eff
    )

    hdf5_matrix(res$filename, res$path)
})
