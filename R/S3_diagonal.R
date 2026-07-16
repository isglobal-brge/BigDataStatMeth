# S3_diagonal.R
#
# S3 generics diag_op() and diag_scale() for HDF5Matrix objects.
# Also updates Ops.HDF5Matrix to auto-detect diagonal operands and
# route them through $diag_op() instead of the element-wise path.
#
# Replaces Ops.HDF5Matrix defined in S3_arithmetic.R.
# R loads files alphabetically, so S3_diagonal.R loads after
# S3_arithmetic.R and its definition takes effect.

# ── diag_op() ────────────────────────────────────────────────────────────────

#' Diagonal-vector operation on an HDF5Matrix
#'
#' @description
#' Applies an element-wise binary operation between an
#' \code{HDF5Matrix} and a diagonal vector (a 1-row or 1-column
#' \code{HDF5Matrix}). The vector is broadcast across each row of the
#' matrix.
#'
#' The standard arithmetic operators (\code{+}, \code{-}, \code{*},
#' \code{/}) dispatch automatically to this function when one operand is
#' a 1-row or 1-column \code{HDF5Matrix}.
#'
#' @param x    An \code{HDF5Matrix} (matrix, \eqn{m \times n}).
#' @param diag An \code{HDF5Matrix} with one row or one column.
#' @param op   Character. One of \code{"+"}, \code{"-"}, \code{"*"},
#'   \code{"/"}.
#' @param outgroup  Character or \code{NULL}. HDF5 group where the result is
#'   stored. Default \code{"OUTPUT"}.
#' @param outdataset Character or \code{NULL}. Dataset name for the result.
#' @param \dots Additional arguments passed to \code{x$diag_op()}:
#'   \code{outgroup}, \code{outdataset}, \code{overwrite},
#'   \code{threads}, \code{compression}.
#' @return A new \code{HDF5Matrix}.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' M   <- hdf5_create_matrix(tmp, "data/M", data = matrix(rnorm(10000), 100, 100))
#' d   <- hdf5_create_matrix(tmp, "data/d", data = matrix(rnorm(10000), 100, 100))
#' R1  <- diag_op(M, d, "*")    # scale each column
#' R2  <- M * d                  # same via operator auto-dispatch
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{diag_scale}}
#' @export
diag_op <- function(x, diag, op = "+", ...) UseMethod("diag_op")

#' @rdname diag_op
#' @export
diag_op.HDF5Matrix <- function(x, diag, 
                               op = "+", 
                               outgroup   = NULL,
                               outdataset = NULL, ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")
    x$diag_op(diag, op = op, ...)
}


# ── diag_scale() ─────────────────────────────────────────────────────────────

#' Scalar diagonal operation on an HDF5Matrix
#'
#' @description
#' Applies a scalar arithmetic operation to the diagonal elements of an
#' \code{HDF5Matrix}. Off-diagonal elements are not modified.
#' Delegates to \code{bdDiag_scalar_hdf5()}.
#'
#' @param x      An \code{HDF5Matrix}.
#' @param scalar Numeric scalar.
#' @param op     Operation: \code{"add"}, \code{"subtract"},
#'   \code{"multiply"} (default), or \code{"divide"}.
#' @param outgroup  Character or \code{NULL}. HDF5 group where the result is
#'   stored. Default \code{"OUTPUT"}.
#' @param outdataset Character or \code{NULL}. Dataset name for the result.
#' @param \dots  Additional arguments passed to \code{x$diag_scale()}:
#'   \code{outgroup}, \code{outdataset}, \code{overwrite},
#'   \code{threads}, \code{compression}.
#' @param overwrite Logical. If \code{TRUE}, overwrite the existing result
#'   dataset. Default \code{FALSE}.
#' @return A new \code{HDF5Matrix}.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' M   <- hdf5_create_matrix(tmp, "data/M", data = diag(5))
#' R   <- diag_scale(M, scalar = 3, op = "multiply")
#' hdf5_close_all()
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{diag_op}}
#' @export
diag_scale <- function(x, scalar, op = "multiply", ...) UseMethod("diag_scale")

#' @rdname diag_scale
#' @export
diag_scale.HDF5Matrix <- function(x, scalar, 
                                  op = "multiply", 
                                  outgroup   = NULL,
                                  outdataset = NULL,
                                  overwrite = FALSE, ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")
    x$diag_scale(scalar = scalar, op = op, overwrite = overwrite, ...)
}


# ── Ops.HDF5Matrix ───────────────────────────────────────────────────────────
#
# The single definition of Ops.HDF5Matrix (element-wise arithmetic plus
# 1xn / nx1 diagonal-broadcast auto-detection) lives in S3_arithmetic.R,
# together with its roxygen documentation. It calls diag_op() (defined above)
# for the broadcast path.
