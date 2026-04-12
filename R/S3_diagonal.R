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
#' \code{\link{HDF5Matrix}} and a diagonal vector (a 1-row or 1-column
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
#' @param \dots Additional arguments passed to \code{x$diag_op()}:
#'   \code{outgroup}, \code{outdataset}, \code{overwrite},
#'   \code{threads}, \code{compression}.
#' @return A new \code{HDF5Matrix}.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' M   <- hdf5_create_matrix(tmp, "data/M", data = matrix(1:12, 3, 4))
#' d   <- hdf5_create_matrix(tmp, "data/d", data = matrix(c(1,2,3,4), 1, 4))
#' R1  <- diag_op(M, d, "*")    # scale each column
#' R2  <- M * d                  # same via operator auto-dispatch
#' close(M); close(d); close(R1); close(R2)
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{diag_scale}}
#' @export
diag_op <- function(x, diag, op = "+", ...) UseMethod("diag_op")

#' @rdname diag_op
#' @export
diag_op.HDF5Matrix <- function(x, diag, op = "+", ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")
    x$diag_op(diag, op = op, ...)
}


# ── diag_scale() ─────────────────────────────────────────────────────────────

#' Scalar diagonal operation on an HDF5Matrix
#'
#' @description
#' Applies a scalar arithmetic operation to the diagonal elements of an
#' \code{\link{HDF5Matrix}}. Off-diagonal elements are not modified.
#' Delegates to \code{bdDiag_scalar_hdf5()}.
#'
#' @param x      An \code{HDF5Matrix}.
#' @param scalar Numeric scalar.
#' @param op     Operation: \code{"add"}, \code{"subtract"},
#'   \code{"multiply"} (default), or \code{"divide"}.
#' @param \dots  Additional arguments passed to \code{x$diag_scale()}:
#'   \code{outgroup}, \code{outdataset}, \code{overwrite},
#'   \code{threads}, \code{compression}.
#' @return A new \code{HDF5Matrix}.
#'
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".h5")
#' M   <- hdf5_create_matrix(tmp, "data/M", data = diag(5))
#' R   <- diag_scale(M, scalar = 3, op = "multiply")
#' close(M); close(R)
#' unlink(tmp)
#' }
#'
#' @seealso \code{\link{diag_op}}
#' @export
diag_scale <- function(x, scalar, op = "multiply", ...) UseMethod("diag_scale")

#' @rdname diag_scale
#' @export
diag_scale.HDF5Matrix <- function(x, scalar, op = "multiply", overwrite = FALSE, ...) {
    if (!x$is_valid()) stop("HDF5Matrix is closed or invalid")
    x$diag_scale(scalar = scalar, op = op, overwrite = overwrite, ...)
}


# ── Ops.HDF5Matrix — extended with diagonal auto-detection ───────────────────
#
# Supersedes the definition in S3_arithmetic.R.
# When one operand is a 1-row or 1-column HDF5Matrix the diagonal path
# is used; otherwise the standard element-wise path is used.

#' @export
Ops.HDF5Matrix <- function(e1, e2) {

    if (!inherits(e1, "HDF5Matrix") || !inherits(e2, "HDF5Matrix"))
        stop("Both operands must be HDF5Matrix objects")
    if (!e1$is_valid()) stop("e1 is closed or invalid")
    if (!e2$is_valid()) stop("e2 is closed or invalid")

    d1 <- dim(e1)
    d2 <- dim(e2)
    e1_is_diag <- (d1[1] == 1L || d1[2] == 1L)
    e2_is_diag <- (d2[1] == 1L || d2[2] == 1L)

    if (e1_is_diag || e2_is_diag) {
        mat  <- if (e1_is_diag) e2 else e1
        diag <- if (e1_is_diag) e1 else e2
        return(mat$diag_op(diag,
                           op          = .Generic,
                           threads     = .get_option("threads"),
                           compression = .get_option("compression",
                                                     default = NULL)))
    }

    switch(.Generic,
        "+" = e1$add(e2,
                     paral       = .get_option("paral"),
                     block_size  = .get_option("block_size"),
                     threads     = .get_option("threads"),
                     compression = .get_option("compression", default = NULL)),

        "-" = e1$subtract(e2,
                          paral       = .get_option("paral"),
                          block_size  = .get_option("block_size"),
                          threads     = .get_option("threads"),
                          compression = .get_option("compression", default = NULL)),

        "*" = e1$multiply_ew(e2,
                             paral       = .get_option("paral"),
                             block_size  = .get_option("block_size"),
                             threads     = .get_option("threads"),
                             compression = .get_option("compression", default = NULL)),

        "/" = e1$divide_ew(e2,
                           paral       = .get_option("paral"),
                           block_size  = .get_option("block_size"),
                           threads     = .get_option("threads"),
                           compression = .get_option("compression", default = NULL)),

        stop(paste0("operator '", .Generic,
                    "' is not supported for HDF5Matrix objects"))
    )
}
