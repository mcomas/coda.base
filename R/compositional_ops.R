#' Perturbation of compositional data
#'
#' The perturbation operation combines two compositions by component-wise
#' multiplication and then applies closure to ensure the result remains a
#' valid composition.
#'
#' @param X A numeric vector, matrix or data.frame containing compositions.
#' @param Y A numeric vector, matrix or data.frame with the same number of
#'   parts as \code{X}. If one input is a matrix or data.frame and the other is
#'   a vector, the vector is applied to every row. If both inputs are matrices
#'   or data.frames, they must have the same dimensions, or \code{Y} may contain
#'   a single composition to be applied to all rows of \code{X}.
#' @return An object with the same format as \code{X} containing the perturbed
#'   compositions, except that vector \code{X} with matrix or data.frame
#'   \code{Y} returns the same rectangular format as \code{Y}.
#'
#' @details
#' Perturbation is the analogue of addition in the simplex. Each part of
#' \code{X} is multiplied by the corresponding part of \code{Y}, and the result
#' is closed with \code{\link{closure}} so that each composition has constant
#' sum.
#'
#' @examples
#' x <- c(a = 1, b = 2, c = 3)
#' y <- c(a = 1, b = 1, c = 2)
#' perturbation(x, y)
#'
#' X <- rbind(
#'   c(1, 2, 3),
#'   c(4, 5, 6)
#' )
#' perturbation(X, c(1, 1, 2))
#' perturbation(c(1, 1, 2), X)
#'
#' @export
perturbation <- function(X, Y) {
  if (is.matrix(X)) {
    if (!is.numeric(X)) {
      stop("Composition must be numeric.", call. = FALSE)
    }

    Ymat <- perturbation_matrix(X, Y)
    P <- closure(X * Ymat)
    suppressWarnings(row.names(P) <- row.names(X))
    colnames(P) <- colnames(X)
    P

  } else if (is.atomic(X) && !is.list(X)) {
    if (!is.numeric(X)) {
      stop("Composition must be numeric.", call. = FALSE)
    }

    if (is.matrix(Y)) {
      if (!is.numeric(Y)) {
        stop("Composition must be numeric.", call. = FALSE)
      }
      Xmat <- perturbation_recycled_vector_matrix(X, Y)
      perturbation(Xmat, Y)

    } else if (inherits(Y, "data.frame")) {
      if (!all(sapply(Y, is.numeric))) {
        stop("All parts must be numeric.", call. = FALSE)
      }
      class_type <- class(Y)
      Ymat <- as.matrix(Y)
      Xmat <- perturbation_recycled_vector_matrix(X, Ymat)
      mP <- perturbation(Xmat, Ymat)
      P <- as.data.frame(mP)
      class(P) <- class_type
      suppressWarnings(row.names(P) <- row.names(Y))
      P

    } else {
      P <- Recall(matrix(X, nrow = 1), Y)
      P[1, ]
    }

  } else if (inherits(X, "data.frame")) {
    if (!all(sapply(X, is.numeric))) {
      stop("All parts must be numeric.", call. = FALSE)
    }
    class_type <- class(X)
    mP <- Recall(as.matrix(X), Y)
    P <- as.data.frame(mP)
    class(P) <- class_type
    suppressWarnings(row.names(P) <- row.names(X))
    P

  } else {
    stop("'X' must be a numeric vector, matrix or data.frame.", call. = FALSE)
  }
}

#' Recycle a perturbation vector to match a matrix operand
#'
#' @param x A numeric vector composition.
#' @param Y A numeric matrix used as the reference composition.
#' @return A numeric matrix with the same dimensions as \code{Y}.
#' @noRd
perturbation_recycled_vector_matrix <- function(x, Y) {
  if (length(x) != ncol(Y)) {
    stop("'X' must have length equal to the number of parts in 'Y'.", call. = FALSE)
  }

  Xmat <- matrix(x, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  row.names(Xmat) <- row.names(Y)
  colnames(Xmat) <- names(x)
  if (is.null(colnames(Xmat))) {
    colnames(Xmat) <- colnames(Y)
  }
  Xmat
}

#' Prepare a perturbation operand
#'
#' @param X A numeric matrix used as the reference composition.
#' @param Y A perturbation operand.
#' @return A numeric matrix with the same dimensions as \code{X}.
#' @noRd
perturbation_matrix <- function(X, Y) {
  if (inherits(Y, "data.frame")) {
    Y <- as.matrix(Y)
  }

  if (is.matrix(Y)) {
    if (!is.numeric(Y)) {
      stop("Composition must be numeric.", call. = FALSE)
    }

    if (all(dim(Y) == dim(X))) {
      Ymat <- Y
    } else if (nrow(Y) == 1 && ncol(Y) == ncol(X)) {
      Ymat <- matrix(rep(Y, each = nrow(X)), nrow = nrow(X), ncol = ncol(X), byrow = FALSE)
    } else {
      stop("'Y' must have the same dimension as 'X' or be a single composition vector.", call. = FALSE)
    }
  } else if (is.atomic(Y) && !is.list(Y)) {
    if (length(Y) == ncol(X)) {
      Ymat <- matrix(Y, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
    } else {
      stop("'Y' must have length equal to the number of parts in 'X'.", call. = FALSE)
    }
  } else {
    stop("'Y' must be a numeric vector, matrix or data.frame.", call. = FALSE)
  }

  colnames(Ymat) <- colnames(X)
  Ymat
}

#' Powering of compositional data
#'
#' The powering operation raises each part of a composition to a scalar exponent
#' and then applies closure to re-normalize the result as a composition.
#'
#' @param X A numeric vector, matrix or data.frame containing compositions.
#' @param alpha A numeric scalar or vector. If \code{X} is a matrix or
#'   data.frame, \code{alpha} may have length 1, length equal to the number of
#'   rows of \code{X}, or length equal to the number of parts of \code{X}. If
#'   \code{X} is a vector and \code{alpha} has length greater than 1, one
#'   powered composition is returned for each value of \code{alpha}.
#' @return An object with the same format as \code{X} containing the powered
#'   compositions, except that vector \code{X} with vector \code{alpha} returns
#'   a matrix.
#'
#' @details
#' Powering is the analogue of scalar multiplication in the simplex. Each part
#' is raised to \code{alpha}, and the result is closed with \code{\link{closure}}.
#' When \code{alpha} has one value per row, each composition is powered by its
#' corresponding value. When it has one value per part, each part receives its
#' corresponding exponent. For vector \code{X} and vector \code{alpha}, each row
#' of the result is \code{X} powered by the corresponding element of
#' \code{alpha}.
#'
#' @examples
#' x <- c(a = 1, b = 2, c = 3)
#' powering(x, 2)
#' powering(x, c(1, 2))
#'
#' X <- rbind(
#'   c(1, 2, 3),
#'   c(4, 5, 6)
#' )
#' powering(X, c(1, 2))
#'
#' @export
powering <- function(X, alpha) {
  if (is.matrix(X)) {
    if (!is.numeric(X)) {
      stop("Composition must be numeric.", call. = FALSE)
    }
    if (!is.numeric(alpha)) {
      stop("'alpha' must be numeric.", call. = FALSE)
    }

    P <- closure(X ^ powering_matrix(X, alpha))
    suppressWarnings(row.names(P) <- row.names(X))
    colnames(P) <- colnames(X)
    P

  } else if (is.atomic(X) && !is.list(X)) {
    if (!is.numeric(X)) {
      stop("Composition must be numeric.", call. = FALSE)
    }
    if (!is.numeric(alpha)) {
      stop("'alpha' must be numeric.", call. = FALSE)
    }

    if (length(alpha) == 1L) {
      P <- Recall(matrix(X, nrow = 1), alpha)
      P[1, ]
    } else {
      Xmat <- matrix(X, nrow = length(alpha), ncol = length(X), byrow = TRUE)
      colnames(Xmat) <- names(X)
      Recall(Xmat, alpha)
    }

  } else if (inherits(X, "data.frame")) {
    if (!all(sapply(X, is.numeric))) {
      stop("All parts must be numeric.", call. = FALSE)
    }
    class_type <- class(X)
    mP <- Recall(as.matrix(X), alpha)
    P <- as.data.frame(mP)
    class(P) <- class_type
    suppressWarnings(row.names(P) <- row.names(X))
    P

  } else {
    stop("'X' must be a numeric vector, matrix or data.frame.", call. = FALSE)
  }
}

#' Prepare a powering exponent matrix
#'
#' @param X A numeric matrix used as the reference composition.
#' @param alpha A numeric scalar or vector of exponents.
#' @return A numeric matrix with the same dimensions as \code{X}.
#' @noRd
powering_matrix <- function(X, alpha) {
  if (length(alpha) == 1L) {
    matrix(alpha, nrow = nrow(X), ncol = ncol(X))
  } else if (length(alpha) == nrow(X)) {
    matrix(alpha, nrow = nrow(X), ncol = ncol(X), byrow = FALSE)
  } else if (length(alpha) == ncol(X)) {
    matrix(alpha, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  } else {
    stop("'alpha' must be a scalar or a numeric vector with length equal to the number of rows or parts of 'X'.", call. = FALSE)
  }
}
