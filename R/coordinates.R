#' Coordinates of compositions with respect to a basis
#'
#' @description
#' Compute coordinates of a composition or a compositional data set with respect
#' to a given log-ratio basis.
#'
#' The `basis` argument can be either:
#' \itemize{
#'   \item a character string identifying a predefined coordinate system, or
#'   \item a matrix whose columns define a system of log-contrasts.
#' }
#'
#' The predefined options are:
#' \itemize{
#'   \item `"ilr"`: isometric log-ratio coordinates,
#'   \item `"olr"`: orthonormal log-ratio coordinates,
#'   \item `"clr"`: centered log-ratio coordinates,
#'   \item `"alr"`: additive log-ratio coordinates,
#'   \item `"pw"`: pairwise log-ratios,
#'   \item `"pc"`: principal component log-ratio coordinates,
#'   \item `"pb"`: principal balance coordinates,
#'   \item `"cdp"`: CoDaPack default balances.
#' }
#'
#' @param X A compositional data set. It can be a numeric matrix, a data frame,
#'   or a numeric vector.
#' @param basis Basis used to compute the coordinates. Either a character string
#'   naming a predefined basis or a matrix with log-ratio basis vectors in
#'   columns.
#' @param ... components of the composition
#'
#' @return
#' Coordinates of `X` with respect to the given `basis`. The returned object has
#' the same general type as the input when possible.
#'
#' @seealso
#' \code{\link{ilr_basis}}, \code{\link{alr_basis}}, \code{\link{clr_basis}},
#' \code{\link{sbp_basis}}, \code{\link{composition}}
#'
#' @examples
#' coordinates(1:5)
#'
#' B <- ilr_basis(5)
#' coordinates(1:5, B)
#'
#' X <- rbind(1:5, 2:6)
#' coordinates(X, "clr")
#'
#' @export
coordinates <- function(X, basis = "ilr") {
  if (is.matrix(X)) {
    if (!is.numeric(X)) {
      stop("Composition must be numeric.", call. = FALSE)
    }

    if (is.character(basis)) {
      basis <- match.arg(
        basis,
        c("ilr", "olr", "clr", "alr", "pw", "pc", "pb", "cdp")
      )

      COORD <- switch(
        basis,
        ilr = ilr(X),
        olr = olr(X),
        clr = clr(X),
        alr = alr(X),
        pw  = pw(X),
        pc  = pc(X),
        pb  = pb(X),
        cdp = cdp(X)
      )
    } else if (is.matrix(basis) || inherits(basis, "sparseMatrix")) {
      COORD <- coordinates.matrix(X, basis)
      colnames(COORD) <- colnames(basis)
      if (is.null(colnames(COORD))) {
        colnames(COORD) <- sprintf("h%d", seq_len(ncol(COORD)))
      }
    } else {
      stop("'basis' must be either a character string or a matrix.", call. = FALSE)
    }

  } else if (is.atomic(X) && !is.list(X)) {
    COORD <- Recall(matrix(X, nrow = 1), basis)
    COORD <- COORD[1, ]

  } else if (inherits(X, "data.frame")) {
    if (!all(sapply(X, is.numeric))) {
      stop("All parts must be numeric.", call. = FALSE)
    }
    class_type <- class(X)
    mCOORD <- Recall(as.matrix(X), basis)
    COORD <- as.data.frame(mCOORD)
    class(COORD) <- class_type

  } else {
    stop("'X' must be a numeric vector, matrix or data.frame.", call. = FALSE)
  }

  suppressWarnings(row.names(COORD) <- row.names(X))
  COORD
}

#' Compositions from coordinates with respect to a basis
#'
#' @description
#' Reconstruct a composition from coordinates with respect to a given basis.
#'
#' @param H Coordinates of a composition. It can be a numeric matrix, a data
#'   frame, or a numeric vector.
#' @param basis Basis used to interpret the coordinates. Either a character
#'   string naming a predefined basis or a matrix.
#'
#' @return
#' A composition corresponding to the given coordinates.
#'
#' @seealso
#' \code{\link{coordinates}}, \code{\link{ilr_basis}}, \code{\link{alr_basis}},
#' \code{\link{clr_basis}}, \code{\link{sbp_basis}}
#'
#' @export
composition <- function(H, basis = "ilr") {
  if (is.matrix(H)) {
    if (!is.numeric(H)) {
      stop("Coordinates must be numeric.", call. = FALSE)
    }

    if (is.character(basis)) {
      basis <- match.arg(basis, c("ilr", "olr", "alr", "clr"))

      if (basis == "ilr" || basis == "olr") {
        D <- ncol(H) + 1
        P <- exp(H %*% t(ilr_basis(D)))

      } else if (basis == "alr") {
        P <- cbind(exp(H), 1)

      } else if (basis == "clr") {
        P <- exp(H)
      }
      colnames(P) <- paste0('c', 1:ncol(P))
    } else if (is.matrix(basis) || inherits(basis, "sparseMatrix")) {
      P <- exp(as.matrix(H) %*% pinv(as.matrix(basis)))

      if (!is.null(rownames(basis))) {
        colnames(P) <- rownames(basis)
      }

    } else {
      stop("'basis' must be either a character string or a matrix.", call. = FALSE)
    }

    P <- P / rowSums(P)

  } else if (is.atomic(H) && !is.list(H)) {
    P <- Recall(matrix(H, nrow = 1), basis)[1, ]

  } else if (inherits(H, "data.frame")) {
    if (!all(sapply(H, is.numeric))) {
      stop("All coordinates must be numeric.", call. = FALSE)
    }

    class_type <- class(H)
    P <- as.data.frame(Recall(as.matrix(H), basis))
    class(P) <- class_type

  } else {
    stop("'H' must be a numeric vector, matrix or data.frame.", call. = FALSE)
  }

  suppressWarnings(row.names(P) <- row.names(H))
  P
}

#' @rdname coordinates
#' @export
coord <- function(..., basis = "ilr") {
  largs <- list(...)
  lpars <- as.list(substitute(largs))

  inum <- sapply(
    seq_along(lpars),
    function(i) {
      is.numeric(lpars[[i]]) &&
        !is.matrix(lpars[[i]]) &&
        length(lpars[[i]]) > 1
    }
  )

  if (sum(inum) == 0) {
    if (length(largs) == 1 &&
        (is.matrix(largs[[1]]) || is.data.frame(largs[[1]]) || is.vector(largs[[1]]))) {
      coordinates(largs[[1]], basis = basis)
    } else {
      if (length(largs) >= 2 && is.character(largs[[2]])) {
        coordinates(largs[[1]], basis = largs[[2]])
      } else {
        stop("Please specify the basis as the second argument.", call. = FALSE)
      }
    }
  } else {
    if (sum(inum) == 1) {
      stop(
        "Compositions should have at least two parts. If you want to calculate ",
        "the coordinates of one composition, please use function coordinates().",
        call. = FALSE
      )
    }
    if (sum(inum) > 1 && sum(inum) < length(lpars)) {
      stop("All components should be numeric.", call. = FALSE)
    }
    coordinates(cbind(...), basis = basis)
  }
}

#' @rdname coordinates
#' @export
alr_c <- function(X) {
  coordinates(X, "alr")
}

#' @rdname coordinates
#' @export
clr_c <- function(X) {
  coordinates(X, "clr")
}

#' @rdname coordinates
#' @export
ilr_c <- function(X) {
  coordinates(X, "ilr")
}

#' @rdname coordinates
#' @export
olr_c <- function(X) {
  coordinates(X, "olr")
}

#' @rdname composition
#' @export
comp <- composition

alr <- function(X) {
  COORD <- alr_coordinates(X, ncol(X))
  colnames(COORD) <- paste0("alr", seq_len(ncol(COORD)))
  COORD
}

ilr <- function(X) {
  COORD <- coordinates(X, ilr_basis(ncol(X)))
  colnames(COORD) <- paste0("ilr", seq_len(ncol(COORD)))
  COORD
}

olr <- function(X) {
  COORD <- coordinates(X, ilr_basis(ncol(X)))
  colnames(COORD) <- paste0("olr", seq_len(ncol(COORD)))
  COORD
}

clr <- function(X) {
  COORD <- clr_coordinates(X)
  colnames(COORD) <- paste0("clr", seq_len(ncol(COORD)))
  COORD
}

pc <- function(X) {
  B <- ilr_basis(ncol(X))
  B <- B %*% svd(scale(log(X) %*% B, scale = FALSE))$v
  COORD <- matrix_coordinates(X, B)
  colnames(COORD) <- paste0("pc", seq_len(ncol(B)))
  COORD
}

cdp <- function(X) {
  B <- cdp_basis(ncol(X))
  COORD <- matrix_coordinates(X, B)
  colnames(COORD) <- colnames(B)
  COORD
}

pb <- function(X) {
  B <- pb_basis(X, method = "exact")
  COORD <- matrix_coordinates(X, B)
  colnames(COORD) <- colnames(B)
  COORD
}

pw <- function(X) {
  B <- pairwise_basis(ncol(X))
  COORD <- matrix_coordinates(X, B)
  colnames(COORD) <- colnames(B)
  COORD
}

coordinates.matrix <- function(X, basis) {
  if (inherits(basis, "sparseMatrix")) {
    sparse_coordinates(X, basis)
  } else {
    matrix_coordinates(X, basis)
  }
}
