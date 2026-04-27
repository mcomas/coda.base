get_part_names <- function(x) {
  if (is.vector(x)) return(names(x))
  colnames(x)
}

#' Isometric and orthonormal log-ratio bases
#'
#' Construct an isometric log-ratio (ilr) basis for a composition with
#' \eqn{D} parts. The ilr basis is an orthonormal basis of the clr-plane and
#' provides \eqn{D - 1} coordinates. The same basis is sometimes referred to as
#' an orthonormal log-ratio (olr) basis.
#'
#' For `type = "default"`, the function returns the standard Helmert-type ilr
#' basis. Alternative constructions are available through `type = "pivot"` and
#' `type = "cdp"`.
#'
#' The default basis vectors are:
#' \deqn{
#' h_i = \sqrt{\frac{i}{i+1}}
#' \log \frac{\sqrt[i]{\prod_{j=1}^i x_j}}{x_{i+1}},
#' \qquad i = 1, \ldots, D - 1
#' }{
#' h[i] = sqrt(i / (i + 1)) * (log(x[1] * ... * x[i]) / i - log(x[i + 1]))
#' }
#'
#' @param dim Number of parts. It can be:
#'   \itemize{
#'   \item a single integer,
#'   \item a matrix or data frame, in which case the number of columns is used,
#'   \item a character vector of part names, in which case its length is used.
#'   }
#' @param type Type of ilr basis to construct. Available options are:
#'   \itemize{
#'   \item `"default"`: standard Helmert-type ilr basis,
#'   \item `"pivot"`: pivot balance basis,
#'   \item `"cdp"`: CoDaPack default basis.
#'   }
#'
#' @return A matrix with \eqn{D} rows and \eqn{D - 1} columns representing an
#'   orthonormal log-ratio basis.
#'
#' @references
#' Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal, C.
#' (2003). \emph{Isometric logratio transformations for compositional data
#' analysis}. Mathematical Geology, \strong{35}(3), 279--300.
#'
#' @examples
#' ilr_basis(5)
#' ilr_basis(alimentation[, 1:9])
#' ilr_basis(c("a", "b", "c", "d"), type = "pivot")
#'
#' @export
ilr_basis <- function(dim, type = "default") {
  parts <- get_part_names(dim)

  if (is.character(dim)) {
    parts <- dim
    dim <- length(dim)
  }

  D <- check_dim(dim)

  if (type == "cdp") {
    B <- cdp_basis_(D)
  } else {
    B <- ilr_basis_default(D)
    if (type == "pivot") {
      B <- (-B)[, ncol(B):1, drop = FALSE][nrow(B):1, ]
    }
  }

  rownames(B) <- sprintf("c%d", seq_len(nrow(B)))
  if (!is.null(parts)) rownames(B) <- parts
  colnames(B) <- sprintf("ilr%d", seq_len(ncol(B)))
  B
}

#' @rdname ilr_basis
#' @export
olr_basis <- function(dim, type = "default") {
  B <- ilr_basis(dim, type)
  colnames(B) <- sprintf("olr%d", seq_len(ncol(B)))
  B
}

#' CoDaPack default ilr basis
#'
#' Construct the default isometric log-ratio basis used in CoDaPack.
#'
#' @param dim Number of parts. It can be a single integer, a matrix or data
#'   frame, or a character vector of part names.
#'
#' @return A matrix with \eqn{D} rows and \eqn{D - 1} columns containing the
#'   CoDaPack default ilr basis.
#'
#' @examples
#' cdp_basis(5)
#' cdp_basis(c("a", "b", "c", "d"))
#'
#' @export
cdp_basis <- function(dim) {
  parts <- get_part_names(dim)

  if (is.character(dim)) {
    parts <- dim
    dim <- length(dim)
  }

  D <- check_dim(dim)
  B <- cdp_basis_(D)

  rownames(B) <- paste0("c", seq_len(D))
  if (!is.null(parts)) rownames(B) <- parts
  colnames(B) <- paste0("ilr", seq_len(ncol(B)))
  B
}

#' Centered log-ratio basis
#'
#' Construct the transformation matrix associated with centered log-ratio (clr)
#' coordinates.
#'
#' CLR coordinates are linearly dependent and lie in the \eqn{D - 1}
#' dimensional clr-plane.
#'
#' @param dim Number of parts. It can be a single integer, a matrix or data
#'   frame, or a character vector of part names.
#'
#' @return A square matrix defining the clr coordinate system.
#'
#' @references
#' Aitchison, J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman & Hall, London.
#'
#' @examples
#' B <- clr_basis(5)
#' clr_coordinates <- coordinates(c(1, 2, 3, 4, 5), B)
#' sum(clr_coordinates) < 1e-15
#'
#' @export
clr_basis <- function(dim) {
  parts <- get_part_names(dim)

  if (is.character(dim)) {
    parts <- dim
    dim <- length(dim)
  }

  D <- check_dim(dim)
  B <- clr_basis_default(D)

  colnames(B) <- sprintf("clr%d", seq_len(ncol(B)))
  rownames(B) <- sprintf("c%d", seq_len(nrow(B)))
  if (!is.null(parts)) rownames(B) <- parts
  B
}

#' Additive log-ratio basis
#'
#' Construct the transformation matrix associated with additive log-ratio (alr)
#' coordinates.
#'
#' @param dim Number of parts. It can be a single integer, a matrix or data
#'   frame, or a character vector of part names.
#' @param denominator Part used as denominator. By default, the last part is
#'   used.
#' @param numerator Parts used as numerators. By default, all parts except the
#'   denominator are used, preserving their original order.
#'
#' @return A matrix defining the alr coordinate system.
#'
#' @references
#' Aitchison, J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman & Hall, London.
#'
#' @examples
#' alr_basis(5)
#' alr_basis(5, 3)
#' alr_basis(5, 3, c(1, 5, 2, 4))
#'
#' @export
alr_basis <- function(dim, denominator = NULL, numerator = NULL) {
  parts <- get_part_names(dim)

  if (is.character(dim)) {
    parts <- dim
    dim <- length(dim)
  }

  D <- check_dim(dim)

  if (is.null(denominator)) denominator <- D
  if (is.null(numerator)) numerator <- which(seq_len(D) != denominator)

  res <- alr_basis_default(D)
  res <- cbind(res, 0)

  if (D != denominator) {
    res[c(denominator, D), ] <- res[c(D, denominator), , drop = FALSE]
    res[, c(denominator, D)] <- res[, c(D, denominator), drop = FALSE]
  }

  B <- res[, numerator, drop = FALSE]

  colnames(B) <- sprintf("alr%d", seq_len(ncol(B)))
  rownames(B) <- sprintf("c%d", seq_len(nrow(B)))
  if (!is.null(parts)) rownames(B) <- parts
  B
}

#' Principal component log-ratio basis
#'
#' Construct an ilr basis rotated according to the principal components of the
#' log-ratio coordinates of a compositional data set.
#'
#' @param X Compositional data set.
#'
#' @return A matrix whose columns define a principal-component-oriented ilr
#'   basis.
#'
#' @export
pc_basis <- function(X) {
  B <- ilr_basis(ncol(X))
  Xilr <- coordinates(X, B)
  B <- B %*% svd(scale(Xilr, scale = FALSE))$v

  parts <- get_part_names(X)
  if (is.null(parts)) {
    parts <- paste0("c", seq_len(nrow(B)))
  }

  rownames(B) <- parts
  colnames(B) <- paste0("pc", seq_len(ncol(B)))
  as.matrix(B)
}

#' Canonical-correlation log-ratio basis
#'
#' Construct an ilr basis rotated according to canonical correlations between a
#' compositional response data set and an explanatory data set.
#'
#' @param Y A compositional data set.
#' @param X An explanatory data set.
#'
#' @return A matrix whose columns define a canonical-correlation-oriented ilr
#'   basis.
#'
#' @export
cc_basis <- function(Y, X) {
  Y <- as.matrix(Y)
  X <- cbind(X)

  B <- ilr_basis(ncol(Y))
  cc <- stats::cancor(coordinates(Y), X)
  B <- B %*% cc$xcoef

  parts <- colnames(Y)
  if (is.null(parts)) {
    parts <- paste0("c", seq_len(nrow(B)))
  }

  rownames(B) <- parts
  colnames(B) <- paste0("cc", seq_len(ncol(B)))
  B
}

#' Basis from a sequential binary partition
#'
#' Construct a balance basis from a sequential binary partition (SBP) or from a
#' more general collection of balances.
#'
#' The argument `sbp` can be specified in two ways:
#' \itemize{
#' \item as a list of formulas, where each formula defines the numerator and the
#'   denominator groups of a balance,
#' \item as a matrix with one column per balance and one row per part. Positive
#'   entries indicate parts in the numerator, negative entries indicate parts in
#'   the denominator, and zeros indicate unused parts.
#' }
#'
#' @param sbp A list of formulas or a matrix describing balances.
#' @param data Optional compositional data set used to extract part names when
#'   `sbp` is given as a list of formulas.
#' @param fill Logical; if `TRUE`, complete the supplied balances to obtain a
#'   full basis.
#' @param silent Logical; if `FALSE`, report whether the resulting balances form
#'   a basis, and whether they are orthogonal or orthonormal.
#'
#' @return A matrix whose columns are balances.
#'
#' @examples
#' X <- data.frame(
#'   a = 1:2, b = 2:3, c = 4:5,
#'   d = 5:6, e = 10:11, f = 100:101, g = 1:2
#' )
#'
#' # Sequential SBP construction
#' sbp_basis(list(
#'   b1 = a ~ b + c + d + e + f + g,
#'   b2 = b ~ c + d + e + f + g,
#'   b3 = c ~ d + e + f + g,
#'   b4 = d ~ e + f + g,
#'   b5 = e ~ f + g,
#'   b6 = f ~ g
#' ), data = X)
#'
#' # Chain construction
#' sbp_basis(list(
#'   b1 = a ~ b,
#'   b2 = b1 ~ c,
#'   b3 = b2 ~ d,
#'   b4 = b3 ~ e,
#'   b5 = b4 ~ f,
#'   b6 = b5 ~ g
#' ), data = X)
#'
#' # Non-orthogonal system of balances
#' sbp_basis(list(
#'   b1 = a + b + c ~ e + f + g,
#'   b2 = d ~ a + b + c,
#'   b3 = d ~ e + g,
#'   b4 = a ~ e + b,
#'   b5 = b ~ f,
#'   b6 = c ~ g
#' ), data = X)
#'
#' # Direct construction from a contrast matrix
#' sbp_basis(cbind(
#'   c( 1,  1, -1, -1),
#'   c( 1, -1,  1, -1),
#'   c( 1, -1, -1,  1)
#' ))
#'
#' @export
sbp_basis <- function(sbp, data = NULL, fill = FALSE, silent = FALSE) {
  if (is.matrix(sbp)) {
    RES <- apply(sign(sbp), 2, function(x) {
      s1 <- x == 1
      s2 <- x == -1
      I1 <- sum(s1)
      I2 <- sum(s2)

      if (I1 == 0 || I2 == 0) {
        stop("Each balance must contain at least one positive and one negative part.")
      }

      l <- +1 / I1 * sqrt(I1 * I2 / (I1 + I2))
      r <- -1 / I2 * sqrt(I1 * I2 / (I1 + I2))

      bal <- rep(0, length(x))
      bal[s1] <- l
      bal[s2] <- r
      bal
    })

    if (is.null(rownames(sbp))) {
      rownames(RES) <- paste0("c", seq_len(nrow(sbp)))
    } else {
      rownames(RES) <- rownames(sbp)
    }

    if (is.null(colnames(sbp))) {
      colnames(RES) <- paste0("b", seq_len(ncol(sbp)))
    } else {
      colnames(RES) <- colnames(sbp)
    }

  } else {
    if (is.character(data)) {
      part_names <- data
    } else {
      if (!is.data.frame(data) &&
          !is.environment(data) &&
          ((is.matrix(data) && !is.null(colnames(data))) ||
           !is.null(attr(data, "class")))) {
        data <- as.data.frame(data)
      } else if (is.array(data)) {
        stop("'data' must be a data.frame or a matrix with column names")
      }
      part_names <- colnames(data)
    }

    if (!all(unlist(lapply(sbp, all.vars)) %in% c(part_names, names(sbp)))) {
      stop("Balances should be columns of 'data'")
    }

    nms <- setdiff(names(sbp), "")

    if (length(nms) > 0) {
      substitutions <- lapply(sbp, all.vars)
      substitutions <- substitutions[nms]

      while (!all(is.na(substitutions)) &&
             !all(unlist(substitutions) %in% part_names)) {
        for (nm in nms) {
          substitutions <- lapply(substitutions, function(subs) {
            I <- match(nm, subs)
            if (!is.na(I)) {
              c(subs[-I], substitutions[[nm]])
            } else {
              subs
            }
          })
        }
      }
    }

    sbp_split <- function(part) {
      RIGHT <- attr(stats::terms(part), "term.labels")
      LEFT <- setdiff(all.vars(part), RIGHT)

      if (length(nms) > 0) {
        for (nm in nms) {
          I <- match(nm, RIGHT)
          if (!is.na(I)) {
            RIGHT <- c(RIGHT[-I], substitutions[[nm]])
          }
          I <- match(nm, LEFT)
          if (!is.na(I)) {
            LEFT <- c(LEFT[-I], substitutions[[nm]])
          }
        }
      }

      list(LEFT, RIGHT)
    }

    sbp_clean <- lapply(sbp, sbp_split)

    RES <- sapply(sbp_clean, function(balance) {
      I1 <- length(balance[[1]])
      I2 <- length(balance[[2]])

      l <- +1 / I1 * sqrt(I1 * I2 / (I1 + I2))
      r <- -1 / I2 * sqrt(I1 * I2 / (I1 + I2))

      bal <- stats::setNames(rep(0, length(part_names)), part_names)
      bal[balance[[1]]] <- l
      bal[balance[[2]]] <- r
      bal
    })

    rownames(RES) <- part_names
    colnames(RES) <- names(sbp)
  }

  if (fill) {
    RES_F <- Recall(fill_sbp(sign(RES)))
    rownames(RES_F) <- rownames(RES)
    colnames(RES_F) <- paste0(".b", seq_len(ncol(RES_F)))
    colnames(RES_F)[seq_len(ncol(RES))] <- colnames(RES)
    return(RES_F)
  }

  if (!silent) {
    if (qr(RES)$rank != (nrow(RES) - 1)) {
      warning("Given partition is not a basis")
    } else {
      Z <- t(RES) %*% RES
      tol <- 1e-10

      off <- Z
      diag(off) <- 0

      if (any(abs(off) > tol)) {
        warning("Given basis is not orthogonal")
      } else if (any(abs(diag(Z) - 1) > tol)) {
        warning("Given basis is not orthonormal")
      }
    }
  }

  RES
}

#' Principal balance basis
#'
#' Construct a basis of principal balances for a compositional data set.
#'
#' Several methods are available:
#' \itemize{
#' \item `"exact"`: exact computation of principal balances,
#' \item `"exact2"`: exact computation using incremental Gray-code updates,
#' \item `"constrained"`: constrained approximation based on a target criterion,
#' \item `"cluster"`: approximation based on hierarchical clustering.
#' }
#'
#' @param X Compositional data set.
#' @param method Method used to construct the principal balances. One of
#'   `"exact"`, `"exact2"`, `"constrained"`, or `"cluster"`.
#' @param constrained.criterion Criterion used by the constrained method.
#'   Either `"variance"` (default) or `"angle"`.
#' @param cluster.method Linkage criterion passed to
#'   \code{\link[stats]{hclust}} when `method = "cluster"`.
#' @param ordering Logical; if `TRUE`, reorder balances by decreasing explained
#'   variance.
#' @param ... Additional arguments passed to \code{\link[stats]{hclust}}.
#'
#' @return A matrix whose columns are principal balances.
#'
#' @references
#' Martín-Fernández, J. A., Pawlowsky-Glahn, V., Egozcue, J. J., &
#' Tolosana-Delgado, R. (2018). \emph{Advances in Principal Balances for
#' Compositional Data}. Mathematical Geosciences, 50, 273--298.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(exp(rnorm(5 * 100)), nrow = 100, ncol = 5)
#'
#' v1 <- apply(coordinates(X, "pc"), 2, var)
#' v2 <- apply(coordinates(X, pb_basis(X, method = "exact")), 2, var)
#' v3 <- apply(coordinates(X, pb_basis(X, method = "constrained")), 2, var)
#' v4 <- apply(coordinates(X, pb_basis(X, method = "cluster")), 2, var)
#'
#' barplot(
#'   rbind(v1, v2, v3, v4),
#'   beside = TRUE,
#'   ylim = c(0, 2),
#'   legend = c(
#'     "Principal Components",
#'     "PB (Exact method)",
#'     "PB (Constrained)",
#'     "PB (Ward approximation)"
#'   ),
#'   names = paste0("Comp.", 1:4),
#'   args.legend = list(cex = 0.8),
#'   ylab = "Variance"
#' )
#'
#' @export
pb_basis <- function(X, method,
                     constrained.criterion = "variance",
                     cluster.method = "ward.D2",
                     ordering = TRUE, ...) {
  X <- as.matrix(X)

  if (!all(X > 0)) {
    stop("All components must be strictly positive.", call. = FALSE)
  }

  if (method %in% c("constrained", "exact", "exact2")) {
    if (method == "exact") {
      B <- find_PB(X)
    }
    if (method == "exact2") {
      B <- find_PB2(X)
    }

    if (method == "constrained") {
      if (constrained.criterion %in% c("angle", "variance")) {
        if (constrained.criterion == "angle") {
          B <- constrained_pb(X, angle = TRUE)
        }
        if (constrained.criterion == "variance") {
          B <- constrained_pb(X, angle = FALSE)
        }
      } else {
        stop("Error: 'constrained.criterion' is not available")
      }
    }

  } else if (method == "cluster") {
    hh <- stats::hclust(
      stats::as.dist(variation_array(X)),
      method = cluster.method,
      ...
    )

    B <- matrix(0, ncol = nrow(hh$merge), nrow = ncol(X))

    for (i in seq_len(nrow(hh$merge))) {
      if (hh$merge[i, 1] < 0 && hh$merge[i, 2] < 0) {
        B[-hh$merge[i, ], i] <- c(-1, +1)
      } else {
        if (hh$merge[i, 1] > 0) {
          B[B[, hh$merge[i, 1]] != 0, i] <- -1
        } else {
          B[-hh$merge[i, 1], i] <- -1
        }

        if (hh$merge[i, 2] > 0) {
          B[B[, hh$merge[i, 2]] != 0, i] <- +1
        } else {
          B[-hh$merge[i, 2], i] <- +1
        }
      }
    }

    B <- sbp_basis(B[, nrow(hh$merge):1, drop = FALSE])

  } else {
    stop(sprintf("Method '%s' does not exist", method))
  }

  if (ordering) {
    B <- B[, order(apply(coordinates(X, B), 2, stats::var), decreasing = TRUE), drop = FALSE]
  }

  parts <- colnames(X)
  if (is.null(parts)) {
    parts <- paste0("c", seq_len(nrow(B)))
  }

  rownames(B) <- parts
  colnames(B) <- paste0("pb", seq_len(ncol(B)))
  B
}

#' Pairwise log-ratio generating system
#'
#' Construct the system of all pairwise log-ratios between parts.
#'
#' @param dim Number of parts. It can be a single integer, a matrix or data
#'   frame, or a character vector of part names.
#'
#' @return A matrix, or a sparse matrix for large dimensions, whose columns
#'   represent all pairwise log-ratio generators.
#'
#' @export
pairwise_basis <- function(dim) {
  parts <- get_part_names(dim)

  if (is.character(dim)) {
    parts <- dim
    dim <- length(dim)
  }

  D <- check_dim(dim)
  I <- utils::combn(D, 2)

  if (D > 100) {
    B <- Matrix::sparseMatrix(
      i = c(I[1, ], I[2, ]),
      j = c(seq_len(ncol(I)), seq_len(ncol(I))),
      x = rep(c(1, -1), each = ncol(I))
    )
  } else {
    B <- apply(I, 2, function(i) {
      b <- rep(0, D)
      b[i] <- c(1, -1)
      b
    })
  }

  colnames(B) <- paste0("pw", apply(I, 2, paste, collapse = "_"))
  rownames(B) <- paste0("c", seq_len(D))
  if (!is.null(parts)) rownames(B) <- parts
  B
}

cdp_basis_ <- function(dim,
                       wR = seq_len(ceiling(dim / 2)),
                       wL = seq(ceiling(dim / 2) + 1, dim)) {
  R <- length(wR)
  L <- length(wL)
  D <- R + L

  v <- rep(0, dim)
  v[wR] <- +sqrt(L / R / D)
  v[wL] <- -sqrt(R / L / D)

  if (R == 1 && L == 1) {
    return(v)
  }

  if (R == 1) {
    return(cbind(
      v,
      Recall(
        dim,
        wR = wL[1:ceiling(L / 2)],
        wL = wL[ceiling(L / 2) + 1:floor(L / 2)]
      )
    ))
  }

  if (L == 1) {
    return(cbind(
      v,
      Recall(
        dim,
        wR = wR[1:ceiling(R / 2)],
        wL = wR[ceiling(R / 2) + 1:floor(R / 2)]
      )
    ))
  }

  cbind(
    v,
    Recall(
      dim,
      wR = wR[1:ceiling(R / 2)],
      wL = wR[ceiling(R / 2) + 1:floor(R / 2)]
    ),
    Recall(
      dim,
      wR = wL[1:ceiling(L / 2)],
      wL = wL[ceiling(L / 2) + 1:floor(L / 2)]
    )
  )
}

check_dim <- function(dim) {
  if (!is.null(ncol(dim))) {
    dim <- ncol(dim)
  }

  if (is.vector(dim) && length(dim) > 1) {
    dim <- length(dim)
  }

  if (!is.numeric(dim) || length(dim) != 1 || is.na(dim)) {
    stop("Dimension should be a single numeric value.", call. = FALSE)
  }

  if (dim <= 1) {
    stop("Dimension should be at least 2.", call. = FALSE)
  }

  as.integer(dim)
}
