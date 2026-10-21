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

# Parse only part names, additions and parentheses; model formula operators do
# not describe groups of compositional parts.
.sbp_group <- function(expr) {
  pending <- list(expr)
  parts <- character()
  while (length(pending)) {
    last <- length(pending)
    node <- pending[[last]]
    pending[[last]] <- NULL
    if (is.symbol(node) && nzchar(as.character(node))) {
      parts <- c(parts, as.character(node))
    } else if (is.call(node) && identical(node[[1]], as.name("+")) &&
               length(node) == 3L) {
      pending <- c(pending, list(node[[3]], node[[2]]))
    } else if (is.call(node) && identical(node[[1]], as.name("(")) &&
               length(node) == 2L) {
      pending <- c(pending, list(node[[2]]))
    } else {
      stop("Balance groups must contain part names joined by '+'.")
    }
  }
  unique(parts)
}

.sbp_formula_matrix <- function(sbp, data) {
  if (!is.list(sbp) || !length(sbp) ||
      !all(vapply(sbp, inherits, logical(1), "formula"))) {
    stop("'sbp' must be a numeric matrix or a non-empty list of formulas.")
  }
  K <- length(sbp)
  balance_names <- names(sbp)
  if (is.null(balance_names)) balance_names <- rep("", K)
  if (anyNA(balance_names) || anyDuplicated(balance_names[nzchar(balance_names)])) {
    stop("Balance names must be non-missing and unique when supplied.")
  }
  named <- nzchar(balance_names)
  groups <- lapply(sbp, function(formula) {
    if (length(formula) != 3L) {
      stop("Each balance formula must have a numerator and a denominator.")
    }
    list(.sbp_group(formula[[2]]), .sbp_group(formula[[3]]))
  })
  variables <- lapply(groups, function(g) unique(c(g[[1]], g[[2]])))
  references <- balance_names[named]
  if (is.null(data)) {
    part_names <- setdiff(unique(unlist(variables, use.names = FALSE)), references)
  } else if (is.character(data) && is.null(dim(data))) {
    part_names <- data
  } else if (is.matrix(data) || is.data.frame(data)) {
    part_names <- colnames(data)
  } else {
    stop("'data' must be a character vector, a data.frame or a matrix with column names.")
  }
  if (is.null(part_names) || anyNA(part_names) || any(!nzchar(part_names)) ||
      anyDuplicated(part_names)) {
    stop("Part names must be non-missing, non-empty and unique.")
  }
  if (any(references %in% part_names)) {
    stop("Balance names must not also be part names.")
  }
  unknown <- setdiff(unlist(variables, use.names = FALSE), c(part_names, references))
  if (length(unknown)) {
    stop("Unknown parts or balance references: ", paste(unknown, collapse = ", "))
  }

  # Resolve the dependency graph once, without recursion, so long chains and
  # forward references work and cycles fail before expanding any cyclic group.
  dependencies <- lapply(variables, function(vars) {
    indices <- match(vars, balance_names, nomatch = 0L)
    indices[indices != 0L]
  })
  remaining <- lengths(dependencies)
  dependents <- vector("list", K)
  for (i in seq_len(K)) {
    for (j in dependencies[[i]]) dependents[[j]] <- c(dependents[[j]], i)
  }
  queue <- integer(K)
  ready <- which(remaining == 0L)
  queue_end <- length(ready)
  queue[seq_len(queue_end)] <- ready
  head <- 1L
  while (head <= queue_end) {
    i <- queue[head]
    head <- head + 1L
    next_indices <- dependents[[i]]
    remaining[next_indices] <- remaining[next_indices] - 1L
    ready <- next_indices[remaining[next_indices] == 0L]
    if (length(ready)) {
      queue[queue_end + seq_along(ready)] <- ready
      queue_end <- queue_end + length(ready)
    }
  }
  if (queue_end != K) {
    stop("Circular balance references; unresolved balances: ",
         paste(balance_names[remaining > 0L], collapse = ", "))
  }
  if (length(part_names) < 2L) stop("At least two parts are required.")

  expanded <- vector("list", K)
  signs <- matrix(0, length(part_names), K)
  expand <- function(tokens) {
    refs <- match(tokens, balance_names, nomatch = 0L)
    unique(unlist(lapply(seq_along(tokens), function(j) {
      if (refs[j] != 0L) expanded[[refs[j]]] else match(tokens[j], part_names)
    }), use.names = FALSE))
  }
  for (i in queue) {
    numerator <- expand(groups[[i]][[1]])
    denominator <- expand(groups[[i]][[2]])
    if (!length(numerator) || !length(denominator) ||
        any(numerator %in% denominator)) {
      stop("Each balance must have non-empty, disjoint numerator and denominator groups (balance ",
           i, ").")
    }
    signs[numerator, i] <- 1
    signs[denominator, i] <- -1
    expanded[[i]] <- c(numerator, denominator)
  }
  # Keep supplied names intact, even when they coincide with generated names.
  if (any(!named)) {
    generated <- make.unique(c(references, paste0("b", which(!named))))
    balance_names[!named] <- generated[length(references) + seq_len(sum(!named))]
  }
  dimnames(signs) <- list(part_names, balance_names)
  signs
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
#' Formula groups contain part names joined by `+`, optionally enclosed in
#' parentheses. Repeated parts within a group are counted once; numerator and
#' denominator groups must be non-empty and disjoint. A named balance can be
#' referenced in another formula to include all its parts, regardless of their
#' signs. Forward references are allowed, but circular references and ambiguous
#' names shared by balances and parts are rejected.
#'
#' With `fill = TRUE`, the balances must be compatible with a hierarchical SBP:
#' each pair must involve disjoint sets of parts, or one balance must lie wholly
#' within one group of the other. Completion preserves the supplied columns and
#' appends balances to form a full orthonormal basis. General, non-hierarchical
#' collections of balances remain supported with `fill = FALSE`.
#'
#' @param sbp A non-empty list of formulas or a numeric matrix with finite
#'   entries describing balances. At least two parts are required. A matrix
#'   with zero columns can be supplied to construct a full basis with
#'   `fill = TRUE`.
#' @param data Optional compositional data set with column names, or a character
#'   vector of part names, used when `sbp` is a list of formulas. If `data = NULL`,
#'   part names are inferred from the formulas themselves.
#' @param fill Logical; if `TRUE`, complete the supplied balances to obtain a
#'   full orthonormal basis. Only compatible hierarchical partitions can be
#'   completed.
#' @param silent Logical; if `FALSE`, report whether the resulting balances form
#'   a basis, and whether they are orthogonal or orthonormal.
#'
#' @return A numeric matrix whose columns are balances. Unnamed supplied balances
#'   receive names `b1`, `b2`, etc.; appended balances receive names `.b1`, `.b2`,
#'   etc. Names are made unique without changing supplied balance names.
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
#' # Formula parts can be inferred when data is not supplied
#' sbp_basis(list(
#'   b1 = a ~ b,
#'   b2 = b1 ~ c
#' ))
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
  if (!is.logical(fill) || length(fill) != 1L || is.na(fill) ||
      !is.logical(silent) || length(silent) != 1L || is.na(silent)) {
    stop("'fill' and 'silent' must each be TRUE or FALSE.")
  }
  if (is.matrix(sbp)) {
    if (!is.numeric(sbp) || any(!is.finite(sbp))) {
      stop("'sbp' must be a numeric matrix with finite entries.")
    }
    if (nrow(sbp) < 2L) stop("At least two parts are required.")
    signs <- sign(sbp)
    if (is.null(rownames(signs))) rownames(signs) <- paste0("c", seq_len(nrow(signs)))
    if (is.null(colnames(signs)) && ncol(signs) > 0L) {
      colnames(signs) <- paste0("b", seq_len(ncol(signs)))
    }
  } else {
    signs <- .sbp_formula_matrix(sbp, data)
  }
  if (anyNA(rownames(signs)) || any(!nzchar(rownames(signs))) ||
      anyDuplicated(rownames(signs))) {
    stop("Part names must be non-missing, non-empty and unique.")
  }
  if (anyNA(colnames(signs)) || any(!nzchar(colnames(signs))) ||
      anyDuplicated(colnames(signs))) {
    stop("Balance names must be non-missing, non-empty and unique.")
  }
  if (fill) {
    if (any(colSums(signs > 0) == 0 | colSums(signs < 0) == 0)) {
      stop("Each balance must contain at least one positive and one negative part.")
    }
    parts <- rownames(signs)
    balance_names <- colnames(signs)
    K <- ncol(signs)
    signs <- fill_sbp(signs)
    extra <- ncol(signs) - K
    dimnames(signs) <- list(parts, make.unique(c(
      balance_names, if (extra > 0L) paste0(".b", K + seq_len(extra))
    )))
  }

  # Normalize once for both input formats, after completing the sign matrix.
  RES <- matrix(0, nrow(signs), ncol(signs), dimnames = dimnames(signs))
  for (j in seq_len(ncol(signs))) {
    positive <- signs[, j] > 0
    negative <- signs[, j] < 0
    p <- sum(positive)
    n <- sum(negative)
    if (p == 0L || n == 0L) {
      stop("Each balance must contain at least one positive and one negative part.")
    }
    total <- as.double(p) + n
    RES[positive, j] <- sqrt(n / (p * total))
    RES[negative, j] <- -sqrt(p / (n * total))
  }

  # Completion constructs an orthonormal SBP. For other full bases, an
  # identity Gram matrix proves full rank and avoids an unnecessary QR.
  if (!silent && !fill) {
    if (ncol(RES) != nrow(RES) - 1L) {
      warning("Given partition is not a basis")
    } else {
      Z <- crossprod(RES)
      tol <- 1e-10
      off <- Z
      diag(off) <- 0
      orthogonal <- !any(abs(off) > tol)
      unit <- !any(abs(diag(Z) - 1) > tol)
      if (!orthogonal || !unit) {
        if (qr(RES)$rank != nrow(RES) - 1L) {
          warning("Given partition is not a basis")
        } else if (!orthogonal) {
          warning("Given basis is not orthogonal")
        } else {
          warning("Given basis is not orthonormal")
        }
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
