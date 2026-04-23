#' Tabu search for a partial principal balance on grouped parts
#'
#' Finds a single grouped balance by tabu search over a partition of selected
#' parts. The search is carried out on groups of parts defined by \code{lI},
#' using the variance criterion induced by the covariance matrix of the
#' log-transformed composition.
#'
#' The initial grouped split is obtained either from the constrained principal
#' balance of the grouped subcomposition, or from a user-supplied grouped
#' balance.
#'
#' @param X A numeric matrix with strictly positive finite entries. Rows are
#'   observations and columns are compositional parts.
#' @param lI A list defining a partition of a subset of the columns of
#'   \code{X}. Each element of \code{lI} is an integer vector giving the
#'   indices of the parts belonging to the same group. Indices must be valid
#'   column positions of \code{X}, and no column may appear in more than one
#'   group.
#' @param iter Integer. Maximum number of tabu search iterations.
#' @param tabu_size Integer. Maximum size of the tabu list.
#' @param ini Initial grouped split. If \code{NULL} (default), the initial
#'   solution is obtained from the constrained principal balance of the grouped
#'   subcomposition. Otherwise, \code{ini} must be an integer/numeric vector in
#'   \eqn{\{-1,0,1\}} of length \code{length(lI)} defining the initial grouped
#'   split. Negative entries indicate groups on the left side, positive entries
#'   indicate groups on the right side, and zeros indicate inactive groups.
#' @param debug Logical. If \code{TRUE}, progress information is printed during
#'   the search.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{iter_best}}{Iteration at which the best solution was found.}
#'   \item{\code{tabu_size}}{Effective tabu list size when the best solution
#'   was found.}
#'   \item{\code{steps}}{Objective values along the visited search path.}
#'   \item{\code{dim}}{Dimension of the grouped problem, equal to
#'   \code{length(lI) - 1}.}
#'   \item{\code{lI}}{The input grouping structure.}
#'   \item{\code{variance}}{Variance criterion of the best grouped balance
#'   found.}
#'   \item{\code{balance_raw}}{Integer vector in \eqn{\{-1,0,1\}} describing
#'   the best grouped split found. Negative entries indicate groups on the left
#'   side of the balance, positive entries indicate groups on the right side,
#'   and zeros indicate inactive groups.}
#'   \item{\code{balance}}{The corresponding one-column balance basis obtained
#'   from \code{balance_raw} using \code{\link[coda.base]{sbp_basis}}.}
#' }
#'
#' @details
#' The objective function is the variance of the balance associated with a split
#' between a left and a right set of active groups. At each tabu iteration,
#' candidate neighbors are obtained by:
#' \itemize{
#'   \item removing one active group from the current split,
#'   \item adding one inactive group to the left side,
#'   \item adding one inactive group to the right side.
#' }
#'
#' When \code{ini = NULL}, the grouped composition used for initialization is
#' obtained by replacing each group with the product of its parts, and the
#' constrained principal balance of this grouped composition is used as the
#' initial grouped split. When \code{ini} is a vector, it is used directly as
#' the initial grouped split.
#'
#' @seealso \code{\link{pb_tabu_search}},
#'   \code{\link[coda.base]{pb_basis}},
#'   \code{\link[coda.base]{sbp_basis}}
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rexp(200), ncol = 4)
#' lI <- list(1, 2, c(3, 4))
#'
#' res1 <- partial_pb_tabu_search(
#'   X = X,
#'   lI = lI,
#'   iter = 20,
#'   tabu_size = 3
#' )
#'
#' res2 <- partial_pb_tabu_search(
#'   X = X,
#'   lI = lI,
#'   iter = 20,
#'   tabu_size = 3,
#'   ini = c(-1, 1, 0)
#' )
#'
#' res1$variance
#' res1$balance_raw
#' res1$balance
#'
#' @export
partial_pb_tabu_search <- function(X, lI, iter, tabu_size,
                                   ini = NULL,
                                   debug = FALSE) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix.")
  }
  if (ncol(X) < 2) {
    stop("X must have at least two columns.")
  }
  if (any(!is.finite(X)) || any(X <= 0)) {
    stop("X must contain strictly positive finite values.")
  }

  if (!is.list(lI) || length(lI) < 2) {
    stop("lI must be a list with at least two groups.")
  }
  if (!all(vapply(lI, function(x) is.numeric(x) && length(x) >= 1, logical(1)))) {
    stop("Each element of lI must be a non-empty numeric vector of column indices.")
  }

  idx <- unlist(lI, use.names = FALSE)

  if (length(idx) < 2) {
    stop("lI must contain at least two column indices in total.")
  }
  if (any(!is.finite(idx)) || any(idx != as.integer(idx))) {
    stop("All indices in lI must be finite integers.")
  }
  idx <- as.integer(idx)

  if (any(idx < 1L | idx > ncol(X))) {
    stop("All indices in lI must be valid column indices of X.")
  }
  if (anyDuplicated(idx)) {
    stop("lI must define a partition of a subset of the columns of X, without duplicates.")
  }

  if (length(iter) != 1 || !is.numeric(iter) || !is.finite(iter) || iter < 1) {
    stop("iter must be a positive integer.")
  }
  if (length(tabu_size) != 1 || !is.numeric(tabu_size) ||
      !is.finite(tabu_size) || tabu_size < 1) {
    stop("tabu_size must be a positive integer.")
  }

  iter <- as.integer(iter)
  tabu_size <- as.integer(tabu_size)

  X_sub <- X[, idx, drop = FALSE]
  M <- cov(log(X_sub))
  G <- length(lI)

  lI_sub <- vector("list", G)
  offset <- 0L
  for (g in seq_len(G)) {
    ng <- length(lI[[g]])
    lI_sub[[g]] <- seq.int(offset + 1L, offset + ng)
    offset <- offset + ng
  }

  if (is.null(ini)) {
    Xb <- do.call(cbind, lapply(lI_sub, function(I) {
      Reduce(`*`, lapply(I, function(i) X_sub[, i]))
    }))
    B0 <- pb_basis(Xb, method = "constrained")
    BAL0 <- as.integer(sign(B0[, 1]))
  } else {
    if (!is.numeric(ini)) {
      stop("If 'ini' is not NULL, it must be a numeric vector in {-1,0,1}.")
    }
    if (length(ini) != G) {
      stop("User-supplied 'ini' must have length equal to length(lI).")
    }
    if (!all(is.finite(ini)) || !all(ini %in% c(-1, 0, 1))) {
      stop("User-supplied 'ini' must contain only values in {-1,0,1}.")
    }

    BAL0 <- as.integer(ini)

    if (!any(BAL0 < 0L) || !any(BAL0 > 0L)) {
      stop("User-supplied 'ini' must contain at least one -1 and one 1.")
    }
  }

  res <- partial_pb_tabu_search_cpp(
    M = M,
    lI = lI_sub,
    bal0 = BAL0,
    iter = iter,
    tabu_size = tabu_size,
    debug = debug
  )

  res$lI <- lI
  res$balance_raw <- as.integer(res$balance_raw)
  res$balance <- sbp_basis(cbind(res$balance_raw), silent = TRUE)

  res
}

#' Tabu search approximation to a sequential binary partition
#'
#' Builds a sequential binary partition (SBP) by repeatedly applying grouped
#' tabu search to select balances over the current sets of parts. At each step,
#' the best candidate split is retained and the remaining candidate subproblems
#' are explored until an SBP with \eqn{D - 1} balances is obtained.
#'
#' This function provides a heuristic approximation to a principal balance
#' basis. The first balance is searched on the full set of parts, and the
#' subsequent balances are obtained by recursively refining the best currently
#' available split.
#'
#' All partial searches are initialized with the constrained principal balance
#' of the corresponding grouped composition.
#'
#' @param X A numeric matrix with strictly positive finite entries. Rows are
#'   observations and columns are compositional parts.
#' @param iter Integer. Maximum number of tabu search iterations used in each
#'   partial search.
#' @param debug Logical. If \code{TRUE}, progress information is printed during
#'   the partial searches.
#'
#' @return An integer matrix representing a sequential binary partition. Rows
#' correspond to the original parts of \code{X} and columns correspond to
#' balances. Entries are in \eqn{\{-1,0,1\}}. The returned matrix has attribute
#' \code{"max_steps"}, giving the largest iteration index at which a best
#' partial solution was found among all partial searches performed.
#'
#' @details
#' The procedure starts from the trivial grouping where each part forms its own
#' singleton group. After each partial tabu search, up to three candidate
#' subproblems may be generated from the selected solution:
#' \itemize{
#'   \item the split between active and inactive groups,
#'   \item the left active branch,
#'   \item the right active branch.
#' }
#'
#' All generated candidates are stored, and at each stage the candidate with
#' the largest variance criterion is selected for inclusion in the SBP and for
#' further refinement.
#'
#' This is a heuristic search strategy and does not guarantee a globally
#' optimal SBP.
#'
#' @seealso \code{\link{partial_pb_tabu_search}},
#'   \code{\link[coda.base]{sbp_basis}}
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rexp(500), ncol = 5)
#'
#' SBP <- pb_tabu_search(X, iter = 30)
#' SBP
#' attr(SBP, "max_steps")
#'
#' @export
pb_tabu_search <- function(X, iter = 100, debug = FALSE) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix.")
  }
  if (ncol(X) < 2) {
    stop("X must have at least two columns.")
  }
  if (any(!is.finite(X)) || any(X <= 0)) {
    stop("X must contain strictly positive finite values.")
  }
  if (length(iter) != 1 || !is.numeric(iter) || !is.finite(iter) || iter < 1) {
    stop("iter must be a positive integer.")
  }

  iter <- as.integer(iter)
  max_steps <- 0L

  best <- partial_pb_tabu_search(
    X = X,
    lI = lapply(seq_len(ncol(X)), identity),
    iter = iter,
    tabu_size = ncol(X),
    debug = debug
  )

  max_steps <- max(max_steps, as.integer(best$iter_best))

  PB <- list(best)
  SOLS <- list()

  for (k in seq_len(max(0, ncol(X) - 2))) {

    if (any(best$balance_raw == 0L)) {
      lI_top <- c(
        list(unlist(best$lI[best$balance_raw != 0L])),
        lapply(which(best$balance_raw == 0L), function(i) best$lI[[i]])
      )

      top <- partial_pb_tabu_search(
        X = X,
        lI = lI_top,
        iter = iter,
        tabu_size = length(lI_top),
        debug = debug
      )

      max_steps <- max(max_steps, as.integer(top$iter_best))
      SOLS[[length(SOLS) + 1L]] <- top
    }

    if (sum(best$balance_raw < 0L) > 1L) {
      lI_left <- best$lI[best$balance_raw < 0L]

      left <- partial_pb_tabu_search(
        X = X,
        lI = lI_left,
        iter = iter,
        tabu_size = length(lI_left),
        debug = debug
      )

      max_steps <- max(max_steps, as.integer(left$iter_best))
      SOLS[[length(SOLS) + 1L]] <- left
    }

    if (sum(best$balance_raw > 0L) > 1L) {
      lI_right <- best$lI[best$balance_raw > 0L]

      right <- partial_pb_tabu_search(
        X = X,
        lI = lI_right,
        iter = iter,
        tabu_size = length(lI_right),
        debug = debug
      )

      max_steps <- max(max_steps, as.integer(right$iter_best))
      SOLS[[length(SOLS) + 1L]] <- right
    }

    if (length(SOLS) == 0L) {
      break
    }

    i.best <- which.max(vapply(SOLS, function(sol) sol$variance, numeric(1)))
    best <- SOLS[[i.best]]
    PB[[length(PB) + 1L]] <- best
    SOLS <- SOLS[-i.best]
  }

  SBP <- do.call(cbind, lapply(PB, function(pb) {
    b <- rep(0L, ncol(X))
    b[unlist(pb$lI[pb$balance_raw > 0L])] <- 1L
    b[unlist(pb$lI[pb$balance_raw < 0L])] <- -1L
    b
  }))

  attr(SBP, "max_steps") <- max_steps
  SBP
}

#' Generate a random composition with a prescribed first principal balance
#'
#' Simulates a random composition whose coordinate system is constructed from
#' a sequential binary partition induced by a given first balance. The supplied
#' balance is completed to a full orthonormal basis using
#' \code{\link[coda.base]{sbp_basis}} with \code{fill = TRUE}.
#'
#' @param principal_balance An integer or numeric vector in
#'   \eqn{\{-1,0,1\}} defining the first balance of the SBP.
#' @param n Integer. Number of observations to generate.
#' @param sd1 Numeric value used to scale the first latent coordinate before
#'   rotating the simulated coordinates.
#'
#' @return A composition matrix with \code{n} rows and
#'   \code{length(principal_balance)} columns.
#'
#' @details
#' Standard normal latent coordinates are first generated in dimension
#' \eqn{D - 1}, where \eqn{D} is the number of parts. Their sample covariance
#' matrix is then diagonalized, and the associated eigenvectors are used to
#' rotate the latent coordinates before mapping them back to the simplex using
#' the basis induced by \code{principal_balance}.
#'
#' This function is mainly intended for examples, simulation studies, and
#' experiments where a specific first balance structure is desired.
#'
#' @seealso \code{\link[coda.base]{sbp_basis}},
#'   \code{\link[coda.base]{composition}}
#' @export
random_composition_with_fixed_pb <- function(principal_balance, n = 100, sd1 = 5) {
  if (!is.numeric(principal_balance)) {
    stop("principal_balance must be a numeric vector.")
  }
  if (!all(principal_balance %in% c(-1, 0, 1))) {
    stop("principal_balance must contain only -1, 0, and 1.")
  }
  if (!any(principal_balance == -1) || !any(principal_balance == 1)) {
    stop("principal_balance must contain at least one -1 and one 1.")
  }
  if (length(n) != 1 || !is.numeric(n) || !is.finite(n) || n < 1) {
    stop("n must be a positive integer.")
  }

  n <- as.integer(n)
  D <- length(principal_balance)
  d <- D - 1L

  H <- matrix(rnorm(n * d), nrow = n)
  H[,1] = H[,1] * sd1
  S <- cov(H)
  EIG <- eigen(S)

  PB <- sbp_basis(cbind(as.integer(principal_balance)), fill = TRUE)
  X <- composition(H %*% EIG$vectors, basis = PB)

  X
}
