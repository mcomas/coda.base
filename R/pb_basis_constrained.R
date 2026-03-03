#' Recursive constrained principal balances on subcompositions
#'
#' Recursively construct balances on selected subcompositions, optionally
#' enforcing groups of variables to remain together through constraints.
#'
#' @param X Compositional data set.
#' @param variables Indices of the variables currently considered.
#' @param constraints Optional list of groups of variables to be constrained
#'   together during the recursive search.
#' @param angle Logical; if `TRUE`, use the angle criterion instead of the
#'   variance criterion when computing constrained balances.
#'
#' @return A list of balance vectors.
#'
#' @keywords internal
pb_subcomposition <- function(X,
                              variables = seq_len(ncol(X)),
                              constraints = NULL,
                              angle = FALSE) {
  Xc <- X

  if (!is.null(constraints)) {
    for (constraint in constraints) {
      Xc[, constraint] <- exp(rowMeans(log(Xc[, constraint])))
    }
  }

  pb <- matrix(0, nrow = ncol(X), ncol = 1)
  pb[variables, 1] <- get_balance_using_pc(Xc[, variables], angle)[, 1]

  lpb <- list(pb)

  if (sum(pb < 0) > 1) {
    if (max(apply(Xc[, pb < 0, drop = FALSE], 1, var)) > 0) {
      lpb <- c(
        lpb,
        Recall(X, variables = which(pb < 0), constraints = constraints, angle = angle)
      )
    }
  }

  if (sum(pb > 0) > 1) {
    if (max(apply(Xc[, pb > 0, drop = FALSE], 1, var)) > 0) {
      lpb <- c(
        lpb,
        Recall(X, variables = which(pb > 0), constraints = constraints, angle = angle)
      )
    }
  }

  sel <- rep(FALSE, ncol(X))
  sel[variables] <- TRUE

  if (sum((pb == 0) & sel) > 0) {
    if (is.null(constraints)) {
      constraints <- list()
    }

    constraints[[length(constraints) + 1]] <- (pb != 0) & sel

    lpb <- c(
      lpb,
      Recall(X, variables = variables, constraints = constraints, angle = angle)
    )
  }

  lpb
}

#' Constrained principal balance basis
#'
#' Compute a basis of constrained principal balances recursively.
#'
#' @param X Compositional data set.
#' @param angle Logical; if `TRUE`, use the angle criterion instead of the
#'   variance criterion.
#'
#' @return A matrix whose columns are constrained principal balances.
#'
#' @keywords internal
constrained_pb <- function(X, angle = FALSE) {
  l_pb <- pb_subcomposition(X, angle = angle)
  sapply(l_pb, identity)
}
