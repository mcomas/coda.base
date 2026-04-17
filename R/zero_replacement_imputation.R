#' Replacement of missing values and below-detection zeros in compositional data
#'
#' @description
#' Performs imputation of missing values and/or values below the detection limit
#' in compositional data using an EM algorithm assuming normality on the simplex.
#'
#' @param X A compositional data set: numeric matrix or data frame where rows
#'   represent observations and columns represent parts.
#' @param DL An optional matrix or vector of detection limits. If `NULL`, the
#'   minimum non-zero value in each column of `X` is used.
#' @param dl_prop A numeric value between 0 and 1 used for initialization in
#'   the EM algorithm.
#' @param eps Convergence tolerance.
#' @param parameters Logical; if `TRUE`, return additional estimated parameters.
#' @param debug Logical; if `TRUE`, print the log-likelihood at each iteration.
#' @param maxit Maximum number of iterations
#'
#' @return
#' If `parameters = FALSE`, the imputed object with the same format as `X`
#' (`matrix`, `data.frame` or `tibble`) and preserving original names.
#' If `parameters = TRUE`, a list with the estimated clr mean, clr covariance,
#' and imputed clr coordinates.
#'
#' @export
coda_replacement <- function(X, DL = NULL, dl_prop = 0.65,
                             eps = 1e-4, parameters = FALSE,
                             debug = FALSE, maxit = 500) {

  X_orig <- X
  x_is_tibble <- inherits(X_orig, "tbl_df")
  x_is_df <- is.data.frame(X_orig)
  x_is_matrix <- is.matrix(X_orig)

  X_mat <- as.matrix(X_orig)
  tX <- t(unname(X_mat))

  if (is.null(DL)) {
    DL <- apply(tX, 1, function(x) min(x[!is.na(x) & x > 0]))
  }

  if (is.vector(DL)) {
    tDL <- matrix(DL, nrow = nrow(tX), ncol = ncol(tX))
  } else if (is.data.frame(DL) || is.matrix(DL)) {
    tDL <- t(as.matrix(DL))
  } else {
    stop("Detection limit parameter (DL) must be a vector, a matrix or NULL.")
  }

  res <- c_coda_replacement(
    tX = tX,
    tDL = tDL,
    dl_prop = dl_prop,
    eps = eps,
    parameters = parameters,
    debug = debug,
    maxit = maxit
  )

  if (!parameters) {
    res_mat <- as.matrix(res)

    if (!all(dim(res_mat) == dim(X_mat))) {
      stop("The imputed object returned by `c_coda_replacement()` has incompatible dimensions.")
    }

    dimnames(res_mat) <- dimnames(X_mat)

    if (x_is_tibble) {
      res <- tibble::as_tibble(res_mat, .name_repair = "minimal")
    } else if (x_is_df) {
      res <- as.data.frame(res_mat, stringsAsFactors = FALSE, check.names = FALSE)
      rownames(res) <- rownames(X_orig)
    } else if (x_is_matrix) {
      res <- res_mat
    } else {
      res <- res_mat
    }
  }

  res
}
