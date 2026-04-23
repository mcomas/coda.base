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
#' (`matrix` or `data.frame`, preserving data-frame subclasses when possible)
#' and preserving original names.
#' If `parameters = TRUE`, a list with the estimated clr mean, clr covariance,
#' and imputed clr coordinates.
#'
#' @examples
#' X <- matrix(c(
#'   0.00, 0.30, 0.70,
#'   0.20,   NA, 0.80,
#'   0.40, 0.60, 0.00,
#'   0.25, 0.25, 0.50,
#'   0.10, 0.30, 0.60
#' ), ncol = 3, byrow = TRUE)
#' colnames(X) <- c("sand", "silt", "clay")
#'
#' DL <- c(0.05, 0.05, 0.05)
#'
#' X_imp <- coda_replacement(X, DL = DL, maxit = 20)
#' X_imp
#'
#' set.seed(10)
#' X <- composition(matrix(rnorm(3*10), ncol = 3))
#' X[sample(c(TRUE, FALSE), 4*10, replace = TRUE, c(1, 3))] <- NA
#' params <- coda_replacement(X, parameters = TRUE, debug = TRUE)
#' names(params)
#' params$clr_mu
#' composition(params$clr_h)
#' @export
coda_replacement <- function(X, DL = NULL, dl_prop = 0.65,
                             eps = 1e-4, parameters = FALSE,
                             debug = FALSE, maxit = 500) {

  X_orig <- X

  if (is.matrix(X_orig)) {
    if (!is.numeric(X_orig)) {
      stop("Composition must be numeric.", call. = FALSE)
    }
    X_mat <- X_orig
  } else if (inherits(X_orig, "data.frame")) {
    if (!all(sapply(X_orig, is.numeric))) {
      stop("All parts must be numeric.", call. = FALSE)
    }
    class_type <- class(X_orig)
    X_mat <- as.matrix(X_orig)
  } else {
    stop("'X' must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (ncol(X_mat) < 2) {
    stop("'X' must contain at least two parts.", call. = FALSE)
  }
  if (any(!is.na(X_mat) & !is.finite(X_mat))) {
    stop("'X' must contain only finite values and missing values.", call. = FALSE)
  }
  if (any(X_mat < 0, na.rm = TRUE)) {
    stop("'X' must contain non-negative values.", call. = FALSE)
  }

  if (length(dl_prop) != 1 || !is.numeric(dl_prop) ||
      !is.finite(dl_prop) || dl_prop <= 0 || dl_prop >= 1) {
    stop("'dl_prop' must be a number between 0 and 1.", call. = FALSE)
  }
  if (length(eps) != 1 || !is.numeric(eps) || !is.finite(eps) || eps <= 0) {
    stop("'eps' must be a positive number.", call. = FALSE)
  }
  if (length(maxit) != 1 || !is.numeric(maxit) || !is.finite(maxit) || maxit < 1) {
    stop("'maxit' must be a positive integer.", call. = FALSE)
  }
  if (length(parameters) != 1 || !is.logical(parameters) || is.na(parameters)) {
    stop("'parameters' must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(debug) != 1 || !is.logical(debug) || is.na(debug)) {
    stop("'debug' must be TRUE or FALSE.", call. = FALSE)
  }

  maxit <- as.integer(maxit)
  tX <- t(unname(X_mat))

  if (is.null(DL)) {
    DL <- apply(tX, 1, function(x) {
      positive <- x[!is.na(x) & x > 0]
      if (length(positive) == 0L) NA_real_ else min(positive)
    })

    if (anyNA(DL) || any(!is.finite(DL))) {
      stop(
        "Could not estimate detection limits: each part must contain at least ",
        "one positive observed value, or 'DL' must be supplied.",
        call. = FALSE
      )
    }
  }

  if (is.vector(DL)) {
    if (length(DL) != nrow(tX)) {
      stop("'DL' must have length ncol(X) or the same dimensions as 'X'.", call. = FALSE)
    }
    tDL <- matrix(DL, nrow = nrow(tX), ncol = ncol(tX))
  } else if (is.data.frame(DL) || is.matrix(DL)) {
    tDL <- t(as.matrix(DL))
  } else {
    stop("Detection limit parameter (DL) must be a vector, a matrix or NULL.")
  }

  if (!is.numeric(tDL)) {
    stop("Detection limits must be numeric.", call. = FALSE)
  }
  if (!all(dim(tDL) == dim(tX))) {
    stop("'DL' must have length ncol(X) or the same dimensions as 'X'.", call. = FALSE)
  }
  if (any(!is.finite(tDL)) || any(tDL <= 0)) {
    stop("Detection limits must be strictly positive finite values.", call. = FALSE)
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

    if (inherits(X_orig, "data.frame")) {
      res <- as.data.frame(res_mat, stringsAsFactors = FALSE, check.names = FALSE)
      class(res) <- class_type
    } else {
      res <- res_mat
    }
  }

  res
}
