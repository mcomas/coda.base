balance_make_d_moves <- function(b0, d, ensure_balance = TRUE, max_tries = 1000) {
  if (!is.numeric(b0) || !all(b0 %in% c(-1, 0, 1))) {
    stop("b0 must be a numeric vector with entries in {-1, 0, 1}.")
  }

  n <- length(b0)

  if (length(d) != 1 || !is.numeric(d) || !is.finite(d) || d < 0) {
    stop("d must be a non-negative integer.")
  }

  d <- as.integer(d)
  vals <- c(-1L, 0L, 1L)

  for (tries in seq_len(max_tries)) {
    b <- as.integer(b0)

    if (d > 0L) {
      for (k in seq_len(d)) {
        j <- sample.int(n, 1)
        b[j] <- sample(vals[vals != b[j]], 1)
      }
    }

    if (!ensure_balance || (any(b == -1L) && any(b == 1L))) {
      return(b)
    }
  }

  stop("Could not generate a feasible balance after max_tries attempts.")
}
