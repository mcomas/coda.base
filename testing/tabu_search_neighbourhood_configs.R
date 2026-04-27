## Compare neighbourhood configurations for partial principal-balance tabu search.
##
## Run from the package root:
##   source("testing/tabu_search_neighbourhood_configs.R")

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(coda.base)
}

make_configs <- function() {
  list(
    original = list(
      remove_active = TRUE,
      add_left = TRUE,
      add_right = TRUE,
      flip_side = FALSE,
      swap_zero = FALSE,
      swap_sides = FALSE
    ),
    flip = list(
      remove_active = TRUE,
      add_left = TRUE,
      add_right = TRUE,
      flip_side = TRUE,
      swap_zero = FALSE,
      swap_sides = FALSE
    ),
    swap_zero = list(
      remove_active = TRUE,
      add_left = TRUE,
      add_right = TRUE,
      flip_side = FALSE,
      swap_zero = TRUE,
      swap_sides = FALSE
    ),
    swap_sides = list(
      remove_active = TRUE,
      add_left = TRUE,
      add_right = TRUE,
      flip_side = FALSE,
      swap_zero = FALSE,
      swap_sides = TRUE
    ),
    swaps = list(
      remove_active = TRUE,
      add_left = TRUE,
      add_right = TRUE,
      flip_side = FALSE,
      swap_zero = TRUE,
      swap_sides = TRUE
    ),
    all = list(
      remove_active = TRUE,
      add_left = TRUE,
      add_right = TRUE,
      flip_side = TRUE,
      swap_zero = TRUE,
      swap_sides = TRUE
    ),
    swaps_only = list(
      remove_active = FALSE,
      add_left = FALSE,
      add_right = FALSE,
      flip_side = FALSE,
      swap_zero = TRUE,
      swap_sides = TRUE
    )
  )
}

balance_code <- function(b) {
  paste(ifelse(b < 0, "-", ifelse(b > 0, "+", ".")), collapse = "")
}

parse_balance_code <- function(code) {
  x <- strsplit(code, "", fixed = TRUE)[[1]]
  ifelse(x == "-", -1L, ifelse(x == "+", 1L, 0L))
}

same_balance <- function(a, b) {
  identical(as.integer(a), as.integer(b)) ||
    identical(as.integer(a), as.integer(-b))
}

make_problem <- function(seed, n = 100, D = 12, type = c("rlnorm", "fixed_pb")) {
  type <- match.arg(type)
  set.seed(seed)

  if (type == "rlnorm") {
    return(matrix(rlnorm(n * D), nrow = n, ncol = D))
  }

  pb1 <- sample(c(-1L, 0L, 1L), D, replace = TRUE)
  while (!any(pb1 < 0L) || !any(pb1 > 0L)) {
    pb1 <- sample(c(-1L, 0L, 1L), D, replace = TRUE)
  }

  random_composition_with_fixed_pb(pb1, n = n)
}

run_config <- function(X,
                       config_name,
                       config,
                       iter = 100,
                       tabu_size = ncol(X),
                       min_parts = 2,
                       max_parts = NULL) {
  lI <- lapply(seq_len(ncol(X)), identity)
  t0 <- proc.time()[["elapsed"]]

  res <- tryCatch(
    do.call(
      partial_pb_tabu_search,
      c(
        list(
          X = X,
          lI = lI,
          min_parts = min_parts,
          max_parts = max_parts,
          iter = iter,
          tabu_size = tabu_size
        ),
        config
      )
    ),
    error = function(e) e
  )

  elapsed <- proc.time()[["elapsed"]] - t0

  if (inherits(res, "error")) {
    return(data.frame(
      config = config_name,
      ok = FALSE,
      error = conditionMessage(res),
      variance = NA_real_,
      iter_best = NA_integer_,
      tabu_size_best = NA_integer_,
      elapsed = elapsed,
      n_active = NA_integer_,
      n_left = NA_integer_,
      n_right = NA_integer_,
      balance = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    config = config_name,
    ok = TRUE,
    error = NA_character_,
    variance = res$variance,
    iter_best = res$iter_best,
    tabu_size_best = res$tabu_size,
    elapsed = elapsed,
    n_active = sum(res$balance_raw != 0L),
    n_left = sum(res$balance_raw < 0L),
    n_right = sum(res$balance_raw > 0L),
    balance = balance_code(res$balance_raw),
    stringsAsFactors = FALSE
  )
}

run_one_problem <- function(seed,
                            n = 100,
                            D = 12,
                            iter = 100,
                            tabu_size = D,
                            min_parts = 2,
                            max_parts = NULL,
                            type = "rlnorm",
                            exact_max_D = 15,
                            configs = make_configs()) {
  X <- make_problem(seed = seed, n = n, D = D, type = type)

  constrained <- pb_basis(X, method = "constrained")
  constrained_raw <- as.integer(sign(constrained[, 1]))
  var_constrained <- as.numeric(var(coordinates(X, constrained[, 1, drop = FALSE])))

  exact_raw <- rep(NA_integer_, D)
  var_exact <- NA_real_
  if (D <= exact_max_D) {
    exact <- pb_basis(X, method = "exact")
    exact_raw <- as.integer(sign(exact[, 1]))
    var_exact <- as.numeric(var(coordinates(X, exact[, 1, drop = FALSE])))
  }

  out <- do.call(
    rbind,
    Map(
      function(config_name, config) {
        run_config(
          X = X,
          config_name = config_name,
          config = config,
          iter = iter,
          tabu_size = tabu_size,
          min_parts = min_parts,
          max_parts = max_parts
        )
      },
      names(configs),
      configs
    )
  )

  n_out <- nrow(out)
  out$seed <- rep(seed, n_out)
  out$n <- rep(n, n_out)
  out$D <- rep(D, n_out)
  out$type <- rep(type, n_out)
  out$min_parts <- rep(min_parts, n_out)
  out$max_parts <- rep(if (is.null(max_parts)) D else max_parts, n_out)
  out$var_constrained <- rep(var_constrained, n_out)
  out$var_exact <- rep(var_exact, n_out)
  out$gain_vs_constrained <- out$variance - var_constrained
  out$gap_vs_exact <- if (is.na(var_exact)) {
    rep(NA_real_, n_out)
  } else {
    var_exact - out$variance
  }
  out$same_as_constrained <- vapply(
    out$balance,
    function(code) {
      if (is.na(code)) {
        return(NA)
      }
      same_balance(parse_balance_code(code), constrained_raw)
    },
    logical(1)
  )
  out$same_as_exact <- if (all(is.na(exact_raw))) {
    rep(NA, n_out)
  } else {
    vapply(
      out$balance,
      function(code) {
        if (is.na(code)) {
          return(NA)
        }
        same_balance(parse_balance_code(code), exact_raw)
      },
      logical(1)
    )
  }

  out[order(!out$ok, -out$variance, out$elapsed), ]
}

run_neighbourhood_benchmark <- function(seeds = 1:30,
                                        n = 100,
                                        D = 12,
                                        iter = 100,
                                        tabu_size = D,
                                        min_parts = 2,
                                        max_parts = NULL,
                                        type = "rlnorm",
                                        exact_max_D = 15,
                                        configs = make_configs()) {
  res <- do.call(
    rbind,
    lapply(seeds, function(seed) {
      cat(sprintf("seed=%d, D=%d, iter=%d\n", seed, D, iter))
      run_one_problem(
        seed = seed,
        n = n,
        D = D,
        iter = iter,
        tabu_size = tabu_size,
        min_parts = min_parts,
        max_parts = max_parts,
        type = type,
        exact_max_D = exact_max_D,
        configs = configs
      )
    })
  )

  rownames(res) <- NULL
  res
}

summarise_neighbourhood_benchmark <- function(res) {
  configs <- unique(res$config)
  best_by_seed <- tapply(res$variance, res$seed, function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    max(x, na.rm = TRUE)
  })

  summary <- do.call(
    rbind,
    lapply(configs, function(config) {
      x <- res[res$config == config, , drop = FALSE]
      x_ok <- x[x$ok, , drop = FALSE]
      if (nrow(x_ok) == 0L) {
        return(data.frame(
          config = config,
          ok_rate = 0,
          first_error = x$error[which(!is.na(x$error))[1]],
          mean_variance = NA_real_,
          mean_gain_vs_constrained = NA_real_,
          mean_gap_vs_exact = NA_real_,
          wins = 0L,
          mean_iter_best = NA_real_,
          mean_elapsed = mean(x$elapsed),
          mean_active = NA_real_,
          same_as_constrained = NA_real_,
          same_as_exact = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      data.frame(
        config = config,
        ok_rate = mean(x$ok),
        first_error = NA_character_,
        mean_variance = mean(x_ok$variance),
        mean_gain_vs_constrained = mean(x_ok$gain_vs_constrained),
        mean_gap_vs_exact = mean(x_ok$gap_vs_exact, na.rm = TRUE),
        wins = sum(x_ok$variance == best_by_seed[as.character(x_ok$seed)]),
        mean_iter_best = mean(x_ok$iter_best),
        mean_elapsed = mean(x$elapsed),
        mean_active = mean(x_ok$n_active),
        same_as_constrained = mean(x_ok$same_as_constrained, na.rm = TRUE),
        same_as_exact = mean(x_ok$same_as_exact, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )

  summary[order(-summary$mean_variance, summary$mean_elapsed), ]
}

## Example run. Adjust seeds, D, iter and type for heavier experiments.
RES <- run_neighbourhood_benchmark(
  seeds = 1:20,
  n = 100,
  D = 12,
  iter = 100,
  tabu_size = 12,
  max_parts = NULL,
  type = "rlnorm"
)

SUMMARY <- summarise_neighbourhood_benchmark(RES)
print(SUMMARY)

## Quick visual checks when running interactively.
if (interactive()) {
  op <- par(mfrow = c(1, 2))
  boxplot(variance ~ config, data = RES, las = 2, main = "Variance")
  boxplot(elapsed ~ config, data = RES, las = 2, main = "Elapsed time")
  par(op)
}
