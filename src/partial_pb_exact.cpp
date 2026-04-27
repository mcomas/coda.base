#include <RcppArmadillo.h>
#include <vector>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

std::vector<arma::uvec> list_to_groups_partial_exact_cpp(const Rcpp::List& lI) {
  std::vector<arma::uvec> groups;
  groups.reserve(lI.size());

  for (int i = 0; i < lI.size(); ++i) {
    Rcpp::IntegerVector idx = lI[i];
    arma::uvec g(idx.size());

    for (int j = 0; j < idx.size(); ++j) {
      if (idx[j] < 1) {
        Rcpp::stop("All indices in lI must be >= 1.");
      }
      g[j] = static_cast<arma::uword>(idx[j] - 1);
    }

    groups.push_back(g);
  }

  return groups;
}

arma::mat grouped_covariance_partial_exact_cpp(
    const arma::mat& M,
    const std::vector<arma::uvec>& groups) {
  const arma::uword G = groups.size();
  arma::mat Mg(G, G, arma::fill::zeros);

  for (arma::uword i = 0; i < G; ++i) {
    for (arma::uword j = i; j < G; ++j) {
      double val = arma::accu(M.submat(groups[i], groups[j]));
      Mg(i, j) = val;
      Mg(j, i) = val;
    }
  }

  return Mg;
}

arma::vec group_sizes_partial_exact_cpp(const std::vector<arma::uvec>& groups) {
  arma::vec out(groups.size());
  for (arma::uword i = 0; i < groups.size(); ++i) {
    out[i] = static_cast<double>(groups[i].n_elem);
  }
  return out;
}

double evaluate_group_balance_partial_exact_cpp(
    const arma::ivec& bal,
    const arma::mat& Mg,
    const arma::vec& group_sizes) {
  double nL = 0.0;
  double nR = 0.0;

  for (arma::uword i = 0; i < bal.n_elem; ++i) {
    if (bal[i] < 0) {
      nL += group_sizes[i];
    } else if (bal[i] > 0) {
      nR += group_sizes[i];
    }
  }

  if (nL <= 0.0 || nR <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  }

  double sumLL = 0.0;
  double sumRR = 0.0;
  double sumLR = 0.0;

  for (arma::uword i = 0; i < bal.n_elem; ++i) {
    if (bal[i] == 0) continue;

    for (arma::uword j = 0; j < bal.n_elem; ++j) {
      if (bal[j] == 0) continue;

      if (bal[i] < 0 && bal[j] < 0) {
        sumLL += Mg(i, j);
      } else if (bal[i] > 0 && bal[j] > 0) {
        sumRR += Mg(i, j);
      } else if (bal[i] > 0 && bal[j] < 0) {
        sumLR += Mg(i, j);
      }
    }
  }

  return ((nR / nL) * sumLL + (nL / nR) * sumRR - 2.0 * sumLR) / (nL + nR);
}

void update_best_partial_exact_cpp(
    const arma::ivec& bal,
    const arma::mat& Mg,
    const arma::vec& group_sizes,
    double& best_variance,
    arma::ivec& best_balance) {
  double variance = evaluate_group_balance_partial_exact_cpp(bal, Mg, group_sizes);

  if (variance > best_variance) {
    best_variance = variance;
    best_balance = bal;
  }
}

void enumerate_signs_gray_partial_exact_cpp(
    const std::vector<int>& support,
    const arma::mat& Mg,
    const arma::vec& group_sizes,
    double& best_variance,
    arma::ivec& best_balance) {
  const int k = static_cast<int>(support.size());
  const int free_bits = k - 1;
  const unsigned long long n_patterns = 1ULL << free_bits;
  arma::ivec bal(Mg.n_rows, arma::fill::zeros);

  bal[support[0]] = -1;

  for (unsigned long long p = 1ULL; p < n_patterns; ++p) {
    unsigned long long gray = p ^ (p >> 1);

    for (int j = 0; j < free_bits; ++j) {
      const bool bit = ((gray >> j) & 1ULL) != 0ULL;
      bal[support[j + 1]] = bit ? 1 : -1;
    }

    update_best_partial_exact_cpp(
      bal,
      Mg,
      group_sizes,
      best_variance,
      best_balance
    );
  }
}

void enumerate_supports_partial_exact_cpp(
    const int G,
    const int k,
    const int start,
    std::vector<int>& support,
    const arma::mat& Mg,
    const arma::vec& group_sizes,
    double& best_variance,
    arma::ivec& best_balance) {
  if (static_cast<int>(support.size()) == k) {
    enumerate_signs_gray_partial_exact_cpp(
      support,
      Mg,
      group_sizes,
      best_variance,
      best_balance
    );
    return;
  }

  const int need = k - static_cast<int>(support.size());
  for (int i = start; i <= G - need; ++i) {
    support.push_back(i);
    enumerate_supports_partial_exact_cpp(
      G,
      k,
      i + 1,
      support,
      Mg,
      group_sizes,
      best_variance,
      best_balance
    );
    support.pop_back();
  }
}

// [[Rcpp::export]]
Rcpp::List partial_pb_exact_cpp(const arma::mat& M,
                                const Rcpp::List& lI,
                                const int min_parts,
                                const int max_parts) {
  if (M.n_rows != M.n_cols) {
    Rcpp::stop("M must be a square matrix.");
  }
  if (min_parts < 2) {
    Rcpp::stop("min_parts must be >= 2.");
  }
  if (max_parts < min_parts) {
    Rcpp::stop("max_parts must be >= min_parts.");
  }

  std::vector<arma::uvec> groups = list_to_groups_partial_exact_cpp(lI);
  const int G = static_cast<int>(groups.size());

  if (G < 2) {
    Rcpp::stop("lI must contain at least two groups.");
  }
  if (max_parts > G) {
    Rcpp::stop("max_parts cannot exceed length(lI).");
  }

  for (std::size_t g = 0; g < groups.size(); ++g) {
    if (groups[g].n_elem == 0) {
      Rcpp::stop("Groups in lI must be non-empty.");
    }
    if (arma::max(groups[g]) >= M.n_cols) {
      Rcpp::stop("An index in lI exceeds the dimension of M.");
    }
  }

  arma::mat Mg = grouped_covariance_partial_exact_cpp(M, groups);
  arma::vec sizes = group_sizes_partial_exact_cpp(groups);

  double best_variance = -std::numeric_limits<double>::infinity();
  arma::ivec best_balance(G, arma::fill::zeros);
  std::vector<int> support;

  for (int k = min_parts; k <= max_parts; ++k) {
    support.clear();
    enumerate_supports_partial_exact_cpp(
      G,
      k,
      0,
      support,
      Mg,
      sizes,
      best_variance,
      best_balance
    );
    Rcpp::checkUserInterrupt();
  }

  if (!std::isfinite(best_variance)) {
    Rcpp::stop("No feasible balance was found.");
  }

  return Rcpp::List::create(
    Rcpp::Named("dim") = G - 1,
    Rcpp::Named("variance") = best_variance,
    Rcpp::Named("balance_raw") = best_balance
  );
}
