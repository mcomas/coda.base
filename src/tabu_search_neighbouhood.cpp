#include <RcppArmadillo.h>
#include <deque>
#include <unordered_set>
#include <string>
#include <vector>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

std::string balance_key_neighbourhood_cpp(const arma::ivec& bal) {
  std::string key;
  key.reserve(bal.n_elem * 2);

  for (arma::uword i = 0; i < bal.n_elem; ++i) {
    if (bal[i] == -1) key += "-1";
    else if (bal[i] == 0) key += "0";
    else key += "1";
  }

  return key;
}

std::string format_balance_neighbourhood_cpp(const arma::ivec& bal) {
  std::string out;
  out.reserve(bal.n_elem * 2);

  for (arma::uword i = 0; i < bal.n_elem; ++i) {
    if (i > 0) out += "|";

    if (bal[i] == -1) out += "-";
    else if (bal[i] == 0) out += ".";
    else out += "+";
  }

  return out;
}

std::vector<arma::uvec> list_to_groups_neighbourhood_cpp(const Rcpp::List& lI) {
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

void expand_balance_indices_neighbourhood_cpp(
    const arma::ivec& bal,
    const std::vector<arma::uvec>& groups,
    arma::uvec& iL,
    arma::uvec& iR) {
  std::size_t nL = 0;
  std::size_t nR = 0;

  for (arma::uword j = 0; j < bal.n_elem; ++j) {
    if (bal[j] == -1) nL += groups[j].n_elem;
    else if (bal[j] == 1) nR += groups[j].n_elem;
  }

  iL.set_size(nL);
  iR.set_size(nR);

  std::size_t posL = 0;
  std::size_t posR = 0;

  for (arma::uword j = 0; j < bal.n_elem; ++j) {
    if (bal[j] == -1) {
      const arma::uvec& g = groups[j];
      iL.subvec(posL, posL + g.n_elem - 1) = g;
      posL += g.n_elem;
    } else if (bal[j] == 1) {
      const arma::uvec& g = groups[j];
      iR.subvec(posR, posR + g.n_elem - 1) = g;
      posR += g.n_elem;
    }
  }
}

double evaluate_balance_grouped_neighbourhood_cpp(
    const arma::ivec& bal,
    const arma::mat& M,
    const std::vector<arma::uvec>& groups) {
  arma::uvec iL, iR;
  expand_balance_indices_neighbourhood_cpp(bal, groups, iL, iR);

  const double nL = static_cast<double>(iL.n_elem);
  const double nR = static_cast<double>(iR.n_elem);

  if (nL == 0.0 || nR == 0.0) {
    return -std::numeric_limits<double>::infinity();
  }

  const double sumLL = arma::accu(M.submat(iL, iL));
  const double sumRR = arma::accu(M.submat(iR, iR));
  const double sumRL = arma::accu(M.submat(iR, iL));

  return ((nR / nL) * sumLL + (nL / nR) * sumRR - 2.0 * sumRL) / (nL + nR);
}

int active_groups_neighbourhood_cpp(const arma::ivec& bal) {
  int out = 0;

  for (arma::uword i = 0; i < bal.n_elem; ++i) {
    if (bal[i] != 0) {
      ++out;
    }
  }

  return out;
}

bool is_valid_balance_neighbourhood_cpp(const arma::ivec& bal,
                                        const std::vector<arma::uvec>& groups,
                                        const int min_parts,
                                        const int max_parts) {
  (void)groups;
  bool has_left = false;
  bool has_right = false;
  const int active = active_groups_neighbourhood_cpp(bal);

  for (arma::uword i = 0; i < bal.n_elem; ++i) {
    if (bal[i] == -1) has_left = true;
    else if (bal[i] == 1) has_right = true;
  }

  if (!has_left || !has_right) {
    return false;
  }

  return active >= min_parts && active <= max_parts;
}

void add_candidate_neighbourhood_cpp(std::vector<arma::ivec>& out,
                                     const arma::ivec& cand,
                                     const std::vector<arma::uvec>& groups,
                                     const int min_parts,
                                     const int max_parts) {
  if (is_valid_balance_neighbourhood_cpp(cand, groups, min_parts, max_parts)) {
    out.push_back(cand);
  }
}

std::vector<arma::ivec> neighbours_neighbourhood_cpp(
    const arma::ivec& bal,
    const bool remove_active,
    const bool add_left,
    const bool add_right,
    const bool flip_side,
    const bool swap_zero,
    const bool swap_sides,
    const std::vector<arma::uvec>& groups,
    const int min_parts,
    const int max_parts) {
  std::vector<arma::ivec> out;
  const arma::uword D = bal.n_elem;

  for (arma::uword i = 0; i < D; ++i) {
    if (remove_active && bal[i] != 0) {
      arma::ivec cand = bal;
      cand[i] = 0;
      add_candidate_neighbourhood_cpp(out, cand, groups, min_parts, max_parts);
    }

    if (add_left && bal[i] == 0) {
      arma::ivec cand = bal;
      cand[i] = -1;
      add_candidate_neighbourhood_cpp(out, cand, groups, min_parts, max_parts);
    }

    if (add_right && bal[i] == 0) {
      arma::ivec cand = bal;
      cand[i] = 1;
      add_candidate_neighbourhood_cpp(out, cand, groups, min_parts, max_parts);
    }

    if (flip_side && bal[i] != 0) {
      arma::ivec cand = bal;
      cand[i] = -cand[i];
      add_candidate_neighbourhood_cpp(out, cand, groups, min_parts, max_parts);
    }
  }

  if (swap_zero) {
    for (arma::uword i = 0; i < D; ++i) {
      if (bal[i] == 0) continue;

      for (arma::uword j = 0; j < D; ++j) {
        if (bal[j] != 0) continue;

        arma::ivec cand = bal;
        cand[j] = bal[i];
        cand[i] = 0;
        add_candidate_neighbourhood_cpp(out, cand, groups, min_parts, max_parts);
      }
    }
  }

  if (swap_sides) {
    for (arma::uword i = 0; i < D; ++i) {
      if (bal[i] != -1) continue;

      for (arma::uword j = 0; j < D; ++j) {
        if (bal[j] != 1) continue;

        arma::ivec cand = bal;
        cand[i] = 1;
        cand[j] = -1;
        add_candidate_neighbourhood_cpp(out, cand, groups, min_parts, max_parts);
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::List partial_pb_tabu_search_cpp(
    const arma::mat& M,
    const Rcpp::List& lI,
    const arma::ivec& bal0,
    const int iter,
    const int tabu_size,
    const bool remove_active = true,
    const bool add_left = true,
    const bool add_right = true,
    const bool flip_side = false,
    const bool swap_zero = false,
    const bool swap_sides = false,
    const int min_parts = 2,
    const int max_parts = -1,
    const bool debug = false) {
  if (M.n_rows != M.n_cols) {
    Rcpp::stop("M must be a square matrix.");
  }
  if (iter < 1) {
    Rcpp::stop("iter must be >= 1.");
  }
  if (tabu_size < 0) {
    Rcpp::stop("tabu_size must be >= 0.");
  }
  if (min_parts < 2) {
    Rcpp::stop("min_parts must be >= 2.");
  }
  if (max_parts != -1 && max_parts < 2) {
    Rcpp::stop("max_parts must be >= 2.");
  }
  if (bal0.n_elem != static_cast<arma::uword>(lI.size())) {
    Rcpp::stop("Length of bal0 must match length(lI).");
  }
  if (!remove_active && !add_left && !add_right && !flip_side &&
      !swap_zero && !swap_sides) {
    Rcpp::stop("At least one neighbourhood type must be active.");
  }

  std::vector<arma::uvec> groups = list_to_groups_neighbourhood_cpp(lI);

  for (std::size_t g = 0; g < groups.size(); ++g) {
    if (groups[g].n_elem == 0) {
      Rcpp::stop("Groups in lI must be non-empty.");
    }
    if (arma::max(groups[g]) >= M.n_cols) {
      Rcpp::stop("An index in lI exceeds the dimension of M.");
    }
  }

  int effective_max_parts = max_parts;
  if (effective_max_parts == -1) {
    effective_max_parts = static_cast<int>(groups.size());
  }
  if (effective_max_parts < min_parts) {
    Rcpp::stop("max_parts must be >= min_parts.");
  }
  if (effective_max_parts > static_cast<int>(groups.size())) {
    Rcpp::stop("max_parts cannot exceed the number of groups in lI.");
  }

  arma::ivec BAL = bal0;
  arma::vec steps(iter, arma::fill::zeros);

  std::deque<std::string> tabu_queue;
  std::unordered_set<std::string> tabu_set;

  arma::ivec BEST = BAL;
  int BEST_TABU_SIZE = 0;
  double BEST_EV = evaluate_balance_grouped_neighbourhood_cpp(BAL, M, groups);
  int BEST_ITER = 0;

  if (!std::isfinite(BEST_EV)) {
    Rcpp::stop("Initial balance is invalid: both sides must be non-empty.");
  }
  if (active_groups_neighbourhood_cpp(BAL) > effective_max_parts) {
    Rcpp::stop("Initial balance has more active groups than max_parts.");
  }
  if (active_groups_neighbourhood_cpp(BAL) < min_parts) {
    Rcpp::stop("Initial balance has fewer active groups than min_parts.");
  }

  if (debug) {
    Rcpp::Rcout << "Current variance: " << BEST_EV << "\n";
    Rcpp::Rcout << "Balance: " << format_balance_neighbourhood_cpp(BEST) << "\n\n";
  }

  if (tabu_size == 0) {
    std::vector<double> local_steps;
    int local_iter = 0;

    while (true) {
      double current_ev = evaluate_balance_grouped_neighbourhood_cpp(BAL, M, groups);
      local_steps.push_back(current_ev);

      std::vector<arma::ivec> BAL_N = neighbours_neighbourhood_cpp(
        BAL,
        remove_active,
        add_left,
        add_right,
        flip_side,
        swap_zero,
        swap_sides,
        groups,
        min_parts,
        effective_max_parts
      );

      if (BAL_N.empty()) {
        break;
      }

      double best_neigh_ev = -std::numeric_limits<double>::infinity();
      int best_neigh_idx = -1;

      for (std::size_t k = 0; k < BAL_N.size(); ++k) {
        double ev_k = evaluate_balance_grouped_neighbourhood_cpp(BAL_N[k], M, groups);

        if (ev_k > best_neigh_ev) {
          best_neigh_ev = ev_k;
          best_neigh_idx = static_cast<int>(k);
        }
      }

      if (best_neigh_idx < 0 || best_neigh_ev <= current_ev) {
        break;
      }

      BAL = BAL_N[best_neigh_idx];
      ++local_iter;

      if (debug) {
        Rcpp::Rcout << "Iteration " << local_iter
                    << " - current variance: " << best_neigh_ev << "\n";
        Rcpp::Rcout << "Balance: " << format_balance_neighbourhood_cpp(BAL) << "\n\n";
      }

      BEST_ITER = local_iter;
      BEST = BAL;
      BEST_EV = best_neigh_ev;
    }

    return Rcpp::List::create(
      Rcpp::Named("iter_best") = BEST_ITER,
      Rcpp::Named("tabu_size") = 0,
      Rcpp::Named("steps") = local_steps,
      Rcpp::Named("dim") = static_cast<int>(lI.size()) - 1,
      Rcpp::Named("variance") = BEST_EV,
      Rcpp::Named("balance_raw") = BEST
    );
  }

  for (int i = 0; i < iter; ++i) {
    steps[i] = evaluate_balance_grouped_neighbourhood_cpp(BAL, M, groups);

    std::string current_key = balance_key_neighbourhood_cpp(BAL);

    if (static_cast<int>(tabu_queue.size()) == tabu_size) {
      std::string old_key = tabu_queue.front();
      tabu_queue.pop_front();
      tabu_set.erase(old_key);
    }

    tabu_queue.push_back(current_key);
    tabu_set.insert(current_key);

    std::vector<arma::ivec> BAL_N = neighbours_neighbourhood_cpp(
      BAL,
      remove_active,
      add_left,
      add_right,
      flip_side,
      swap_zero,
      swap_sides,
      groups,
      min_parts,
      effective_max_parts
    );

    if (BAL_N.empty()) {
      break;
    }

    double best_neigh_ev = -std::numeric_limits<double>::infinity();
    int best_neigh_idx = -1;

    for (std::size_t k = 0; k < BAL_N.size(); ++k) {
      std::string key_k = balance_key_neighbourhood_cpp(BAL_N[k]);

      if (tabu_set.find(key_k) != tabu_set.end()) {
        continue;
      }

      double ev_k = evaluate_balance_grouped_neighbourhood_cpp(BAL_N[k], M, groups);

      if (ev_k > best_neigh_ev) {
        best_neigh_ev = ev_k;
        best_neigh_idx = static_cast<int>(k);
      }
    }

    if (best_neigh_idx < 0) {
      break;
    }

    BAL = BAL_N[best_neigh_idx];
    double BAL_EV = best_neigh_ev;

    if (debug) {
      Rcpp::Rcout << "Iteration " << (i + 1)
                  << " - current variance: " << BAL_EV << "\n";
      Rcpp::Rcout << "Balance: " << format_balance_neighbourhood_cpp(BAL) << "\n\n";
    }

    if (BAL_EV > BEST_EV) {
      BEST_ITER = i + 1;
      BEST = BAL;
      BEST_TABU_SIZE = static_cast<int>(tabu_queue.size());
      BEST_EV = BAL_EV;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("iter_best") = BEST_ITER,
    Rcpp::Named("tabu_size") = BEST_TABU_SIZE,
    Rcpp::Named("steps") = steps,
    Rcpp::Named("dim") = static_cast<int>(lI.size()) - 1,
    Rcpp::Named("variance") = BEST_EV,
    Rcpp::Named("balance_raw") = BEST
  );
}
