#include <RcppArmadillo.h>
#include <deque>
#include <unordered_set>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------
// Build a unique key for a balance vector so it can be stored
// in the tabu list.
// Example: c(-1, 0, 1, 1) -> "-1011"
// ---------------------------------------------------------------
std::string balance_key_cpp(const arma::ivec& bal) {
  std::string key;
  key.reserve(bal.n_elem * 2);

  for (arma::uword i = 0; i < bal.n_elem; ++i) {
    if (bal[i] == -1) key += "-1";
    else if (bal[i] == 0) key += "0";
    else key += "1";
  }

  return key;
}

// ---------------------------------------------------------------
// Format a balance vector for debug printing.
// Example: c(-1, 0, 1) -> "-|.|+"
// ---------------------------------------------------------------
std::string format_balance_cpp(const arma::ivec& bal) {
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

// ---------------------------------------------------------------
// Convert an R list of integer vectors into a C++ vector of uvec.
// Indices in R are 1-based, so they are converted to 0-based.
// ---------------------------------------------------------------
std::vector<arma::uvec> list_to_groups_cpp(const Rcpp::List& lI) {
  std::vector<arma::uvec> groups;
  groups.reserve(lI.size());

  for (int i = 0; i < lI.size(); ++i) {
    IntegerVector idx = lI[i];
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

// ---------------------------------------------------------------
// Expand the active left/right groups of a balance into the original
// variable indices of X.
//
// groups contains the mapping from grouped variables to original
// variable indices.
// ---------------------------------------------------------------
void expand_balance_indices_cpp(const arma::ivec& bal,
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

// ---------------------------------------------------------------
// Evaluate the variance criterion of a grouped balance using the
// covariance matrix M over the original variables.
//
// Returns -Inf if one side is empty.
// ---------------------------------------------------------------
double evaluate_balance_grouped_cpp(const arma::ivec& bal,
                                    const arma::mat& M,
                                    const std::vector<arma::uvec>& groups) {
  arma::uvec iL, iR;
  expand_balance_indices_cpp(bal, groups, iL, iR);

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

// ---------------------------------------------------------------
// Generate all neighbors of a balance:
//
// 1) Remove one active group (-1 or +1) -> 0
// 2) Add one inactive group to the left side -> -1
// 3) Add one inactive group to the right side -> +1
//
// This matches the logic of the R implementation.
// ---------------------------------------------------------------
std::vector<arma::ivec> neighbours_cpp(const arma::ivec& bal) {
  std::vector<arma::ivec> out;
  const arma::uword D = bal.n_elem;

  int n_left = 0;
  int n_right = 0;
  int n_zero = 0;

  for (arma::uword i = 0; i < D; ++i) {
    if (bal[i] == -1) ++n_left;
    else if (bal[i] == 1) ++n_right;
    else ++n_zero;
  }

  out.reserve((n_left + n_right) + 2 * n_zero);

  // Remove one active group
  for (arma::uword i = 0; i < D; ++i) {
    if (bal[i] != 0) {
      arma::ivec cand = bal;
      cand[i] = 0;

      int cand_left = n_left - (bal[i] == -1 ? 1 : 0);
      int cand_right = n_right - (bal[i] ==  1 ? 1 : 0);

      if (cand_left > 0 && cand_right > 0) {
        out.push_back(cand);
      }
    }
  }

  // Add one inactive group to the left side
  for (arma::uword i = 0; i < D; ++i) {
    if (bal[i] == 0) {
      arma::ivec cand = bal;
      cand[i] = -1;
      out.push_back(cand);
    }
  }

  // Add one inactive group to the right side
  for (arma::uword i = 0; i < D; ++i) {
    if (bal[i] == 0) {
      arma::ivec cand = bal;
      cand[i] = 1;
      out.push_back(cand);
    }
  }

  return out;
}

// ---------------------------------------------------------------
// Build the grouped log-composition covariance matrix.
//
// For each group in lI, the grouped log-coordinate is the sum of the
// log-coordinates of the original parts in that group. Therefore,
// the grouped covariance matrix can be obtained by summing the
// corresponding blocks of M.
// ---------------------------------------------------------------
arma::mat grouped_covariance_cpp(const arma::mat& M,
                                 const std::vector<arma::uvec>& groups) {
  const arma::uword G = groups.size();
  arma::mat Mg(G, G, fill::zeros);

  for (arma::uword a = 0; a < G; ++a) {
    for (arma::uword b = a; b < G; ++b) {
      double val = arma::accu(M.submat(groups[a], groups[b]));
      Mg(a, b) = val;
      Mg(b, a) = val;
    }
  }

  return Mg;
}

// ---------------------------------------------------------------
// Compute a guide vector from the first principal component of the
// grouped covariance matrix.
//
// The sign of the eigenvector is arbitrary. This does not matter
// because it only reverses the guide ordering symmetrically.
// ---------------------------------------------------------------
arma::vec guide_scores_from_pc1_cpp(const arma::mat& M,
                                    const std::vector<arma::uvec>& groups) {
  arma::mat Mg = grouped_covariance_cpp(M, groups);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Mg);

  return eigvec.col(eigvec.n_cols - 1);
}

// ---------------------------------------------------------------
// Compute the guide score of a candidate balance as the inner product
// between the discrete balance vector and the continuous guide vector.
// Larger values indicate stronger alignment with the guide direction.
// ---------------------------------------------------------------
double guide_score_cpp(const arma::ivec& bal, const arma::vec& guide) {
  double out = 0.0;
  for (arma::uword j = 0; j < bal.n_elem; ++j) {
    out += static_cast<double>(bal[j]) * guide[j];
  }
  return out;
}

// ---------------------------------------------------------------
// Main tabu search routine for grouped balances.
//
// Inputs:
//   - M: covariance matrix of log(X) on the original variables
//   - lI: list of groups, each group containing 1-based indices
//         of original variables
//   - bal0: initial grouped balance in {-1, 0, 1}
//   - iter: number of iterations
//   - tabu_size: maximum tabu list size
//   - debug: whether to print progress messages
//
// Outputs:
//   - iter_best
//   - tabu_size
//   - steps
//   - dim
//   - variance
//   - balance_raw
// ---------------------------------------------------------------
Rcpp::List partial_pb_tabu_search_legacy_cpp(const arma::mat& M,
                                             const Rcpp::List& lI,
                                             const arma::ivec& bal0,
                                             const int iter,
                                             const int tabu_size,
                                             const bool debug = false) {

  if (M.n_rows != M.n_cols) {
    Rcpp::stop("M must be a square matrix.");
  }
  if (iter < 1) {
    Rcpp::stop("iter must be >= 1.");
  }
  if (tabu_size < 1) {
    Rcpp::stop("tabu_size must be >= 1.");
  }
  if (bal0.n_elem != lI.size()) {
    Rcpp::stop("Length of bal0 must match length(lI).");
  }

  std::vector<arma::uvec> groups = list_to_groups_cpp(lI);

  for (std::size_t g = 0; g < groups.size(); ++g) {
    if (groups[g].n_elem == 0) {
      Rcpp::stop("Groups in lI must be non-empty.");
    }
    if (arma::max(groups[g]) >= M.n_cols) {
      Rcpp::stop("An index in lI exceeds the dimension of M.");
    }
  }

  arma::ivec BAL = bal0;
  arma::vec steps(iter, fill::zeros);

  std::deque<std::string> tabu_queue;
  std::unordered_set<std::string> tabu_set;

  arma::ivec BEST = BAL;
  int BEST_TABU_SIZE = 0;
  double BEST_EV = evaluate_balance_grouped_cpp(BAL, M, groups);
  int BEST_ITER = 0;

  if (!std::isfinite(BEST_EV)) {
    Rcpp::stop("Initial balance is invalid: both sides must be non-empty.");
  }

  if (debug) {
    Rcpp::Rcout << "Current variance: " << BEST_EV << "\n";
    Rcpp::Rcout << "Balance: " << format_balance_cpp(BEST) << "\n\n";
  }

  for (int i = 0; i < iter; ++i) {
    steps[i] = evaluate_balance_grouped_cpp(BAL, M, groups);

    std::string current_key = balance_key_cpp(BAL);

    if (static_cast<int>(tabu_queue.size()) == tabu_size) {
      std::string old_key = tabu_queue.front();
      tabu_queue.pop_front();
      tabu_set.erase(old_key);
    }

    tabu_queue.push_back(current_key);
    tabu_set.insert(current_key);

    std::vector<arma::ivec> BAL_N = neighbours_cpp(BAL);

    if (BAL_N.empty()) {
      break;
    }

    double best_neigh_ev = -std::numeric_limits<double>::infinity();
    int best_neigh_idx = -1;

    for (std::size_t k = 0; k < BAL_N.size(); ++k) {
      std::string key_k = balance_key_cpp(BAL_N[k]);

      if (tabu_set.find(key_k) != tabu_set.end()) {
        continue;
      }

      double ev_k = evaluate_balance_grouped_cpp(BAL_N[k], M, groups);

      if (ev_k > best_neigh_ev) {
        best_neigh_ev = ev_k;
        best_neigh_idx = static_cast<int>(k);
      }
    }

    // All possibilities explored
    if (best_neigh_idx < 0) {
      break;
    }

    BAL = BAL_N[best_neigh_idx];
    double BAL_EV = best_neigh_ev;

    if (debug) {
      Rcpp::Rcout << "Iteration " << (i + 1)
                  << " - current variance: " << BAL_EV << "\n";
      Rcpp::Rcout << "Balance: " << format_balance_cpp(BAL) << "\n\n";
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
