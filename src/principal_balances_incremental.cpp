#include <RcppArmadillo.h>
#include "balance_incremental.h"

// [[Rcpp::depends(RcppArmadillo)]]

void optimise_incremental(BalanceIncremental& balance, const arma::mat& X) {
  MaximumVarianceIncremental ebalance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
}

// [[Rcpp::export]]
arma::mat find_PB2(const arma::mat& X) {
  const int K = X.n_cols;

  arma::mat pb_mat(K, K - 1);
  std::vector<BalanceIncremental> sols;

  BalanceIncremental balance(X.n_cols);
  optimise_incremental(balance, X);
  sols.push_back(balance);

  for (int l = 0; l < K - 1; ++l) {
    double vBestSolution = -std::numeric_limits<double>::infinity();
    int iBestSolution = 0;

    for (unsigned int i = 0; i < sols.size(); ++i) {
      const double v = sols[i].eval();
      if (v > vBestSolution) {
        vBestSolution = v;
        iBestSolution = i;
      }
    }

    pb_mat.col(l) = sols[iBestSolution].getBalance();

    BalanceIncremental top = sols[iBestSolution].top();
    BalanceIncremental left = sols[iBestSolution].left();
    BalanceIncremental right = sols[iBestSolution].right();

    if (top.nodes.size() > 1) {
      optimise_incremental(top, X);
      sols.push_back(top);
    }
    if (left.nodes.size() > 1) {
      optimise_incremental(left, X);
      sols.push_back(left);
    }
    if (right.nodes.size() > 1) {
      optimise_incremental(right, X);
      sols.push_back(right);
    }

    sols[iBestSolution] = sols.back();
    sols.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return pb_mat;
}
