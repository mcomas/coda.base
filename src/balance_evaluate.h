#ifndef BALANCE_EVALUATE_H
#define BALANCE_EVALUATE_H

#include <RcppArmadillo.h>
#include <limits>
#include <map>

class EvaluateBalance {
public:
  EvaluateBalance() = default;
  virtual ~EvaluateBalance() = default;

  virtual double eval(const arma::uvec& L, const arma::uvec& R, int l, int r) {
    return -1.0;
  }
};

class MaximumVariance : public EvaluateBalance {
private:
  arma::mat M;
  arma::vec N;
  double bestScore = -std::numeric_limits<double>::infinity();

public:
  arma::uvec bestL, bestR;

  MaximumVariance() = default;

  MaximumVariance(const std::map<int, arma::uvec>& nodes0, const arma::mat& X) {
    arma::mat S = arma::cov(arma::log(X));

    M = arma::mat(nodes0.size(), nodes0.size());
    N = arma::vec(nodes0.size());

    for (int i = 0; i < static_cast<int>(nodes0.size()); ++i) {
      N(i) = nodes0.at(i).n_elem;
      for (int j = 0; j < static_cast<int>(nodes0.size()); ++j) {
        M(i, j) = arma::accu(S(nodes0.at(i), nodes0.at(j)));
      }
    }
  }

  void init() {
    bestScore = -std::numeric_limits<double>::infinity();
  }

  double eval(const arma::uvec& L, const arma::uvec& R, int l, int r) override {
    double nL = 0.0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(l); ++i) {
      nL += N[L[i]];
    }

    double nR = 0.0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(r); ++i) {
      nR += N[R[i]];
    }

    double variance = 0.0;

    for (int i = 0; i < l; ++i) {
      for (int j = 0; j < l; ++j) {
        variance += (nR / nL) * M(L[i], L[j]);
      }
      for (int j = 0; j < r; ++j) {
        variance -= 2.0 * M(L[i], R[j]);
      }
    }

    for (int i = 0; i < r; ++i) {
      for (int j = 0; j < r; ++j) {
        variance += (nL / nR) * M(R[i], R[j]);
      }
    }

    variance /= (nL + nR);

    if (variance > bestScore) {
      bestScore = variance;
      bestL = arma::uvec(L.head(l));
      bestR = arma::uvec(R.head(r));
    }

    return variance;
  }
};

#endif
