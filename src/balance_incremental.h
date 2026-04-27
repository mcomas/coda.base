#ifndef BALANCE_INCREMENTAL_H
#define BALANCE_INCREMENTAL_H

#include <RcppArmadillo.h>
#include "balance.h"

class MaximumVarianceIncremental {
private:
  arma::mat M;
  arma::vec N;
  arma::ivec set_id;

  double nL = 0.0;
  double nR = 0.0;
  double sumLL = 0.0;
  double sumRR = 0.0;
  double sumLR = 0.0;
  double bestScore = -std::numeric_limits<double>::infinity();

  double cross_sum(unsigned int e, int set) const {
    double out = 0.0;
    for (arma::uword i = 0; i < set_id.n_elem; ++i) {
      if (set_id[i] == set) {
        out += M(e, i);
      }
    }
    return out;
  }

  void remove_from_left(unsigned int e) {
    const double crossL = cross_sum(e, 1);
    const double crossR = cross_sum(e, 2);

    sumLL -= 2.0 * crossL - M(e, e);
    sumLR -= crossR;
    nL -= N[e];
    set_id[e] = 0;
  }

  void remove_from_right(unsigned int e) {
    const double crossR = cross_sum(e, 2);
    const double crossL = cross_sum(e, 1);

    sumRR -= 2.0 * crossR - M(e, e);
    sumLR -= crossL;
    nR -= N[e];
    set_id[e] = 0;
  }

  void add_to_left(unsigned int e) {
    const double crossL = cross_sum(e, 1);
    const double crossR = cross_sum(e, 2);

    sumLL += 2.0 * crossL + M(e, e);
    sumLR += crossR;
    nL += N[e];
    set_id[e] = 1;
  }

  void add_to_right(unsigned int e) {
    const double crossR = cross_sum(e, 2);
    const double crossL = cross_sum(e, 1);

    sumRR += 2.0 * crossR + M(e, e);
    sumLR += crossL;
    nR += N[e];
    set_id[e] = 2;
  }

public:
  arma::uvec bestL, bestR;

  MaximumVarianceIncremental() = default;

  MaximumVarianceIncremental(const std::map<int, arma::uvec>& nodes0,
                             const arma::mat& X) {
    arma::mat S = arma::cov(arma::log(X));

    M = arma::mat(nodes0.size(), nodes0.size());
    N = arma::vec(nodes0.size());
    set_id = arma::ivec(nodes0.size(), arma::fill::zeros);

    for (int i = 0; i < static_cast<int>(nodes0.size()); ++i) {
      N(i) = nodes0.at(i).n_elem;
      for (int j = 0; j < static_cast<int>(nodes0.size()); ++j) {
        M(i, j) = arma::accu(S(nodes0.at(i), nodes0.at(j)));
      }
    }
  }

  void initialise(const std::vector<arma::uvec>& P, const int* p) {
    set_id.zeros();
    nL = 0.0;
    nR = 0.0;
    sumLL = 0.0;
    sumRR = 0.0;
    sumLR = 0.0;
    bestScore = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < p[1]; ++i) {
      add_to_left(P[1][i]);
    }
    for (int i = 0; i < p[2]; ++i) {
      add_to_right(P[2][i]);
    }
  }

  void move(unsigned int e, unsigned int s1, unsigned int s2) {
    if (s1 == s2) {
      return;
    }

    if (s1 == 1) {
      remove_from_left(e);
    } else if (s1 == 2) {
      remove_from_right(e);
    } else {
      set_id[e] = 0;
    }

    if (s2 == 1) {
      add_to_left(e);
    } else if (s2 == 2) {
      add_to_right(e);
    } else {
      set_id[e] = 0;
    }
  }

  double eval(const arma::uvec& L, const arma::uvec& R, int l, int r) {
    if (l == 0 || r == 0 || nL <= 0.0 || nR <= 0.0) {
      return -std::numeric_limits<double>::infinity();
    }

    const double variance =
      ((nR / nL) * sumLL + (nL / nR) * sumRR - 2.0 * sumLR) / (nL + nR);

    if (variance > bestScore) {
      bestScore = variance;
      bestL = arma::uvec(L.head(l));
      bestR = arma::uvec(R.head(r));
    }

    return variance;
  }

  double best() const {
    return bestScore;
  }
};

template <class EB>
void f_incremental(int mu, int nu, int sigma,
                   arma::uvec& I, arma::uvec& A,
                   std::vector<arma::uvec>& P, int* p, EB& ebalance);

template <class EB>
void b_incremental(int mu, int nu, int sigma,
                   arma::uvec& I, arma::uvec& A,
                   std::vector<arma::uvec>& P, int* p, EB& ebalance);

template <class EB>
inline void visit_incremental(arma::uvec& /*I*/, arma::uvec& /*A*/,
                              std::vector<arma::uvec>& P, int* p,
                              EB& ebalance) {
  ebalance.eval(P[1], P[2], p[1], p[2]);
}

template <class EB>
inline void update_incremental(unsigned int e, unsigned int s1, unsigned int s2,
                               arma::uvec& I, arma::uvec& /*A*/,
                               std::vector<arma::uvec>& P, int* p,
                               EB& ebalance) {
  p[s1] = p[s1] - 1;
  P[s1][I[e]] = P[s1][p[s1]];
  I[P[s1][I[e]]] = I[e];

  P[s2][p[s2]] = e;
  I[e] = p[s2];
  p[s2] = p[s2] + 1;

  ebalance.move(e, s1, s2);
}

template <class EB>
void f_incremental(int mu, int nu, int sigma,
                   arma::uvec& I, arma::uvec& A,
                   std::vector<arma::uvec>& P, int* p, EB& ebalance) {
  unsigned int e, s1, s2;

  if (mu == 2) {
    visit_incremental<EB>(I, A, P, p, ebalance);
  } else {
    f_incremental<EB>(mu - 1, nu - 1, (mu + sigma) % 2, I, A, P, p, ebalance);
  }

  if (nu == mu) {
    e = mu - 1;
    s1 = 0;
    s2 = mu - 1;
    A[e] = s2;
    update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);

    visit_incremental<EB>(I, A, P, p, ebalance);

    while (A[nu] > 0) {
      e = nu;
      s1 = A[nu];
      s2 = A[nu] - 1;
      A[e] = s2;
      update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);

      visit_incremental<EB>(I, A, P, p, ebalance);
    }
  } else if (nu > mu) {
    if ((mu + sigma) % 2) {
      e = nu - 1;
      s1 = 0;
      s2 = mu - 1;
      A[e] = s2;
      update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);
    } else {
      e = mu - 1;
      s1 = 0;
      s2 = mu - 1;
      A[e] = s2;
      update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);
    }

    if ((A[nu] + sigma) % 2) {
      b_incremental<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
    } else {
      f_incremental<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
    }

    while (A[nu] > 0) {
      e = nu;
      s1 = A[nu];
      s2 = A[nu] - 1;
      A[e] = s2;
      update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);

      if ((A[nu] + sigma) % 2) {
        b_incremental<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
      } else {
        f_incremental<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
      }
    }
  }
}

template <class EB>
void b_incremental(int mu, int nu, int sigma,
                   arma::uvec& I, arma::uvec& A,
                   std::vector<arma::uvec>& P, int* p, EB& ebalance) {
  unsigned int e, s1, s2;

  if (nu == mu) {
    while (A[nu] < static_cast<unsigned int>(mu - 1)) {
      visit_incremental<EB>(I, A, P, p, ebalance);

      e = nu;
      s1 = A[nu];
      s2 = A[nu] + 1;
      A[e] = s2;
      update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);
    }

    visit_incremental<EB>(I, A, P, p, ebalance);

    e = mu - 1;
    s1 = A[mu - 1];
    s2 = 0;
    A[e] = s2;
    update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);
  } else if (nu > mu) {
    if ((A[nu] + sigma) % 2) {
      f_incremental<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
    } else {
      b_incremental<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
    }

    while (A[nu] < static_cast<unsigned int>(mu - 1)) {
      e = nu;
      s1 = A[nu];
      s2 = A[nu] + 1;
      A[e] = s2;
      update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);

      if ((A[nu] + sigma) % 2) {
        f_incremental<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
      } else {
        b_incremental<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
      }
    }

    if ((mu + sigma) % 2) {
      e = nu - 1;
      s1 = A[e];
      s2 = 0;
      A[e] = s2;
      update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);
    } else {
      e = mu - 1;
      s1 = A[e];
      s2 = 0;
      A[e] = 0;
      update_incremental(e - 1, s1, s2, I, A, P, p, ebalance);
    }
  }

  if (mu == 2) {
    visit_incremental<EB>(I, A, P, p, ebalance);
  } else {
    b_incremental<EB>(mu - 1, nu - 1, (mu + sigma) % 2, I, A, P, p, ebalance);
  }
}

class BalanceIncremental {
public:
  int D;

  arma::uvec L, R;
  unsigned int L_length, R_length;

  std::map<int, arma::uvec> nodes;
  int n_nodes;

  MaximumVarianceIncremental ebalance = MaximumVarianceIncremental();
  double ebalance_value = 0.0;

  BalanceIncremental(int D0) {
    D = D0;
    for (int i = 0; i < D; ++i) {
      nodes[i] = arma::uvec(1);
      nodes[i][0] = i;
    }
    n_nodes = D;

    L = arma::uvec(n_nodes);
    L_length = 0;
    R = arma::uvec(n_nodes);
    R_length = 0;
  }

  BalanceIncremental(int D0, const std::map<int, arma::uvec>& nodes0) {
    D = D0;
    nodes = nodes0;
    n_nodes = nodes0.size();

    L = arma::uvec(n_nodes);
    L_length = 0;
    R = arma::uvec(n_nodes);
    R_length = 0;
  }

  explicit BalanceIncremental(const Balance<MaximumVariance>& balance) {
    D = balance.D;
    nodes = balance.nodes;
    n_nodes = balance.n_nodes;

    L = arma::uvec(n_nodes);
    L_length = balance.L_length;
    L.head(L_length) = balance.L.head(balance.L_length);

    R = arma::uvec(n_nodes);
    R_length = balance.R_length;
    R.head(R_length) = balance.R.head(balance.R_length);

    ebalance_value = balance.ebalance_value;
  }

  void setEvaluator(const MaximumVarianceIncremental& ebalance0) {
    ebalance = ebalance0;
  }

  void set(const arma::uvec& uL, const arma::uvec& uR) {
    L_length = uL.n_elem;
    L.head(L_length) = uL;

    R_length = uR.n_elem;
    R.head(R_length) = uR;
  }

  void setWithExhaustiveSearch() {
    const int N = nodes.size();

    std::vector<arma::uvec> P(3);
    P[0] = arma::uvec(N);
    P[1] = arma::uvec(N);
    P[2] = arma::uvec(N);

    int p[3];
    p[0] = N - 2;
    p[1] = 1;
    p[2] = 1;

    for (int i = 0; i < N - 2; ++i) {
      P[0][i] = i;
    }
    P[1][0] = N - 2;
    P[2][0] = N - 1;

    arma::uvec I(N);
    for (int i = 0; i < N - 2; ++i) {
      I(i) = i;
    }
    I(N - 2) = 0;
    I(N - 1) = 0;

    arma::uvec A(N + 1);
    A.fill(0);
    A[N - 1] = 1;
    A[N] = 2;

    ebalance.initialise(P, p);
    f_incremental<MaximumVarianceIncremental>(3, N, 0, I, A, P, p, ebalance);

    set(ebalance.bestL, ebalance.bestR);
    ebalance_value = ebalance.best();
  }

  arma::vec getBalance() const {
    double nL = 0.0;
    double nR = 0.0;

    for (unsigned int i = 0; i < L_length; ++i) {
      nL += nodes.at(L[i]).n_elem;
    }
    for (unsigned int i = 0; i < R_length; ++i) {
      nR += nodes.at(R[i]).n_elem;
    }

    arma::vec b = arma::zeros(D);
    const double l_v = -1.0 / nL * std::sqrt(nL * nR / (nL + nR));
    const double r_v = +1.0 / nR * std::sqrt(nL * nR / (nL + nR));

    for (unsigned int i = 0; i < L_length; ++i) {
      b(nodes.at(L[i])).fill(l_v);
    }
    for (unsigned int i = 0; i < R_length; ++i) {
      b(nodes.at(R[i])).fill(r_v);
    }

    return b;
  }

  double eval() const {
    return ebalance_value;
  }

  BalanceIncremental top() const {
    Balance<MaximumVariance> tmp(D, nodes);
    tmp.L_length = L_length;
    tmp.R_length = R_length;
    tmp.L.head(L_length) = L.head(L_length);
    tmp.R.head(R_length) = R.head(R_length);
    return BalanceIncremental(tmp.top());
  }

  BalanceIncremental left() const {
    Balance<MaximumVariance> tmp(D, nodes);
    tmp.L_length = L_length;
    tmp.L.head(L_length) = L.head(L_length);
    return BalanceIncremental(tmp.left());
  }

  BalanceIncremental right() const {
    Balance<MaximumVariance> tmp(D, nodes);
    tmp.R_length = R_length;
    tmp.R.head(R_length) = R.head(R_length);
    return BalanceIncremental(tmp.right());
  }
};

#endif
