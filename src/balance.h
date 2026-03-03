#ifndef BALANCE_H
#define BALANCE_H

#include <RcppArmadillo.h>
#include "balance_evaluate.h"
#include "balance_optimal.h"

template <class EB>
class Balance {
public:
  int D;

  arma::uvec L, R;
  unsigned int L_length, R_length;

  std::map<int, arma::uvec> nodes;
  int n_nodes;

  EB ebalance = EB();
  double ebalance_value = 0.0;

  Balance(int D0) {
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

  Balance(int D0, const std::map<int, arma::uvec>& nodes0) {
    D = D0;
    nodes = nodes0;
    n_nodes = nodes0.size();

    L = arma::uvec(n_nodes);
    L_length = 0;
    R = arma::uvec(n_nodes);
    R_length = 0;
  }

  void setEvaluator(const EB& ebalance0) {
    ebalance = ebalance0;
  }

  void set(const arma::uvec& uL, const arma::uvec& uR) {
    L_length = uL.n_elem;
    L.head(L_length) = uL;

    R_length = uR.n_elem;
    R.head(R_length) = uR;

    ebalance_value = ebalance.eval(L, R, L_length, R_length);
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

    f<EB>(3, N, 0, I, A, P, p, ebalance);

    set(ebalance.bestL, ebalance.bestR);
  }

  void setWithLogContrast(arma::vec V) {
    const int imin = arma::index_min(V);
    const int imax = arma::index_max(V);

    V(imin) = 0;
    V(imax) = 0;

    arma::uvec ord = arma::sort_index(arma::abs(V), "descend");
    arma::uvec uL(ord.n_elem), uR(ord.n_elem);
    uL[0] = imin;
    uR[0] = imax;

    unsigned int l = 1;
    unsigned int r = 1;

    ebalance.eval(uL, uR, l, r);

    for (int i = 0; i < n_nodes - 2; ++i) {
      if (V(ord[i]) < 0) {
        uL(l++) = ord[i];
      } else {
        uR(r++) = ord[i];
      }

      ebalance.eval(uL, uR, l, r);
    }

    set(ebalance.bestL, ebalance.bestR);
  }

  void setWithLogContrastForceBranch(arma::vec V, const arma::uvec& forced) {
    double vforced = V(forced[0]);
    V(forced[0]) = 0;

    for (unsigned int i = 1; i < forced.n_elem; ++i) {
      vforced += V(forced[i]);
      V(forced[i]) = 0;
    }

    const int imax = arma::index_max(arma::abs(V));
    arma::uvec ord;

    if (V(imax) > 0) {
      ord = arma::sort_index(V, "descend");
    } else {
      ord = arma::sort_index(V, "ascend");
    }

    unsigned int pos = 0;
    arma::uvec uR(ord.n_elem);

    while (pos < ord.n_elem && V(ord[pos]) != 0) {
      uR[pos] = ord[pos];
      ++pos;
      ebalance.eval(forced, uR, forced.n_elem, pos);
    }

    set(ebalance.bestL, ebalance.bestR);
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

  Balance<EB> top() const {
    arma::uvec O = arma::zeros<arma::uvec>(nodes.size());
    O(L.head(L_length)).fill(1);
    O(R.head(R_length)).fill(1);

    arma::uvec uV = arma::find(O == 1);
    arma::uvec uI = arma::find(O == 0);

    const int nI = uI.n_elem;
    std::map<int, arma::uvec> node0;

    for (int i = 0; i < nI; ++i) {
      node0[i] = nodes.at(uI[i]);
    }

    int nV = 0;
    for (unsigned int i = 0; i < L_length; ++i) {
      nV += nodes.at(L[i]).n_elem;
    }
    for (unsigned int i = 0; i < R_length; ++i) {
      nV += nodes.at(R[i]).n_elem;
    }

    arma::uvec V(nV);
    int k = 0;

    for (unsigned int i = 0; i < L_length; ++i) {
      for (unsigned int j = 0; j < nodes.at(L[i]).n_elem; ++j, ++k) {
        V(k) = nodes.at(L[i])[j];
      }
    }
    for (unsigned int i = 0; i < R_length; ++i) {
      for (unsigned int j = 0; j < nodes.at(R[i]).n_elem; ++j, ++k) {
        V(k) = nodes.at(R[i])[j];
      }
    }

    node0[nI] = V;

    return Balance<EB>(D, node0);
  }

  Balance<EB> left() const {
    arma::uvec uL = L.head(L_length);
    std::map<int, arma::uvec> node0;

    for (unsigned int i = 0; i < uL.n_elem; ++i) {
      node0[i] = nodes.at(uL[i]);
    }

    return Balance<EB>(D, node0);
  }

  Balance<EB> right() const {
    arma::uvec uR = R.head(R_length);
    std::map<int, arma::uvec> node0;

    for (unsigned int i = 0; i < uR.n_elem; ++i) {
      node0[i] = nodes.at(uR[i]);
    }

    return Balance<EB>(D, node0);
  }

  void print() const {
    Rcpp::Rcout << "Elements: ";
    for (unsigned int i = 0; i < n_nodes; ++i) {
      Rcpp::Rcout << "{";
      for (unsigned int j = 0; j < nodes.at(i).n_elem; ++j) {
        Rcpp::Rcout << " " << nodes.at(i)[j];
      }
      Rcpp::Rcout << " } ";
    }

    Rcpp::Rcout << "L:";
    for (unsigned int i = 0; i < L_length; ++i) {
      Rcpp::Rcout << " " << L[i];
    }

    Rcpp::Rcout << ", R:";
    for (unsigned int i = 0; i < R_length; ++i) {
      Rcpp::Rcout << " " << R[i];
    }
    Rcpp::Rcout << std::endl;

    if (L_length > 0 && R_length > 0) {
      Rcpp::Rcout << getBalance().t()
                  << "Value = " << ebalance_value << std::endl;
    } else {
      Rcpp::Rcout << "Balance is not defined" << std::endl;
    }
  }
};

#endif
