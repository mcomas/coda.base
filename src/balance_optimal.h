#ifndef BALANCE_OPTIMAL_H
#define BALANCE_OPTIMAL_H

#include <RcppArmadillo.h>
#include <vector>

template <class EB>
void f(int mu, int nu, int sigma,
       arma::uvec& I, arma::uvec& A,
       std::vector<arma::uvec>& P, int* p, EB& ebalance);

template <class EB>
void b(int mu, int nu, int sigma,
       arma::uvec& I, arma::uvec& A,
       std::vector<arma::uvec>& P, int* p, EB& ebalance);

template <class EB>
inline void visit(arma::uvec& /*I*/, arma::uvec& /*A*/,
                  std::vector<arma::uvec>& P, int* p, EB& ebalance) {
  ebalance.eval(P[1], P[2], p[1], p[2]);
}

inline void update(unsigned int e, unsigned int s1, unsigned int s2,
                   arma::uvec& I, arma::uvec& /*A*/,
                   std::vector<arma::uvec>& P, int* p) {
  // Remove element e from set s1
  p[s1] = p[s1] - 1;
  P[s1][I[e]] = P[s1][p[s1]];
  I[P[s1][I[e]]] = I[e];

  // Add element e to set s2
  P[s2][p[s2]] = e;
  I[e] = p[s2];
  p[s2] = p[s2] + 1;
}

template <class EB>
void f(int mu, int nu, int sigma,
       arma::uvec& I, arma::uvec& A,
       std::vector<arma::uvec>& P, int* p, EB& ebalance) {
  unsigned int e, s1, s2;

  if (mu == 2) {
    visit<EB>(I, A, P, p, ebalance);
  } else {
    f<EB>(mu - 1, nu - 1, (mu + sigma) % 2, I, A, P, p, ebalance);
  }

  if (nu == mu) {
    e = mu - 1;
    s1 = 0;
    s2 = mu - 1;
    A[e] = s2;
    update(e - 1, s1, s2, I, A, P, p);

    visit<EB>(I, A, P, p, ebalance);

    while (A[nu] > 0) {
      e = nu;
      s1 = A[nu];
      s2 = A[nu] - 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);

      visit<EB>(I, A, P, p, ebalance);
    }
  } else if (nu > mu) {
    if ((mu + sigma) % 2) {
      e = nu - 1;
      s1 = 0;
      s2 = mu - 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
    } else {
      e = mu - 1;
      s1 = 0;
      s2 = mu - 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
    }

    if ((A[nu] + sigma) % 2) {
      b<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
    } else {
      f<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
    }

    while (A[nu] > 0) {
      e = nu;
      s1 = A[nu];
      s2 = A[nu] - 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);

      if ((A[nu] + sigma) % 2) {
        b<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
      } else {
        f<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
      }
    }
  }
}

template <class EB>
void b(int mu, int nu, int sigma,
       arma::uvec& I, arma::uvec& A,
       std::vector<arma::uvec>& P, int* p, EB& ebalance) {
  unsigned int e, s1, s2;

  if (nu == mu) {
    while (A[nu] < static_cast<unsigned int>(mu - 1)) {
      visit<EB>(I, A, P, p, ebalance);

      e = nu;
      s1 = A[nu];
      s2 = A[nu] + 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
    }

    visit<EB>(I, A, P, p, ebalance);

    e = mu - 1;
    s1 = A[mu - 1];
    s2 = 0;
    A[e] = s2;
    update(e - 1, s1, s2, I, A, P, p);
  } else if (nu > mu) {
    if ((A[nu] + sigma) % 2) {
      f<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
    } else {
      b<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
    }

    while (A[nu] < static_cast<unsigned int>(mu - 1)) {
      e = nu;
      s1 = A[nu];
      s2 = A[nu] + 1;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);

      if ((A[nu] + sigma) % 2) {
        f<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
      } else {
        b<EB>(mu, nu - 1, 0, I, A, P, p, ebalance);
      }
    }

    if ((mu + sigma) % 2) {
      e = nu - 1;
      s1 = A[e];
      s2 = 0;
      A[e] = s2;
      update(e - 1, s1, s2, I, A, P, p);
    } else {
      e = mu - 1;
      s1 = A[e];
      s2 = 0;
      A[e] = 0;
      update(e - 1, s1, s2, I, A, P, p);
    }
  }

  if (mu == 2) {
    visit<EB>(I, A, P, p, ebalance);
  } else {
    b<EB>(mu - 1, nu - 1, (mu + sigma) % 2, I, A, P, p, ebalance);
  }
}

#endif
