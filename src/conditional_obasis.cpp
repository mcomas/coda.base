// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "coda.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::cube c_conditional_obasis(const arma::mat& C) {
  const unsigned int n = C.n_cols;
  const unsigned int D = C.n_rows;
  const unsigned int d = D - 1;

  arma::uvec d0(n, arma::fill::zeros);
  arma::uvec d1(n, arma::fill::zeros);
  arma::uvec D1(n, arma::fill::zeros);

  arma::umat I0(D, n, arma::fill::zeros);
  arma::umat I1(D, n, arma::fill::zeros);

  for (unsigned int k = 0; k < n; ++k) {
    arma::uvec i0 = arma::find(C.col(k) == 0);
    d0(k) = i0.n_elem;
    if (d0(k) > 0) {
      I0.col(k).head(d0(k)) = i0;
    }

    arma::uvec i1 = arma::find(C.col(k) > 0);
    D1(k) = i1.n_elem;
    d1(k) = (D1(k) > 0) ? D1(k) - 1 : 0;
    if (D1(k) > 0) {
      I1.col(k).head(D1(k)) = i1;
    }
  }

  arma::cube B(d, D, n, arma::fill::zeros);

  for (unsigned int k = 0; k < n; ++k) {
    if (d0(k) == 0) {
      B.slice(k) = ilr_basis_default(D1(k)).t();
      continue;
    }

    if (d0(k) > 1) {
      arma::uvec i0 = I0.col(k).head(d0(k));
      for (unsigned int i = 0; i + 1 < d0(k); ++i) {
        const unsigned int m1 = i + 1;
        const unsigned int m2 = i + 2;
        const double l = 1.0 / std::sqrt(static_cast<double>(m1) * m2);
        const double r = -std::sqrt(static_cast<double>(m1) / m2);

        for (unsigned int j = 0; j < m1; ++j) {
          B(i, i0(j), k) = l;
        }
        B(i, i0(m1), k) = r;
      }
    }

    for (unsigned int j = 0; j < D; ++j) {
      if (C(j, k) == 0) {
        B(d0(k) - 1, j, k) =
          1.0 / d0(k) * std::sqrt(static_cast<double>(d0(k) * D1(k)) / D);
      } else {
        B(d0(k) - 1, j, k) =
          -1.0 / D1(k) * std::sqrt(static_cast<double>(d0(k) * D1(k)) / D);
      }
    }

    if (D1(k) > 1) {
      arma::uvec i1 = I1.col(k).head(D1(k));
      for (unsigned int i = 0; i + 1 < D1(k); ++i) {
        const unsigned int m1 = i + 1;
        const unsigned int m2 = i + 2;
        const double l = 1.0 / std::sqrt(static_cast<double>(m1) * m2);
        const double r = -std::sqrt(static_cast<double>(m1) / m2);

        for (unsigned int j = 0; j < m1; ++j) {
          B(d0(k) + i, i1(j), k) = l;
        }
        B(d0(k) + i, i1(m1), k) = r;
      }
    }
  }

  return B;
}
