// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "coda.h"

using namespace Rcpp;

namespace {

inline void fill_default_ilr_rows(arma::cube& B,
                                  unsigned int slice,
                                  unsigned int row_offset,
                                  const arma::uvec& idx) {
  const unsigned int m = idx.n_elem;
  if (m <= 1) return;

  for (unsigned int i = 0; i + 1 < m; ++i) {
    const unsigned int m1 = i + 1;
    const unsigned int m2 = i + 2;

    const double l = 1.0 / std::sqrt(static_cast<double>(m1) * m2);
    const double r = -std::sqrt(static_cast<double>(m1) / m2);

    for (unsigned int j = 0; j < m1; ++j) {
      B(row_offset + i, idx(j), slice) = l;
    }
    B(row_offset + i, idx(m1), slice) = r;
  }
}

} // namespace

// [[Rcpp::export]]
arma::cube c_conditional_obasis(const arma::mat& C) {
  const unsigned int n = C.n_cols;
  const unsigned int D = C.n_rows;
  const unsigned int d = D - 1;

  arma::cube B(d, D, n, arma::fill::zeros);

  for (unsigned int k = 0; k < n; ++k) {
    // Build index sets manually: i0 for zeros, i1 for positives
    arma::uvec i0(D);
    arma::uvec i1(D);
    unsigned int nz0 = 0;
    unsigned int nz1 = 0;

    for (unsigned int j = 0; j < D; ++j) {
      const double cj = C(j, k);
      if (cj == 0.0) {
        i0(nz0++) = j;
      } else if (cj > 0.0) {
        i1(nz1++) = j;
      }
    }

    i0.resize(nz0);
    i1.resize(nz1);

    // If there is no zero block, just use the default basis on the positive block
    if (nz0 == 0) {
      B.slice(k) = ilr_basis_default(nz1).t();
      continue;
    }

    // Internal structure of the zero block
    fill_default_ilr_rows(B, k, 0, i0);

    // Contrast between zero block and positive block
    const double scale = std::sqrt(static_cast<double>(nz0 * nz1) / D);
    const double w0 = scale / nz0;
    const double w1 = -scale / nz1;

    for (unsigned int t = 0; t < nz0; ++t) {
      B(nz0 - 1, i0(t), k) = w0;
    }
    for (unsigned int t = 0; t < nz1; ++t) {
      B(nz0 - 1, i1(t), k) = w1;
    }

    // Internal structure of the positive block
    fill_default_ilr_rows(B, k, nz0, i1);
  }

  return B;
}

// [[Rcpp::export]]
arma::cube c_zero_na_conditional_obasis(const arma::mat& tX) {
  const unsigned int D = tX.n_rows;
  const unsigned int n = tX.n_cols;

  // Precompute transposed default ilr bases:
  // basis_cache[k] has dimension (k-1) x k for k >= 2
  std::vector<arma::mat> basis_cache(D + 1);
  for (unsigned int k = 2; k <= D; ++k) {
    basis_cache[k] = ilr_basis_default(k).t();
  }

  arma::cube B(D - 1, D, n, arma::fill::zeros);

  for (unsigned int i = 0; i < n; ++i) {
    // Partition indices into NA / zero / positive
    arma::uvec idx_na(D);
    arma::uvec idx_zero(D);
    arma::uvec idx_pos(D);

    unsigned int n_na = 0;
    unsigned int n_zero = 0;
    unsigned int n_pos = 0;

    for (unsigned int j = 0; j < D; ++j) {
      const double x = tX(j, i);
      if (std::isnan(x)) {
        idx_na(n_na++) = j;
      } else if (x == 0.0) {
        idx_zero(n_zero++) = j;
      } else {
        idx_pos(n_pos++) = j;
      }
    }

    idx_na.resize(n_na);
    idx_zero.resize(n_zero);
    idx_pos.resize(n_pos);

    // If there are no positive parts, return the default basis on all D parts
    if (n_pos == 0) {
      B.slice(i) = basis_cache[D];
      continue;
    }

    unsigned int row_offset = 0;

    // 1) Internal basis for NA block
    if (n_na > 1) {
      const arma::mat& Bna = basis_cache[n_na];
      for (unsigned int r = 0; r < n_na - 1; ++r) {
        for (unsigned int c = 0; c < n_na; ++c) {
          B(row_offset + r, idx_na(c), i) = Bna(r, c);
        }
      }
    }

    // 2) Contrast: NA block vs positive block
    if (n_na > 0) {
      const double w = std::sqrt(static_cast<double>(n_na * n_pos) / (n_na + n_pos));
      const double w_na = w / n_na;
      const double w_pos = -w / n_pos;

      for (unsigned int c = 0; c < n_na; ++c) {
        B(row_offset + (n_na - 1), idx_na(c), i) = w_na;
      }
      for (unsigned int c = 0; c < n_pos; ++c) {
        B(row_offset + (n_na - 1), idx_pos(c), i) = w_pos;
      }

      row_offset += n_na;
    }

    // 3) Internal basis for zero block
    if (n_zero > 1) {
      const arma::mat& Bzero = basis_cache[n_zero];
      for (unsigned int r = 0; r < n_zero - 1; ++r) {
        for (unsigned int c = 0; c < n_zero; ++c) {
          B(row_offset + r, idx_zero(c), i) = Bzero(r, c);
        }
      }
    }

    // 4) Contrast: zero block vs positive block
    if (n_zero > 0) {
      const double w = std::sqrt(static_cast<double>(n_zero * n_pos) / (n_zero + n_pos));
      const double w_zero = w / n_zero;
      const double w_pos = -w / n_pos;

      for (unsigned int c = 0; c < n_zero; ++c) {
        B(row_offset + (n_zero - 1), idx_zero(c), i) = w_zero;
      }
      for (unsigned int c = 0; c < n_pos; ++c) {
        B(row_offset + (n_zero - 1), idx_pos(c), i) = w_pos;
      }

      row_offset += n_zero;
    }

    // 5) Internal basis for positive block
    if (n_pos > 1) {
      const arma::mat& Bpos = basis_cache[n_pos];
      for (unsigned int r = 0; r < n_pos - 1; ++r) {
        for (unsigned int c = 0; c < n_pos; ++c) {
          B(row_offset + r, idx_pos(c), i) = Bpos(r, c);
        }
      }
    }
  }

  return B;
}
