// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "coda.h"

using namespace Rcpp;

namespace {

inline bool is_finite_double(double x) {
  return std::isfinite(x);
}

arma::vec row_means_na(const arma::mat& X) {
  arma::vec out(X.n_rows, arma::fill::zeros);

  for (unsigned int i = 0; i < X.n_rows; ++i) {
    double s = 0.0;
    unsigned int n = 0;
    for (unsigned int j = 0; j < X.n_cols; ++j) {
      const double x = X(i, j);
      if (is_finite_double(x)) {
        s += x;
        ++n;
      }
    }
    out(i) = (n > 0) ? s / n : arma::datum::nan;
  }

  return out;
}

arma::mat pairwise_cov_complete(const arma::mat& X) {
  const unsigned int D = X.n_rows;
  const unsigned int N = X.n_cols;

  arma::mat S(D, D, arma::fill::zeros);

  for (unsigned int a = 0; a < D; ++a) {
    for (unsigned int b = a; b < D; ++b) {
      double ma = 0.0;
      double mb = 0.0;
      unsigned int n = 0;

      for (unsigned int i = 0; i < N; ++i) {
        const double xa = X(a, i);
        const double xb = X(b, i);
        if (is_finite_double(xa) && is_finite_double(xb)) {
          ma += xa;
          mb += xb;
          ++n;
        }
      }

      double covab = 0.0;
      if (n >= 2) {
        ma /= n;
        mb /= n;

        for (unsigned int i = 0; i < N; ++i) {
          const double xa = X(a, i);
          const double xb = X(b, i);
          if (is_finite_double(xa) && is_finite_double(xb)) {
            covab += (xa - ma) * (xb - mb);
          }
        }
        covab /= (n - 1.0);
      }

      S(a, b) = covab;
      S(b, a) = covab;
    }
  }

  return S;
}

arma::mat strict_positive(const arma::mat& X, const double tol = 1e-8) {
  arma::vec eigval;
  arma::mat eigvec;

  arma::mat Xsym = arma::symmatu(X);
  arma::eig_sym(eigval, eigvec, Xsym);

  if (eigval.min() >= tol) {
    return arma::symmatu(X);
  }

  arma::vec new_eig = arma::clamp(eigval, tol, arma::datum::inf);
  return eigvec * arma::diagmat(new_eig) * eigvec.t();
}

arma::mat sweep_matrix(arma::mat A, const unsigned int k_from, const unsigned int k_to) {
  if (k_from > k_to) return A;

  for (unsigned int k = k_from; k <= k_to; ++k) {
    const double D = A(k, k);

    A.row(k) /= D;

    for (unsigned int i = 0; i < A.n_rows; ++i) {
      if (i == k) continue;

      const double B = A(i, k);
      A.row(i) -= B * A.row(k);
      A(i, k) = -B / D;
    }

    A(k, k) = 1.0 / D;
  }

  return A;
}

double dmvnorm_log_single(const arma::vec& x,
                          const arma::vec& mean,
                          const arma::mat& sigma) {
  arma::mat dec;
  if (!arma::chol(dec, sigma)) {
    return R_NegInf;
  }

  // Solve dec.t() * y = x - mean
  arma::vec tmp = arma::solve(arma::trimatl(dec.t()), x - mean);
  const double rss = arma::dot(tmp, tmp);

  return -arma::sum(arma::log(dec.diag())) -
    0.5 * x.n_elem * std::log(2.0 * arma::datum::pi) -
    0.5 * rss;
}

double mLogLik_zr_cpp(const arma::mat& tX,
                      const arma::mat& tDL,
                      const double dl_prop,
                      const arma::vec& mu,
                      const arma::mat& sigma,
                      const arma::cube& Bc) {
  const unsigned int D = tX.n_rows;
  const unsigned int N = tX.n_cols;
  const unsigned int d = D - 1;

  double total = 0.0;

  for (unsigned int i = 0; i < N; ++i) {
    arma::vec lx(D);
    arma::uvec obs_idx(D);
    unsigned int n_obs = 0;
    unsigned int n_na = 0;
    unsigned int n_zero = 0;

    for (unsigned int j = 0; j < D; ++j) {
      const double x = tX(j, i);

      if (std::isnan(x)) {
        lx(j) = arma::datum::nan;
        ++n_na;
      } else if (x == 0.0) {
        lx(j) = std::log(tDL(j, i)) + std::log(dl_prop);
        obs_idx(n_obs++) = j;
        ++n_zero;
      } else {
        lx(j) = std::log(x);
        obs_idx(n_obs++) = j;
      }
    }

    obs_idx.resize(n_obs);

    const unsigned int n_na_zero = n_na + ((n_zero > 0) ? (n_zero - 1) : 0);

    if (n_na_zero >= d) {
      continue;
    }

    const arma::mat Bi = Bc.slice(i);

    if (n_na_zero == 0) {
      const arma::vec h_obs = Bi * lx;
      const arma::vec mu_obs = Bi * mu;
      const arma::mat sigma_obs = Bi * sigma * Bi.t();

      total += dmvnorm_log_single(h_obs, mu_obs, sigma_obs);
    } else {
      const arma::mat B2 = Bi.rows(n_na_zero, d - 1);
      const arma::vec lx_obs = lx.elem(obs_idx);
      const arma::mat B2_obs = B2.cols(obs_idx);

      const arma::vec h_obs = B2_obs * lx_obs;
      const arma::vec mu_obs = B2 * mu;
      const arma::mat sigma_obs = B2 * sigma * B2.t();

      total += dmvnorm_log_single(h_obs, mu_obs, sigma_obs);
    }
  }

  return total;
}

} // namespace


// [[Rcpp::export]]
SEXP c_coda_replacement(const arma::mat& tX,
                        const arma::mat& tDL,
                        const double dl_prop = 0.65,
                        const double eps = 1e-4,
                        const bool parameters = false,
                        const bool debug = false,
                        const unsigned int maxit = 500) {
  const unsigned int D = tX.n_rows;
  const unsigned int N = tX.n_cols;
  const unsigned int d = D - 1;

  if (tDL.n_rows != D || tDL.n_cols != N) {
    stop("'tDL' must have the same dimensions as 'tX'.");
  }

  arma::mat lX(D, N, arma::fill::zeros);
  arma::uvec n_na(N, arma::fill::zeros);
  arma::uvec n_zero(N, arma::fill::zeros);
  arma::uvec n_na_zero(N, arma::fill::zeros);

  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < D; ++j) {
      const double x = tX(j, i);

      if (std::isnan(x)) {
        lX(j, i) = arma::datum::nan;
        ++n_na(i);
      } else if (x == 0.0) {
        lX(j, i) = std::log(tDL(j, i)) + std::log(dl_prop);
        ++n_zero(i);
      } else {
        lX(j, i) = std::log(x);
      }
    }

    n_na_zero(i) = n_na(i) + ((n_zero(i) > 0) ? (n_zero(i) - 1) : 0);
  }

  arma::mat G = arma::eye(D, D) - (1.0 / D) * arma::ones(D, D);

  arma::vec MU = G * row_means_na(lX);

  arma::mat S0 = pairwise_cov_complete(lX);
  S0 *= static_cast<double>(N - 1) / N;
  arma::mat SIGMA = strict_positive(arma::symmatu(G * S0 * G));

  arma::cube Bc = c_zero_na_conditional_obasis(tX);

  std::vector<arma::mat> pinvBi(N);
  for (unsigned int i = 0; i < N; ++i) {
    pinvBi[i] = arma::pinv(Bc.slice(i));
  }

  arma::mat M1(D, N, arma::fill::zeros);
  arma::mat M2(D, D, arma::fill::zeros);

  bool CONT = true;
  unsigned int it = 0;
  double ll_old = R_NegInf;

  while (CONT && it < maxit) {
    ++it;
    M2.zeros();

    for (unsigned int i = 0; i < N; ++i) {
      arma::vec m1(D, arma::fill::zeros);
      arma::mat m2(D, D, arma::fill::zeros);

      if (n_na_zero(i) >= d) {
        m1 = MU;
        m2 = MU * MU.t() + SIGMA;

      } else if (n_na_zero(i) == 0) {
        arma::vec lx = lX.col(i);
        const double m = arma::mean(lx);
        m1 = lx - m;
        m2 = m1 * m1.t();

      } else {
        const unsigned int q = n_na_zero(i);

        arma::uvec obs_idx(D);
        unsigned int n_obs = 0;
        for (unsigned int j = 0; j < D; ++j) {
          if (!std::isnan(lX(j, i))) {
            obs_idx(n_obs++) = j;
          }
        }
        obs_idx.resize(n_obs);

        const arma::mat Bi = Bc.slice(i);
        const arma::mat& P = pinvBi[i];

        const arma::mat B2 = Bi.rows(q, d - 1);
        const arma::vec lx_col = lX.col(i);
        const arma::vec lx_obs = lx_col.elem(obs_idx);
        const arma::mat B2_obs = B2.cols(obs_idx);
        const arma::vec h2 = B2_obs * lx_obs;

        const arma::vec mu_i = Bi * MU;
        const arma::mat sigma_i = Bi * SIGMA * Bi.t();

        arma::mat aug(d + 1, d + 1, arma::fill::zeros);
        aug.submat(0, 0, d - 1, d - 1) = sigma_i;
        aug.submat(0, d, d - 1, d) = mu_i;
        aug.submat(d, 0, d, d - 1) = mu_i.t();
        aug(d, d) = 1.0;

        arma::mat swept = sweep_matrix(aug, q, d - 1);

        arma::vec rhs((d - q) + 1, arma::fill::ones);
        rhs.head(d - q) = h2;

        arma::mat A = swept.submat(q, 0, d, q - 1);
        arma::vec h1_m1 = A.t() * rhs;

        arma::mat h1_m2 = h1_m1 * h1_m1.t() + swept.submat(0, 0, q - 1, q - 1);

        arma::mat upper = arma::join_rows(h1_m2, h1_m1 * h2.t());
        arma::mat lower = arma::join_rows(h2 * h1_m1.t(), h2 * h2.t());
        arma::mat H = arma::join_cols(upper, lower);

        arma::vec z = arma::join_cols(h1_m1, h2);

        m1 = P * z;
        m2 = P * H * P.t();
      }

      M1.col(i) = m1;
      M2 += m2;
    }

    arma::vec MU_new = arma::mean(M1, 1);
    arma::mat SIGMA_new = strict_positive(arma::symmatu(M2 / N - MU_new * MU_new.t()));

    const double ll_new = mLogLik_zr_cpp(tX, tDL, dl_prop, MU_new, SIGMA_new, Bc);

    if (debug) {
      Rcpp::Rcout << "Iteration: " << it << ", LogLik: " << ll_new << "\n";
    }

    // arma::vec denom = arma::clamp(arma::abs(MU), eps, arma::datum::inf);
    // const double rel_change = arma::max(arma::abs(MU_new - MU) / denom);

    MU = MU_new;
    SIGMA = SIGMA_new;

    if (it > 1 && std::abs(ll_new - ll_old) < eps * (1.0 + std::abs(ll_old))) {
      CONT = false;
    }

    ll_old = ll_new;

    // if (rel_change < eps) {
    //   CONT = false;
    // }
  }
  if (it == maxit) {
    Rcpp::Rcout << "Maximum number of iterations reached (" << maxit << ").\n";
  }
  if (parameters) {
    return List::create(
      _["clr_mu"] = MU,
      _["clr_sigma"] = SIGMA,
      _["clr_h"] = M1.t()
    );
  }

  arma::mat Ximp = arma::exp(M1);
  for (unsigned int i = 0; i < N; ++i) {
    const double s = arma::accu(Ximp.col(i));
    if (s > 0) {
      Ximp.col(i) /= s;
    }
  }

  return wrap(Ximp.t());
}
