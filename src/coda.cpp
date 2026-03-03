#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "coda.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat pinv(const arma::mat& X) {
  return arma::pinv(X);
}

// [[Rcpp::export]]
arma::mat c_variation_array(const arma::mat& X,
                            bool include_means,
                            bool ml_covariance) {
  const unsigned int K = X.n_cols;
  const unsigned int n = X.n_rows;

  arma::mat lX = arma::log(X);
  arma::mat varray(K, K, arma::fill::zeros);
  arma::mat Xcov = arma::cov(lX);

  if (ml_covariance) {
    Xcov = (static_cast<double>(n - 1) / n) * Xcov;
  }

  if (include_means) {
    arma::rowvec Xmeans = arma::mean(lX, 0);
    for (unsigned int i = 0; i < K; ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        varray(i, j) = Xmeans(i) - Xmeans(j);
        varray(j, i) = Xcov(i, i) + Xcov(j, j) - 2.0 * Xcov(i, j);
      }
    }
  } else {
    for (unsigned int i = 0; i < K; ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        const double v = Xcov(i, i) + Xcov(j, j) - 2.0 * Xcov(i, j);
        varray(i, j) = v;
        varray(j, i) = v;
      }
    }
  }

  return varray;
}

// [[Rcpp::export]]
arma::mat alr_basis_default(unsigned int dim) {
  arma::mat B(dim, dim - 1, arma::fill::zeros);

  for (unsigned int i = 0; i < dim - 1; ++i) {
    B(i, i) = 1.0;
  }
  B.row(dim - 1).fill(-1.0);

  return B;
}

// [[Rcpp::export]]
arma::mat clr_basis_default(unsigned int dim) {
  arma::mat B(dim, dim);
  const double dim_inv = -1.0 / static_cast<double>(dim);

  for (unsigned int i = 0; i < dim; ++i) {
    for (unsigned int j = 0; j < dim; ++j) {
      B(i, j) = (i == j) ? 1.0 + dim_inv : dim_inv;
    }
  }

  return B;
}

// [[Rcpp::export]]
arma::mat ilr_basis_default(unsigned int dim) {
  arma::mat B(dim, dim - 1, arma::fill::zeros);

  for (unsigned int i = 0; i < dim - 1; ++i) {
    const unsigned int I1 = i + 1;
    const unsigned int I2 = i + 2;

    const double l = 1.0 / std::sqrt(static_cast<double>(I1) * I2);
    const double r = -std::sqrt(static_cast<double>(I1) / I2);

    double* col = B.colptr(i);

    for (unsigned int j = 0; j < I1; ++j) {
      col[j] = l;
    }
    col[I1] = r;
  }

  return B;
}

// [[Rcpp::export]]
arma::mat ilr_basis_simplex(unsigned int dim) {
  arma::mat B(dim, dim - 1);

  for (unsigned int i = 0; i < dim - 1; ++i) {
    const unsigned int I1 = i + 1;
    const unsigned int I2 = i + 2;

    const double a = std::sqrt(static_cast<double>(I1) / I2);
    const double b = 1.0 / std::sqrt(static_cast<double>(I1) * I2);

    const double l0 = std::exp(b);
    const double r0 = std::exp(-a);
    const double denom = I1 * l0 + r0 + (dim - I2);

    const double l = l0 / denom;
    const double r = r0 / denom;
    const double one = 1.0 / denom;

    double* col = B.colptr(i);

    for (unsigned int j = 0; j < I1; ++j) {
      col[j] = l;
    }
    col[I1] = r;
    for (unsigned int j = I2; j < dim; ++j) {
      col[j] = one;
    }
  }

  return B;
}



// [[Rcpp::export]]
arma::mat clr_coordinates(const arma::mat& X) {
  arma::mat logX = arma::log(X);
  arma::colvec m = arma::mean(logX, 1);
  logX.each_col() -= m;
  return logX;
}

// [[Rcpp::export]]
arma::mat inv_clr_coordinates(const arma::mat& clrX) {
  arma::mat X = arma::exp(clrX);
  arma::colvec sumX = arma::sum(X, 1);
  X.each_col() /= sumX;
  return X;
}

// [[Rcpp::export]]
arma::mat alr_coordinates(const arma::mat& X, unsigned int denominator) {
  if (denominator < 1 || denominator > X.n_cols) {
    Rcpp::stop("denominator out of bounds.");
  }

  const unsigned int jdenom = denominator - 1;

  arma::mat logX = arma::log(X);
  arma::colvec vdenom = logX.col(jdenom);

  arma::uvec idx = arma::regspace<arma::uvec>(0, X.n_cols - 1);
  idx.shed_row(jdenom);

  arma::mat res = logX.cols(idx);
  res.each_col() -= vdenom;

  return res;
}

// [[Rcpp::export]]
arma::mat ilr_coordinates(const arma::mat& X) {
  arma::mat B = ilr_basis_default(X.n_cols);
  return arma::log(X) * B;
}

// [[Rcpp::export]]
arma::mat inv_ilr_coordinates(const arma::mat& ilrX) {
  arma::mat B = ilr_basis_default(ilrX.n_cols + 1);
  arma::mat clrX = ilrX * B.t();
  return inv_clr_coordinates(clrX);
}

// [[Rcpp::export]]
arma::mat matrix_coordinates(const arma::mat& X, const arma::mat& B) {
  return arma::log(X) * B;
}

// [[Rcpp::export]]
arma::mat sparse_coordinates(const arma::mat& X, const arma::sp_mat& B) {
  return arma::log(X) * B;
}

// [[Rcpp::export]]
arma::mat ilr_to_alr(unsigned int dim) {
  const unsigned int k = dim - 1;
  arma::mat B = ilr_basis_default(dim);
  arma::mat ilr_to_alr_mat = B(arma::span(0, k - 1), arma::span(0, k - 1));

  for (unsigned int i = 0; i < k; ++i) {
    ilr_to_alr_mat.col(i) -= B(k, i);
  }

  return ilr_to_alr_mat;
}
