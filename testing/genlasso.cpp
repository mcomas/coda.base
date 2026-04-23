// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

inline double soft_thresh(double x, double t) {
  if (x > t)  return x - t;
  if (x < -t) return x + t;
  return 0.0;
}

// [[Rcpp::export]]
Rcpp::List generalized_lasso_admm_cpp(const arma::mat& X,
                                      const arma::vec& y,
                                      const arma::mat& D,
                                      const double lambda,
                                      const bool intercept = false,
                                      const double rho = 1.0,
                                      const double tol = 1e-4,
                                      const int max_iter = 1000) {

  if (X.n_rows == 0) stop("X no pot estar buit.");
  if (X.n_rows != y.n_rows) stop("nrow(X) ha de coincidir amb length(y).");
  if (X.n_cols != D.n_cols) stop("ncol(X) ha de coincidir amb ncol(D).");
  if (lambda < 0.0) stop("lambda ha de ser >= 0.");
  if (rho <= 0.0) stop("rho ha de ser > 0.");
  if (tol <= 0.0) stop("tol ha de ser > 0.");
  if (max_iter <= 0) stop("max_iter ha de ser > 0.");

  const int n = X.n_rows;
  const int p = X.n_cols;
  const int m = D.n_rows;

  int pDesign = p + (intercept ? 1 : 0);

  // Xdesign
  mat Xd(n, pDesign);
  if (intercept) {
    Xd.col(0).ones();
    Xd.cols(1, pDesign - 1) = X;
  } else {
    Xd = X;
  }

  // D_full (no penalitzem intercept)
  mat Dfull(m, pDesign, fill::zeros);
  if (intercept) {
    Dfull.cols(1, pDesign - 1) = D;
  } else {
    Dfull = D;
  }

  // OBJECTIU: 0.5 ||y - Xd beta||^2 + lambda ||Dfull beta||_1
  mat XtX = Xd.t() * Xd;
  vec XtY = Xd.t() * y;
  mat DtD = Dfull.t() * Dfull;

  mat A = XtX + rho * DtD;

  // Cholesky en forma "lower": A = R * R.t(), R triangular inferior
  mat R;
  bool chol_ok = arma::chol(R, A, "lower");
  if (!chol_ok) {
    stop("No s'ha pogut fer la descomposició de Cholesky de A (potser singular o mal condicionada).");
  }

  vec beta(pDesign, fill::zeros);
  vec z(m, fill::zeros);
  vec u(m, fill::zeros);

  const double lambda_over_rho = lambda / rho;
  const double eps_abs = tol;
  const double eps_rel = tol;

  int iters = 0;

  for (int k = 0; k < max_iter; ++k) {
    iters = k + 1;

    // (1) beta-update:
    // (XtX + rho D'D) beta = XtY + rho D'(z - u)
    vec c   = z - u;
    vec DtC = Dfull.t() * c;
    vec rhs = XtY + rho * DtC;

    // Resol A * beta_new = rhs amb A = R * R.t() (R inferior)
    vec w        = solve(trimatl(R),   rhs);   // R * w = rhs
    vec beta_new = solve(trimatu(R.t()), w);   // R.t() * beta_new = w

    // (2) z-update: v = D beta, z = S_{lambda/rho}(v + u)
    vec v = Dfull * beta_new;
    vec v_plus_u = v + u;

    vec z_new(m);
    for (int i = 0; i < m; ++i) {
      z_new[i] = soft_thresh(v_plus_u[i], lambda_over_rho);
    }

    // (3) u-update
    vec u_new = u + (v - z_new);

    // Residus ADMM (Boyd)
    vec r = v - z_new;                         // primal
    vec s = rho * Dfull.t() * (z_new - z);     // dual

    double r_norm = norm(r, 2);
    double s_norm = norm(s, 2);
    double v_norm = norm(v, 2);
    double z_norm = norm(z_new, 2);
    double u_norm = norm(u_new, 2);

    double eps_pri  = std::sqrt((double)m)       * eps_abs +
      eps_rel * std::max(v_norm, z_norm);
    double eps_dual = std::sqrt((double)pDesign) * eps_abs +
      eps_rel * norm(rho * Dfull.t() * u_new, 2);

    beta = beta_new;
    z    = z_new;
    u    = u_new;

    if (r_norm <= eps_pri && s_norm <= eps_dual) {
      break;
    }
  }

  return Rcpp::List::create(
    _["beta"]      = beta,
    _["z"]         = z,
    _["u"]         = u,
    _["n"]         = n,
    _["p"]         = p,
    _["m"]         = m,
    _["lambda"]    = lambda,
    _["rho"]       = rho,
    _["tol"]       = tol,
    _["iter"]      = iters,
    _["intercept"] = intercept
  );
}
