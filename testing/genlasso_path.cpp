// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

inline double soft_thresh(double x, double t) {
  if (x > t)  return x - t;
  if (x < -t) return x + t;
  return 0.0;
}

//' Lasso/generalized lasso path via ADMM
 //'
 //' Resol per a cada lambda_k del vector `lambdas` el problema
 //'   min_beta  0.5 * ||y - X beta||^2 + lambda_k * ||D beta||_1
 //' (amb opció d'intercept no penalitzat) utilitzant ADMM,
 //' reaprofitant la factorització de (X'X + rho D'D) i
 //' fent warm start de (beta, z, u) entre lambdes.
 //'
 //' @param X matriu n x p
 //' @param y vector n
 //' @param D matriu m x p
 //' @param lambdas vector de penalitzacions (>0)
 //' @param intercept si TRUE, s'afegeix un intercept no penalitzat
 //' @param rho paràmetre ADMM (>0)
 //' @param tol tolerància base (eps_abs i eps_rel)
 //' @param max_iter màxim nombre d'iteracions ADMM per cada lambda
 //' @return Llista amb
 //'   * coef: matriu pDesign x K (una columna per lambda)
 //'           on pDesign = p + (intercept ? 1 : 0)
 //'           (si intercept=TRUE, coef[0,] és l'intercept)
 //'   * iter: enter K amb les iteracions gastades per cada lambda
 //'   * n, p, m, rho, tol, lambdas
 //'
 //' @export
 // [[Rcpp::export]]
 Rcpp::List generalized_lasso_path_admm_cpp(const arma::mat& X,
                                            const arma::vec& y,
                                            const arma::mat& D,
                                            const arma::vec& lambdas,
                                            const bool intercept = false,
                                            const double rho = 1.0,
                                            const double tol = 1e-4,
                                            const int max_iter = 1000) {
   
   if (X.n_rows == 0) stop("X no pot estar buit.");
   if (X.n_rows != y.n_rows) stop("nrow(X) ha de coincidir amb length(y).");
   if (X.n_cols != D.n_cols) stop("ncol(X) ha de coincidir amb ncol(D).");
   if (rho <= 0.0) stop("rho ha de ser > 0.");
   if (tol <= 0.0) stop("tol ha de ser > 0.");
   if (max_iter <= 0) stop("max_iter ha de ser > 0.");
   if (lambdas.n_elem == 0) stop("cal almenys un valor a `lambdas`.");
   
   const int n = X.n_rows;
   const int p = X.n_cols;
   const int m = D.n_rows;
   
   const int K = lambdas.n_elem;
   const int pDesign = p + (intercept ? 1 : 0);
   
   // --- Xdesign (amb o sense intercept) ---
   arma::mat Xd(n, pDesign);
   if (intercept) {
     Xd.col(0).ones();
     Xd.cols(1, pDesign - 1) = X;
   } else {
     Xd = X;
   }
   
   // --- D_full: no penalitzem intercept ---
   arma::mat Dfull(m, pDesign, arma::fill::zeros);
   if (intercept) {
     Dfull.cols(1, pDesign - 1) = D;
   } else {
     Dfull = D;
   }
   
   // --- precomputacions XtX, XtY, D'D, A ---
   // Objectiu: 0.5 ||y - Xd beta||^2 + lambda ||Dfull beta||_1
   arma::mat XtX = Xd.t() * Xd;     // pDesign x pDesign
   arma::vec XtY = Xd.t() * y;      // pDesign
   arma::mat DtD = Dfull.t() * Dfull;
   
   arma::mat A = XtX + rho * DtD;
   
   // Cholesky: A = R * R.t(), R triangular inferior
   arma::mat R;
   bool chol_ok = arma::chol(R, A, "lower");
   if (!chol_ok) {
     stop("No s'ha pogut fer la descomposició de Cholesky de A (potser singular o mal condicionada).");
   }
   
   // --- inicialització per al path: warm start des de lambda[0] ---
   arma::vec beta(pDesign, arma::fill::zeros); // beta (pDesign)
   arma::vec z(m,        arma::fill::zeros);   // D beta
   arma::vec u(m,        arma::fill::zeros);   // dual escalat
   
   arma::mat coef(pDesign, K, arma::fill::zeros);
   arma::ivec iters(K, arma::fill::zeros);
   
   const double eps_abs = tol;
   const double eps_rel = tol;
   
   // recomanable passar lambdes en ordre decreixent per warm start
   for (int k = 0; k < K; ++k) {
     
     double lambda_k = lambdas[k];
     if (lambda_k < 0.0) {
       stop("tots els lambdes han de ser >= 0.");
     }
     
     const double lambda_over_rho = lambda_k / rho;
     int it = 0;
     
     for (; it < max_iter; ++it) {
       
       // (1) beta-update:
       // (XtX + rho D'D) beta^{k+1} = XtY + rho D'(z^k - u^k)
       arma::vec c   = z - u;             // m
       arma::vec DtC = Dfull.t() * c;     // pDesign
       arma::vec rhs = XtY + rho * DtC;   // pDesign
       
       // Resol A * beta_new = rhs amb A = R * R.t(), R inferior
       arma::vec w        = arma::solve(arma::trimatl(R),   rhs);   // R * w = rhs
       arma::vec beta_new = arma::solve(arma::trimatu(R.t()), w);   // R.t() * beta_new = w
       
       // (2) z-update:
       // v = D beta^{k+1}
       // z^{k+1} = S_{lambda/rho}(v + u^k)
       arma::vec v = Dfull * beta_new;
       arma::vec v_plus_u = v + u;
       
       arma::vec z_new(m);
       for (int i = 0; i < m; ++i) {
         z_new[i] = soft_thresh(v_plus_u[i], lambda_over_rho);
       }
       
       // (3) u-update:
       arma::vec u_new = u + (v - z_new);
       
       // Residus ADMM (Boyd):
       arma::vec r = v - z_new;                   // primal
       arma::vec s = rho * Dfull.t() * (z_new - z); // dual
       
       double r_norm = arma::norm(r, 2);
       double s_norm = arma::norm(s, 2);
       double v_norm = arma::norm(v, 2);
       double z_norm = arma::norm(z_new, 2);
       
       double eps_pri  = std::sqrt((double)m)       * eps_abs +
         eps_rel * std::max(v_norm, z_norm);
       double eps_dual = std::sqrt((double)pDesign) * eps_abs +
         eps_rel * arma::norm(rho * Dfull.t() * u_new, 2);
       
       beta = beta_new;
       z    = z_new;
       u    = u_new;
       
       if (r_norm <= eps_pri && s_norm <= eps_dual) {
         break;
       }
     }
     
     coef.col(k) = beta;
     iters[k]    = it + 1; // iteracions realment fetes
   }
   
   return List::create(
     _["coef"]    = coef,
     _["iter"]    = iters,
     _["n"]       = n,
     _["p"]       = p,
     _["m"]       = m,
     _["rho"]     = rho,
     _["tol"]     = tol,
     _["lambdas"] = lambdas,
     _["intercept"] = intercept
   );
 }
