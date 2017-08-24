// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec eigen1(arma::mat M){
  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, M);
  return(eigvec.tail_cols(1));
}

// [[Rcpp::export]]
arma::vec eigen2(arma::mat M){
  arma::vec eigval;
  arma::mat eigvec;

  arma::eigs_sym(eigval, eigvec, arma::sp_mat(M), 1);
  return(eigvec);
}

// [[Rcpp::export]]
arma::vec eigen3(arma::mat M){
  arma::vec b_old = arma::zeros<arma::vec>(M.n_cols);
  arma::vec b_new = arma::ones<arma::vec>(M.n_cols);
  int iter = 0;
  while( max(abs(b_old-b_new)) > 0.001){
    iter++;
    b_old = b_new;
    b_new = normalise(M * b_old);
  }
  Rcpp::Rcout << iter;
  return(b_new);
}

// [[Rcpp::export]]
arma::vec eigen4(arma::mat M){
  int k = M.n_rows;

  arma::vec b = arma::ones<arma::vec>(k);
  double mu_prev = 0, mu = 1;

  while( std::abs(mu - mu_prev) > 0.001 ){
    b = normalise(M * b);
    mu_prev = mu;
    mu = ((arma::mat)(b.t() * M * b))[0];
  }

  arma::mat inv_M = arma::mat(k, k);
  arma::mat MuId = arma::zeros(k, k);
  MuId.diag().fill(mu);

  for(int iter= 0; iter < 10; iter++){
    b = normalise( inv(M - MuId) * b );
    mu = ((arma::mat)(b.t() * M * b))(0,0);
    MuId.diag().fill( mu );
  }
  return(b);
}
