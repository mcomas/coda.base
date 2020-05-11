#include <RcppArmadillo.h>
#include "../src/balance.h"
#include "../src/balance_evaluate.h"
#include "../src/coda.h"

void optimise(Balance<MaximumVariance>& balance, arma::mat &X){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
}

// [[Rcpp::export]]
void testing_01(arma::mat X) {

  Balance<MaximumVariance> balance = Balance<MaximumVariance>(X.n_cols);
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
  balance.print();

  Balance<MaximumVariance> top = balance.top();
  Balance<MaximumVariance> left = balance.left();
  Balance<MaximumVariance> right = balance.right();

  if(top.nodes.size() > 1){
    optimise(top, X);
    top.print();
  }
  if(left.nodes.size() > 1){
    optimise(left, X);
    left.print();
  }
  if(right.nodes.size() > 1){
    optimise(right, X);
    right.print();
  }

}

// [[Rcpp::export]]
void testing_02(arma::mat X) {

  Balance<MaximumVariance> balance = Balance<MaximumVariance>(X.n_cols);
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);

  arma::mat B1 = ilr_basis_default(X.n_cols);
  arma::mat S = cov(log(X) * B1);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym( eigval, eigvec, S);
  arma::vec LC1 = B1 * eigvec.tail_cols(1);
  balance.setWithLogContrast(LC1);
  balance.print();

  // Balance<MaximumVariance> top = balance.top();
  // Balance<MaximumVariance> left = balance.left();
  // Balance<MaximumVariance> right = balance.right();
  //
  // if(top.nodes.size() > 1){
  //   optimise(top, X);
  //   top.print();
  // }
  // if(left.nodes.size() > 1){
  //   optimise(left, X);
  //   left.print();
  // }
  // if(right.nodes.size() > 1){
  //   optimise(right, X);
  //   right.print();
  // }

}


/*** R
set.seed(2)
X = matrix(rlnorm(10*7), ncol = 7)

testing_02(X)
*/
