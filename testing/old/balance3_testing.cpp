#include <RcppArmadillo.h>
#include "../src/balance3.h"
#include "../src/evaluate_balance3.h"
#include "../src/coda.h"

// [[Rcpp::export]]
void testing_01(arma::mat X) {
  std::map<int,arma::uvec> nodes;
  nodes[0] << 0 << 4;
  nodes[1] << 1 << 3;
  nodes[2] << 2;
  nodes[3] << 5;
  nodes[4] << 6;
  Balance3<MaximumVariance3a> balance = Balance3<MaximumVariance3a>(X.n_cols, nodes);
  MaximumVariance3a ebalance = MaximumVariance3a(nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();


  balance.print();

}

// [[Rcpp::export]]
void testing_02(arma::mat X) {

  Balance3<MaximumVariance3a> balance = Balance3<MaximumVariance3a>(X.n_cols);
  MaximumVariance3a ebalance = MaximumVariance3a(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
  balance.print();

}

// [[Rcpp::export]]
void testing_03(arma::mat X) {
  Balance3<MaximumVariance3a> balance1 = Balance3<MaximumVariance3a>(X.n_cols);
  Balance3<MaximumVariance3b> balance2 = Balance3<MaximumVariance3b>(X.n_cols);
  MaximumVariance3a ebalance1 = MaximumVariance3a(balance1.nodes, X);
  MaximumVariance3b ebalance2 = MaximumVariance3b(balance2.nodes, X);
  balance1.setEvaluator(ebalance1);
  balance2.setEvaluator(ebalance2);
  clock_t t0 = clock();
  balance1.setWithExhaustiveSearch();
  clock_t t1 = clock();
  balance2.setWithExhaustiveSearch();
  clock_t t2 = clock();
  Rcpp::Rcout << "Timing MaximumVariance3a:" << (double) (t1-t0) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcpp::Rcout << "Timing MaximumVariance3b:" << (double) (t2-t1) / CLOCKS_PER_SEC * 1000.0 << std::endl;

}

/*** R
set.seed(2)
X = matrix(rlnorm(10*7), ncol = 7)

A6 = testing_01(X)
*/
