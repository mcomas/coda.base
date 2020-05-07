#include <RcppArmadillo.h>
#include "../src/balance2.h"
#include "../src/coda.h"


// [[Rcpp::export]]
void eb_test_01a(arma::mat X, arma::vec V) {
  Balance2 balance = Balance2(X.n_cols);
  balance.addL(0);
  balance.addR(1);
  MaximumDotProduct ebalance1 = MaximumDotProduct(V);
  Rcpp::Rcout << "Current value:" << ebalance1.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance1.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance1.evalIfAddR(&balance, 2) << std::endl;
  Rcpp::Rcout << "Optimised minimum angle:" << balance.iterateLogContrast<MaximumDotProduct>(V, &ebalance1) << std::endl;
  balance.print();
}

// [[Rcpp::export]]
void eb_test_01b(arma::mat X, arma::vec V) {
  Balance2 balance = Balance2(X.n_cols);
  balance.addL(0);
  balance.addR(1);
  MaximumVariance1 ebalance2 = MaximumVariance1(X);
  Rcpp::Rcout << "Current value:" << ebalance2.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance2.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance2.evalIfAddR(&balance, 2) << std::endl;
  Rcpp::Rcout << "Optimised maximum variance (method 1):" << balance.iterateLogContrast<MaximumVariance1>(V, &ebalance2) << std::endl;
  balance.print();
}

// [[Rcpp::export]]
void eb_test_01c(arma::mat X, arma::vec V) {
  Balance2 balance = Balance2(X.n_cols);
  balance.addL(0);
  balance.addR(1);
  MaximumVariance3 ebalance2 = MaximumVariance3(X);
  Rcpp::Rcout << "Current value:" << ebalance2.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance2.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance2.evalIfAddR(&balance, 2) << std::endl;
  Rcpp::Rcout << "Optimised maximum variance (method 3):" << balance.iterateLogContrast<MaximumVariance3>(V, &ebalance2) << std::endl;
  balance.print();
}

// [[Rcpp::export]]
void eb_test_02(arma::mat X, arma::vec V) {
  Balance2 balance = Balance2(X.n_cols);
  balance.addL(0);
  balance.addR(1);
  MaximumDotProduct ebalance1 = MaximumDotProduct(V);
  MaximumVariance1 ebalance2 = MaximumVariance1(X);
  MaximumVariance2 ebalance3 = MaximumVariance2(X);
  MaximumVariance3 ebalance4 = MaximumVariance3(X);
  clock_t t0 = clock();
  Rcpp::Rcout << "Optimised minimum angle:" << balance.iterateLogContrast<MaximumDotProduct>(V, &ebalance1) << std::endl;
  clock_t t1 = clock();
  Rcpp::Rcout << "Optimised maximum variance1:" << balance.iterateLogContrast<MaximumVariance1>(V, &ebalance2) << std::endl;
  clock_t t2 = clock();
  Rcpp::Rcout << "Optimised maximum variance2:" << balance.iterateLogContrast<MaximumVariance2>(V, &ebalance3) << std::endl;
  clock_t t3 = clock();
  Rcpp::Rcout << "Optimised maximum variance3:" << balance.iterateLogContrast<MaximumVariance3>(V, &ebalance4) << std::endl;
  clock_t t4 = clock();
  Rcout << "Timing angle:" << (double) (t1-t0) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance1:" << (double) (t2-t1) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance2:" << (double) (t3-t2) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance3:" << (double) (t4-t3) / CLOCKS_PER_SEC * 1000.0 << std::endl;
}

// [[Rcpp::export]]
void eb_test_02b(arma::mat X, arma::vec V) {
  std::map<int,arma::uvec> nodes;
  nodes[0] << 0 << 4;
  nodes[1] << 1 << 3;
  nodes[2] << 2;
  nodes[3] << 5;
  nodes[4] << 6;
  Balance2 balance = Balance2(X.n_cols, nodes);
  balance.addL(0);
  balance.addR(1);
  MaximumDotProduct ebalance1 = MaximumDotProduct(V);
  MaximumVariance1 ebalance2 = MaximumVariance1(X);
  MaximumVariance2 ebalance3 = MaximumVariance2(X);
  MaximumVariance3 ebalance4 = MaximumVariance3(X);
  clock_t t0 = clock();
  Rcpp::Rcout << "Optimised minimum angle:" << balance.iterateLogContrast<MaximumDotProduct>(V, &ebalance1) << std::endl;
  clock_t t1 = clock();
  Rcpp::Rcout << "Optimised maximum variance1:" << balance.iterateLogContrast<MaximumVariance1>(V, &ebalance2) << std::endl;
  clock_t t2 = clock();
  Rcpp::Rcout << "Optimised maximum variance2:" << balance.iterateLogContrast<MaximumVariance2>(V, &ebalance3) << std::endl;
  clock_t t3 = clock();
  Rcpp::Rcout << "Optimised maximum variance3:" << balance.iterateLogContrast<MaximumVariance3>(V, &ebalance4) << std::endl;
  clock_t t4 = clock();
  Rcout << "Timing angle:" << (double) (t1-t0) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance1:" << (double) (t2-t1) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance2:" << (double) (t3-t2) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance3:" << (double) (t4-t3) / CLOCKS_PER_SEC * 1000.0 << std::endl;
}

// [[Rcpp::export]]
void eb_test_03b(arma::mat X, arma::vec V) {
  std::map<int,arma::uvec> nodes;
  nodes[0] << 0 << 4;
  nodes[1] << 1 << 3;
  nodes[2] << 2;
  nodes[3] << 5;
  nodes[4] << 6;
  Balance2 balance = Balance2(X.n_cols, nodes);
  balance.addL(0);
  balance.addR(1);
  MaximumVariance1 ebalance2 = MaximumVariance1(X);
  balance.print();
  Rcpp::Rcout << "Current value:" << ebalance2.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance2.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance2.evalIfAddR(&balance, 2) << std::endl;
  Rcpp::Rcout << "Optimised maximum variance (method 1):" << balance.iterateLogContrast<MaximumVariance1>(V, &ebalance2) << std::endl;
  balance.print();
}

// [[Rcpp::export]]
void eb_test_03c(arma::mat X, arma::vec V) {
  std::map<int,arma::uvec> nodes;
  nodes[0] << 0 << 4;
  nodes[1] << 1 << 3;
  nodes[2] << 2;
  nodes[3] << 5;
  nodes[4] << 6;
  Balance2 balance = Balance2(X.n_cols, nodes);
  balance.addL(0);
  balance.addR(1);
  MaximumVariance3 ebalance2 = MaximumVariance3(X);
  balance.print();
  Rcpp::Rcout << "Current value:" << ebalance2.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance2.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance2.evalIfAddR(&balance, 2) << std::endl;
  Rcpp::Rcout << "Optimised maximum variance (method 3):" << balance.iterateLogContrast<MaximumVariance3>(V, &ebalance2) << std::endl;
  balance.print();
}

// [[Rcpp::export]]
void eb_test_04c(arma::mat X, arma::vec V) {
  std::map<int,arma::uvec> nodes;
  nodes[0] << 0 << 4;
  nodes[1] << 1 << 3;
  nodes[2] << 2;
  nodes[3] << 5;
  nodes[4] << 6;
  Balance2 balance = Balance2(X.n_cols, nodes);
  balance.addL(0);
  balance.addR(1);
  MaximumVariance3 ebalance2 = MaximumVariance3(X);
  balance.print();
  Rcpp::Rcout << "Current value:" << ebalance2.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance2.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance2.evalIfAddR(&balance, 2) << std::endl;
  arma::mat lX = log(X);
  Rcpp::Rcout << "Optimised maximum variance (balance in):" << balance.approximateLogContrast(V, &lX) << std::endl;
  balance.print();
}

/*** R
set.seed(2)
X = matrix(rlnorm(10*7), ncol = 7)
eb_test_01a(X, eigen(cov(X))$vector[,1,drop=FALSE])
eb_test_01b(X, eigen(cov(X))$vector[,1,drop=FALSE])
eb_test_01c(X, eigen(cov(X))$vector[,1,drop=FALSE])
*/
