#include <RcppArmadillo.h>
#include "balance2.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
void testing_01(arma::mat X, arma::vec V) {
  std::map<int,arma::uvec> nodes;
  nodes[0] << 0 << 4;
  nodes[1] << 1 << 3;
  nodes[2] << 2;
  nodes[3] << 5;
  nodes[4] << 6;
  Balance2 balance = Balance2(7, nodes);
  balance.approximateLogContrast(V);
  // clock_t t1 = clock();
  // Rcout << "Balance:" << std::endl << balance.getBalance();
  // clock_t t2 = clock();
  // Rcout << "Balance Alt.:" << std::endl << balance.getBalance2();
  // clock_t t3 = clock();
  // Rcout << t2-t1 << std::endl << t3-t2 << std::endl;
  balance.print();
}

// [[Rcpp::export]]
void testing_02() {
  std::map<int,arma::uvec> nodes;
  nodes[0] << 0;
  nodes[1] << 1;
  nodes[2] << 2;
  nodes[3] << 3;
  nodes[4] << 4;
  Balance2 balance = Balance2(5, nodes);
  arma::uvec uL, uR;
  uL << 0 << 1;
  uR << 2 << 3;
  balance.setL(uL);
  balance.setR(uR);
  balance.print();
  Rcout << balance.getBalance();
}

// [[Rcpp::export]]
void testing_03() {
  std::map<int,arma::uvec> nodes;
  nodes[0] << 0 << 4;
  nodes[1] << 1 << 3;
  nodes[2] << 2;
  nodes[3] << 5;
  nodes[4] << 6;
  Balance2 balance = Balance2(7, nodes);
  arma::uvec uL, uR;
  uL << 0;
  uR << 2;
  balance.addL(0);
  balance.setR(uR);
  balance.print();
  Rcout << balance.getBalanceIfAddL(3) << std::endl;
  Rcout << balance.getBalanceIfAddR(3) << std::endl;
  arma::uvec iL;
  arma::uvec iR(1);
  iL << 3 << 1; iR(0) = 4;
  Rcout << balance.getBalanceIfAdd(iL, iR);
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
set.seed(1)
X = matrix(rlnorm(10*7), ncol = 7)
PC1 = coda.base::pc_basis(X)[,1,drop=FALSE]
testing_01(X, PC1)
testing_02()
testing_03()
*/
