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

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
set.seed(1)
X = matrix(rlnorm(10*7), ncol = 7)
PC1 = coda.base::pc_basis(X)[,1,drop=FALSE]
testing_01(X, PC1)
*/
