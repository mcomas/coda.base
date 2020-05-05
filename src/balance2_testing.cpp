#include <RcppArmadillo.h>
#include "balance2.h"
#include "coda.h"

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
  Rcout << balance.getBalanceIfAddL(3).t() << std::endl;
  Rcout << balance.getBalanceIfAddR(3).t() << std::endl;
  arma::uvec iL;
  arma::uvec iR(1);
  iL << 3 << 1; iR(0) = 4;
  Rcout << balance.getBalanceIfAdd(iL, iR).t();
}

// [[Rcpp::export]]
arma::vec testing_04(arma::mat X, arma::vec V) {
  int D = X.n_cols;
  Balance2 balance = Balance2(D);
  balance.approximateLogContrast(V);
  balance.print();

  arma::mat B1 = ilr_basis_default(D);
  arma::mat Q, R;
  qr(Q,R,B1.t() * balance.getBalance());
  arma::mat B2 = Q.tail_cols(D-2);
  arma::mat S2 = cov(matrix_coordinates(X, B1 * B2));

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym( eigval, eigvec, S2);
  V = B1 * B2 * eigvec.tail_cols(1);

  Rcout << "Top:" << std::endl;
  Balance2 bal_ = balance.top();
  bal_.print();
  bal_.approximateLogContrast(V);
  bal_.print();

  Rcout << "Left:" << std::endl;
  bal_ = balance.left();
  bal_.print();
  bal_.approximateLogContrast(V);
  bal_.print();

  Rcout << "Right:" << std::endl;
  bal_ = balance.right();
  bal_.print();
  bal_.approximateLogContrast(V);
  bal_.print();

  return(balance.getBalance());
}

/*** R
set.seed(2)
X = matrix(rlnorm(10*7), ncol = 7)
PC1 = coda.base::pc_basis(X)[,1,drop=FALSE]
# testing_01(X, PC1)
# testing_02()
# testing_03()
V = testing_04(X, PC1)
V
*/
