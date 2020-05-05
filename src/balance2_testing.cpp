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
arma::vec testing_04(arma::mat X) {
  int D = X.n_cols;
  arma::mat lX = log(X);
  arma::vec eigval, V;
  arma::mat eigvec;
  arma::mat S;

  arma::mat B1 = ilr_basis_default(D);
  S = cov(lX * B1);
  arma::eig_sym( eigval, eigvec, S);
  V = B1 * eigvec.tail_cols(1);

  Balance2 balance = Balance2(D);
  balance.approximateLogContrast(V);
  balance.print();

  arma::mat Q, R;
  qr(Q, R, B1.t() * balance.getBalance());

  arma::mat B2 = Q.tail_cols(D-2);
  S = cov(lX * B1 * B2);
  arma::eig_sym( eigval, eigvec, S);
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

// [[Rcpp::export]]
arma::mat testing_05(arma::mat X){
  int D = X.n_cols;
  arma::mat lX = log(X);
  arma::vec eigval, V;
  arma::mat eigvec;
  arma::mat S;
  arma::mat Q, R;

  std::vector<Balance2> balance_candidates;
  arma::mat pb_mat = arma::zeros(D,D-1);

  Balance2 balance = Balance2(D);
  balance_candidates.push_back(balance);

  // ITERATION 1
  arma::mat B1 = ilr_basis_default(D);

  S = cov(lX * B1);
  arma::eig_sym( eigval, eigvec, S);
  V = B1 * eigvec.tail_cols(1);

  int best_balance = 0;
  double best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V, &lX);
    balance_candidates[i].print();
    if(dp > best_dp){
      best_balance = i;
    }
  }
  pb_mat.col(0) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);

  Rcout << "Current list:" << std::endl;
  for(Balance2 balance: balance_candidates){
    balance.print();
  }

  // ITERATION 2
  qr(Q, R, B1.t() * pb_mat.head_cols(1));
  arma::mat B2 = Q.tail_cols(D-2);

  S = cov(lX * B1 * B2);
  arma::eig_sym( eigval, eigvec, S);
  V = B1 * B2 * eigvec.tail_cols(1);

  best_balance = 0;
  best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V, &lX);
    balance_candidates[i].print();
    if(dp > best_dp){
      best_dp = dp;
      best_balance = i;
    }
  }
  pb_mat.col(1) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);

  Rcout << "Current list:" << std::endl;
  for(Balance2 balance: balance_candidates){
    balance.print();
  }

  // ITERATION 3
  qr(Q, R, (B1*B2).t() * pb_mat.head_cols(2));
  arma::mat B3 = Q.tail_cols(D-3);

  S = cov(lX * B1 * B2 * B3);
  arma::eig_sym( eigval, eigvec, S);
  V = B1 * B2 * B3 * eigvec.tail_cols(1);

  best_balance = 0;
  best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V, &lX);
    balance_candidates[i].print();
    if(dp > best_dp){
      best_dp = dp;
      best_balance = i;
    }
  }
  pb_mat.col(2) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);

  Rcout << "Current list:" << std::endl;
  for(Balance2 balance: balance_candidates){
    balance.print();
  }

  // ITERATION 4
  qr(Q, R, (B1*B2*B3).t() * pb_mat.head_cols(3));
  arma::mat B4 = Q.tail_cols(D-4);

  S = cov(lX * B1 * B2 * B3 * B4);
  arma::eig_sym( eigval, eigvec, S);
  V = B1 * B2 * B3 * B4 * eigvec.tail_cols(1);

  best_balance = 0;
  best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V, &lX);
    balance_candidates[i].print();
    if(dp > best_dp){
      best_dp = dp;
      best_balance = i;
    }
  }

  pb_mat.col(3) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);

  Rcout << "Current list:" << std::endl;
  for(Balance2 balance: balance_candidates){
    balance.print();
  }

  // ITERATION 5
  qr(Q, R, (B1*B2*B3*B4).t() * pb_mat.head_cols(4));
  arma::mat B5 = Q.tail_cols(D-5);

  S = cov(lX * B1 * B2 * B3 * B4 * B5);
  arma::eig_sym( eigval, eigvec, S);
  V = B1 * B2 * B3 * B4 * B5 * eigvec.tail_cols(1);

  best_balance = 0;
  best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V, &lX);
    balance_candidates[i].print();
    if(dp > best_dp){
      best_dp = dp;
      best_balance = i;
    }
  }

  pb_mat.col(4) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);

  Rcout << "Current list:" << std::endl;
  for(Balance2 balance: balance_candidates){
    balance.print();
  }

  // ITERATION 6
  qr(Q, R, (B1*B2*B3*B4*B5).t() * pb_mat.head_cols(5));
  arma::mat B6 = Q.tail_cols(D-6);

  S = cov(lX * B1 * B2 * B3 * B4 * B5 * B6);
  arma::eig_sym( eigval, eigvec, S);
  V = B1 * B2 * B3 * B4 * B5 * B6 * eigvec.tail_cols(1);

  best_balance = 0;
  best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V, &lX);
    balance_candidates[i].print();
    if(dp > best_dp){
      best_dp = dp;
      best_balance = i;
    }
  }

  pb_mat.col(5) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);

  Rcout << "Current list:" << std::endl;
  for(Balance2 balance: balance_candidates){
    balance.print();
  }

  return(pb_mat);
}

/*** R
set.seed(3)
X = matrix(rlnorm(10*7), ncol = 7)
# testing_01(X, PC1)
# testing_02()
# testing_03()
testing_05(X)
pb_basis(X, method = 'exact')
*/
