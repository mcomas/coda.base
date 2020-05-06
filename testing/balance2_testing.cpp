#include <RcppArmadillo.h>
#include "balance2.h"
#include "coda.h"


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
  arma::mat lX = log(X), eigvec, S, Q, R;
  arma::vec eigval, V;
  std::vector<Balance2> balance_candidates;
  arma::mat pb_mat = arma::zeros(D,D-1), v_mat = arma::zeros(D,D-1);
  Balance2 balance = Balance2(D);


  balance_candidates.push_back(balance);

  // ITERATION 1
  arma::mat B = ilr_basis_default(D);

  S = cov(lX * B);
  arma::eig_sym( eigval, eigvec, S);
  V = B * eigvec.tail_cols(1);

  int best_balance = 0;
  double best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V, &lX);
    // balance_candidates[i].print();
    if(dp > best_dp){
      best_balance = i;
    }
  }
  v_mat.col(0) = arma::mat(V);
  pb_mat.col(0) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);
  for(int k=1; k < D-1; k++){
    qr(Q, R, B.t() * pb_mat.col(k-1));
    // Rcout << Q.tail_cols(D-k-1) << std::endl;
    B = B * Q.tail_cols(D-k-1);
    // Rcout << B;
    S = cov(lX * B);
    arma::eig_sym( eigval, eigvec, S);
    V = B * eigvec.tail_cols(1);

    best_balance = 0;
    best_dp = 0;
    for(int i=0; i < balance_candidates.size(); i++){
      double dp = balance_candidates[i].approximateLogContrast(V, &lX);
      if(dp > best_dp){
        best_dp = dp;
        best_balance = i;
      }
    }
    v_mat.col(k) = arma::mat(V);
    pb_mat.col(k) = balance_candidates[best_balance].getBalance();

    balance = balance_candidates[best_balance].top();
    if(balance.size() > 1) balance_candidates.push_back(balance);
    balance = balance_candidates[best_balance].left();
    if(balance.size() > 1) balance_candidates.push_back(balance);
    balance = balance_candidates[best_balance].right();
    if(balance.size() > 1) balance_candidates.push_back(balance);

    // Balance removed from the partitions available
    balance_candidates.erase(balance_candidates.begin()+best_balance);
    // Rcout << "Current list:" << std::endl;
    // for(Balance2 balance: balance_candidates){
    //   balance.print();
    // }
  }
  // List::create(pb_mat, v_mat)
  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat testing_06(arma::mat X){
  int D = X.n_cols;
  arma::mat lX = log(X), eigvec, S, Q, R;
  arma::vec eigval, V;
  std::vector<Balance2> balance_candidates;
  arma::mat pb_mat = arma::zeros(D,D-1), v_mat = arma::zeros(D,D-1);
  Balance2 balance = Balance2(D);


  balance_candidates.push_back(balance);

  // ITERATION 1
  arma::mat B = ilr_basis_default(D);

  S = cov(lX * B);
  arma::eig_sym( eigval, eigvec, S);
  V = B * eigvec.tail_cols(1);

  int best_balance = 0;
  double best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V, &lX);
    // balance_candidates[i].print();
    if(dp > best_dp){
      best_balance = i;
    }
  }
  v_mat.col(0) = arma::mat(V);
  pb_mat.col(0) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);
  for(int k=1; k < D-1; k++){
    qr(Q, R, B.t() * pb_mat.head_cols(k));
    // Rcout << Q.tail_cols(D-k-1) << std::endl;
    arma::mat B_orth = B * Q.tail_cols(D-k-1);
    // Rcout << B;
    S = cov(lX * B_orth);
    arma::eig_sym( eigval, eigvec, S);
    V = B_orth * eigvec.tail_cols(1);

    best_balance = 0;
    best_dp = 0;
    for(int i=0; i < balance_candidates.size(); i++){
      double dp = balance_candidates[i].approximateLogContrast(V, &lX);
      if(dp > best_dp){
        best_dp = dp;
        best_balance = i;
      }
    }
    v_mat.col(k) = arma::mat(V);
    pb_mat.col(k) = balance_candidates[best_balance].getBalance();

    balance = balance_candidates[best_balance].top();
    if(balance.size() > 1) balance_candidates.push_back(balance);
    balance = balance_candidates[best_balance].left();
    if(balance.size() > 1) balance_candidates.push_back(balance);
    balance = balance_candidates[best_balance].right();
    if(balance.size() > 1) balance_candidates.push_back(balance);

    // Balance removed from the partitions available
    balance_candidates.erase(balance_candidates.begin()+best_balance);
    // Rcout << "Current list:" << std::endl;
    // for(Balance2 balance: balance_candidates){
    //   balance.print();
    // }
  }
  // List::create(pb_mat, v_mat)
  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat testing_07(arma::mat X){
  int D = X.n_cols;
  arma::mat lX = log(X), eigvec, S, Q, R;
  arma::vec eigval, V;
  std::vector<Balance2> balance_candidates;
  arma::mat pb_mat = arma::zeros(D,D-1), v_mat = arma::zeros(D,D-1);
  Balance2 balance = Balance2(D);


  balance_candidates.push_back(balance);

  // ITERATION 1
  arma::mat B = ilr_basis_default(D);

  S = cov(lX * B);
  arma::eig_sym( eigval, eigvec, S);
  V = B * eigvec.tail_cols(1);

  int best_balance = 0;
  double best_dp = 0;
  for(int i=0; i < balance_candidates.size(); i++){
    double dp = balance_candidates[i].approximateLogContrast(V);
    // balance_candidates[i].print();
    if(dp > best_dp){
      best_balance = i;
    }
  }
  v_mat.col(0) = arma::mat(V);
  pb_mat.col(0) = balance_candidates[best_balance].getBalance();

  balance = balance_candidates[best_balance].top();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].left();
  if(balance.size() > 1) balance_candidates.push_back(balance);
  balance = balance_candidates[best_balance].right();
  if(balance.size() > 1) balance_candidates.push_back(balance);

  // Balance removed from the partitions available
  balance_candidates.erase(balance_candidates.begin()+best_balance);
  for(int k=1; k < D-1; k++){
    qr(Q, R, B.t() * pb_mat.head_cols(k));
    // Rcout << Q.tail_cols(D-k-1) << std::endl;
    arma::mat B_orth = B * Q.tail_cols(D-k-1);
    // Rcout << B;
    S = cov(lX * B_orth);
    arma::eig_sym( eigval, eigvec, S);
    V = B_orth * eigvec.tail_cols(1);

    best_balance = 0;
    best_dp = 0;
    for(int i=0; i < balance_candidates.size(); i++){
      double dp = balance_candidates[i].approximateLogContrast(V);
      if(dp > best_dp){
        best_dp = dp;
        best_balance = i;
      }
    }
    v_mat.col(k) = arma::mat(V);
    pb_mat.col(k) = balance_candidates[best_balance].getBalance();

    balance = balance_candidates[best_balance].top();
    if(balance.size() > 1) balance_candidates.push_back(balance);
    balance = balance_candidates[best_balance].left();
    if(balance.size() > 1) balance_candidates.push_back(balance);
    balance = balance_candidates[best_balance].right();
    if(balance.size() > 1) balance_candidates.push_back(balance);

    // Balance removed from the partitions available
    balance_candidates.erase(balance_candidates.begin()+best_balance);
    // Rcout << "Current list:" << std::endl;
    // for(Balance2 balance: balance_candidates){
    //   balance.print();
    // }
  }
  // List::create(pb_mat, v_mat)
  return(pb_mat);
}

/*** R
set.seed(2)
X = matrix(rlnorm(10*7), ncol = 7)
# testing_01(X, PC1)
# testing_02()
# testing_03()

A6 = testing_06(X)
# pb_basis(X, method = 'exact')
*/
