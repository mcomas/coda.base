// #ifndef balance_optimal_H
// #define balance_optimal_H

#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"

void optimise(Balance<MaximumVariance>& balance, arma::mat& X){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
}

void optimise_using_pc(Balance<MaximumVariance>& balance, arma::mat& X){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  if(balance.nodes.size() == 2){
    arma::uvec uL(1), uR(1);
    uL[0] = 0; uR[0] = 1;
    ebalance.eval(uL, uR, 1, 1);
    balance.set(uL, uR);
  }else{
    arma::mat Xsub = arma::mat(X.n_rows, balance.nodes.size());
    for(int i = 0; i < Xsub.n_cols; i++){
      Xsub.col(i) = X.col(balance.nodes[i][0]);
      for(int j = 1; j < balance.nodes[i].n_elem; j++){
        Xsub.col(i) += X.col(balance.nodes[i][j]);
      }
    }
    arma::mat U, V;
    arma::vec s;

    arma::svd_econ(U, s, V, clr_coordinates(Xsub));
    // Rcpp::Rcout << V.col(0).t();
    balance.setWithLogContrast(V.col(0));

  }
}

void optimise_using_pc_forcing_branch(Balance<MaximumVariance>& balance, arma::mat& X, unsigned forced){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  if(balance.nodes.size() == 2){
    arma::uvec uL(1), uR(1);
    uL[0] = 0; uR[0] = 1;
    ebalance.eval(uL, uR, 1, 1);
    balance.set(uL, uR);
  }else{
    arma::mat Xsub = arma::mat(X.n_rows, balance.nodes.size());
    for(int i = 0; i < Xsub.n_cols; i++){
      Xsub.col(i) = X.col(balance.nodes[i][0]);
      for(int j = 1; j < balance.nodes[i].n_elem; j++){
        Xsub.col(i) += X.col(balance.nodes[i][j]);
      }
    }
    arma::mat U, V;
    arma::vec s;

    arma::svd_econ(U, s, V, clr_coordinates(Xsub));
    arma::uvec uF = arma::uvec(1);
    uF[0] = forced;
    balance.setWithLogContrastForceBranch(V.col(0), uF);

  }
}

void optimise_recursively(Balance<MaximumVariance>& balance, arma::mat& X, arma::mat& pb_mat, int *pb_size){

  optimise_using_pc(balance, X);
  pb_mat.col(*pb_size) = balance.getBalance();
  (*pb_size)++;

  Balance<MaximumVariance> top = balance.top();
  if(top.nodes.size() > 1){
    optimise_recursively(top, X, pb_mat, pb_size);
  }
  Balance<MaximumVariance> left = balance.left();
  if(left.nodes.size() > 1){
    optimise_recursively(left, X, pb_mat, pb_size);
  }
  Balance<MaximumVariance> right  = balance.right();
  if(right.nodes.size() > 1){
    optimise_recursively(right, X, pb_mat, pb_size);
  }

}

void optimise_recursively_forcing_parents(Balance<MaximumVariance>& balance, arma::mat& X, arma::mat& pb_mat, int *pb_size){

  optimise_using_pc(balance, X);
  pb_mat.col(*pb_size) = balance.getBalance();
  (*pb_size)++;

  Balance<MaximumVariance> top = balance.top();
  while(top.nodes.size() > 1){
    optimise_using_pc_forcing_branch(top, X, top.nodes.size()-1);
    pb_mat.col(*pb_size) = top.getBalance();
    (*pb_size)++;
    Balance<MaximumVariance> left = top.left();
    if(left.nodes.size() > 1){
      optimise_recursively_forcing_parents(left, X, pb_mat, pb_size);
    }
    Balance<MaximumVariance> right  = top.right();
    if(right.nodes.size() > 1){
      optimise_recursively_forcing_parents(right, X, pb_mat, pb_size);
    }
    top = top.top();
  }
  Balance<MaximumVariance> left = balance.left();
  if(left.nodes.size() > 1){
    optimise_recursively_forcing_parents(left, X, pb_mat, pb_size);
  }
  Balance<MaximumVariance> right  = balance.right();
  if(right.nodes.size() > 1){
    optimise_recursively_forcing_parents(right, X, pb_mat, pb_size);
  }

}


// [[Rcpp::export]]
arma::mat find_PB(arma::mat& X){
  int K = X.n_cols;

  arma::mat pb_mat = arma::mat(K,K-1);

  std::vector<Balance<MaximumVariance>> SOLS;

  Balance<MaximumVariance> balance = Balance<MaximumVariance>(X.n_cols);

  optimise(balance, X);

  SOLS.push_back(balance);

  for(int l=0;l<K-1;l++){
    //Rcpp::Rcout << "Starting step " << l + 1 << " of " << K-1 << std::endl;
    //print_list(SOLS);
    double vBestSolution = 0;
    int iBestSolution = -1;
    for(unsigned int i =0; i< SOLS.size();i++){
      double v = SOLS[i].eval();
      if(v > vBestSolution){
        vBestSolution = v;
        iBestSolution = i;
      }
    }

    pb_mat.col(l) = SOLS[iBestSolution].getBalance();

    Balance<MaximumVariance> top = SOLS[iBestSolution].top();
    Balance<MaximumVariance> left = SOLS[iBestSolution].left();
    Balance<MaximumVariance> right  = SOLS[iBestSolution].right();

    if(top.nodes.size() > 1){
      optimise(top, X);
      SOLS.push_back(top);
    }
    if(left.nodes.size() > 1){
      optimise(left, X);
      SOLS.push_back(left);
    }
    if(right.nodes.size() > 1){
      optimise(right, X);
      SOLS.push_back(right);
    }

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat find_PB_using_pc(arma::mat& X){
  int K = X.n_cols;

  arma::mat pb_mat = arma::mat(K,K-1);

  std::vector<Balance<MaximumVariance>> SOLS;

  Balance<MaximumVariance> balance = Balance<MaximumVariance>(X.n_cols);

  optimise_using_pc(balance, X);

  SOLS.push_back(balance);

  for(int l=0;l<K-1;l++){
    //Rcpp::Rcout << "Starting step " << l + 1 << " of " << K-1 << std::endl;
    //print_list(SOLS);
    double vBestSolution = 0;
    int iBestSolution = -1;
    for(unsigned int i =0; i< SOLS.size();i++){
      double v = SOLS[i].eval();
      if(v > vBestSolution){
        vBestSolution = v;
        iBestSolution = i;
      }
    }

    pb_mat.col(l) = SOLS[iBestSolution].getBalance();

    Balance<MaximumVariance> top = SOLS[iBestSolution].top();
    Balance<MaximumVariance> left = SOLS[iBestSolution].left();
    Balance<MaximumVariance> right  = SOLS[iBestSolution].right();

    if(top.nodes.size() > 1){
      optimise_using_pc(top, X);
      SOLS.push_back(top);
    }
    if(left.nodes.size() > 1){
      optimise_using_pc(left, X);
      SOLS.push_back(left);
    }
    if(right.nodes.size() > 1){
      optimise_using_pc(right, X);
      SOLS.push_back(right);
    }

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat find_PB_using_pc_recursively(arma::mat& X){
  int K = X.n_cols;

  arma::mat pb_mat = arma::zeros(K,K-1);

  std::vector<Balance<MaximumVariance>> SOLS;
  int pb_size = 0;
  Balance<MaximumVariance> root = Balance<MaximumVariance>(X.n_cols);
  optimise_recursively(root, X, pb_mat, &pb_size);


  // Rcpp::Rcout << pb_size << std::endl;
  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat find_PB_using_pc_recursively_forcing_parents(arma::mat& X){
  int K = X.n_cols;

  arma::mat pb_mat = arma::zeros(K,K-1);

  int pb_size = 0;
  Balance<MaximumVariance> root = Balance<MaximumVariance>(X.n_cols);
  optimise_recursively_forcing_parents(root, X, pb_mat, &pb_size);

  return(pb_mat);
}

/*** R
SEED = round(1000*runif(1))
set.seed(SEED)
SEED
D = 10
X = matrix(rlnorm(100*D), ncol = D)
find_PB(X)
find_PB_using_pc(X)
find_PB_using_pc_recursively(X)
find_PB_using_pc_recursively_forcing_parents(X)
*/

// #endif
