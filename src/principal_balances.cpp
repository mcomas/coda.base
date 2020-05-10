// #ifndef balance_optimal_H
// #define balance_optimal_H

#include <RcppArmadillo.h>
#include "balance.h"

void optimise(Balance<MaximumVariance>& balance, arma::mat& X){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
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


// #endif
