#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "balance.h"

std::map<int,arma::uvec> default_nodes(int size){
  std::map<int,arma::uvec> node;
  for(int i=0;i<size;i++){
    node[i] = arma::uvec(1);
    node[i][0] = i;
  }
  return(node);
}

class MaxVariance: public EvaluateBalance {
  arma::mat logX;
  arma::mat M;

public:
  MaxVariance (Balance *bal_, arma::mat X): EvaluateBalance(bal_){
    logX = log(X);
    M = cov(logX);
  }
  double eval(){
    double nL = bal->get_nL();
    double nR = bal->get_nR();
    arma::uvec uL = bal->getL();
    arma::uvec uR = bal->getR();
    double sL = (nR/nL) * arma::accu(M(uL,uL));
    double sR = (nL/nR) * arma::accu(M(uR,uR));
    double sC = - 2*arma::accu(M(uR,uL));

    double variance = (sL + sR + sC) / (nL+nR);
    return variance;
  }

};

//' @export
// [[Rcpp::export]]
arma::vec sbp2_test_1(arma::mat X){

  Balance balance = Balance(default_nodes(X.n_cols));
  MaxVariance h = MaxVariance(&balance, X);
  double score = h.setOptimal();

  arma::vec bal = arma::zeros(X.n_cols);
  bal(balance.getL()).fill(1);
  bal(balance.getR()).fill(-1);
  return bal;

}
