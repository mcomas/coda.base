#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"

std::map<int,arma::uvec> default_nodes(int size){
  std::map<int,arma::uvec> node;
  for(int i=0;i<size;i++){
    node[i] = arma::uvec(1);
    node[i][0] = i;
  }
  return(node);
}

class PrincipalBalance: public EvaluateBalance {
  arma::mat M;
public:
  PrincipalBalance (Balance *bal_, arma::mat X): EvaluateBalance(bal_){
    M = cov(log(X));
  }
  double eval(Balance *bal){
    double nL = bal->L_length; //bal->get_nL();
    double nR = bal->R_length; //bal->get_nR();
    // Rcpp::Rcout << "nL: " << nL << std::endl;
    // Rcpp::Rcout << "nR: " << nR << std::endl;
    arma::uvec uL = bal->getL();
    arma::uvec uR = bal->getR();
    double sL = (nR/nL) * arma::accu(M(uL,uL));
    double sR = (nL/nR) * arma::accu(M(uR,uR));
    double sC = - 2*arma::accu(M(uR,uL));

    double variance = (sL + sR + sC) / (nL+nR);

    // Rcpp::Rcout << variance << std::endl;
    return variance;
  }
};



class PredictiveBalance: public EvaluateBalance {
  arma::mat lY;
  arma::mat x;
public:
  PredictiveBalance (Balance *bal_, arma::mat Y, arma::mat X): EvaluateBalance(bal_){
    lY = log(Y);
    x = X;
  }
  double eval(Balance *bal){
    arma::uvec uL = bal->getL();
    arma::uvec uR = bal->getR();

    arma::mat y =  mean(lY.cols(bal->getL()), 1) - mean(lY.cols(bal->getR()), 1);
    //Rcpp::Rcout << y << std::endl << x << std::endl;
    arma::mat b = cov(y,x)/var(x);
    arma::mat a = mean(y) - b * mean(x);
    arma::mat r = a(0) + b(0) * x - y;
    arma::mat r2 = r.t() * r;
    //arma::mat c = arma::cov(y,x);
    //Rcpp::Rcout <<c  << std::endl;
    return b(0) * b(0) / r2(0);

  }
};

class PredictiveBalance2: public EvaluateBalance {
  arma::mat lY;
  arma::mat x;
public:
  PredictiveBalance2 (Balance *bal_, arma::mat Y, arma::mat X): EvaluateBalance(bal_){
    lY = log(Y);
    arma::mat m = mean(X);
    arma::mat v = var(X);
    x = (X - m(0))/sqrt(v(0));
  }
  double eval(Balance *bal){
    arma::uvec uL = bal->getL();
    arma::uvec uR = bal->getR();

    arma::mat y =  mean(lY.cols(bal->getL()), 1) - mean(lY.cols(bal->getR()), 1);
    //Rcpp::Rcout << y << std::endl << x << std::endl;
    arma::mat b = cov(y,x);
    arma::mat a = mean(y) - b;
    arma::mat r = a(0) + b(0) * x - y;
    arma::mat r2 = r.t() * r;
    //arma::mat c = arma::cov(y,x);
    //Rcpp::Rcout <<c  << std::endl;
    return b(0) * b(0) / r2(0);

  }
};

//' @export
// [[Rcpp::export]]
arma::vec find_testing(arma::mat X, arma::vec v){

  Balance balance = Balance(default_nodes(X.n_cols));
  // arma::uvec uL;
  // uL << 1 << 2;
  // arma::uvec uR;
  // uR << 0 << 4;
  // balance.init(uL, uR);
  arma::uvec uL0 = find(v == max(v));
  arma::uvec uR0 = find(v == min(v));
  balance.init(uL0, uR0);
  PrincipalBalance h = PrincipalBalance(&balance, X);
  h.print_state();
  h.setLocalSearch();
  h.print_state();
  //double score = h.setOptimal();

  arma::vec bal = arma::zeros(X.n_cols);
  bal(balance.getL()).fill(1);
  bal(balance.getR()).fill(-1);
  return bal;

}

//' @export
// [[Rcpp::export]]
arma::vec find_principal_balance(arma::mat X){

  Balance balance = Balance(default_nodes(X.n_cols));
  PrincipalBalance h = PrincipalBalance(&balance, X);
  double score = h.setOptimal();

  arma::vec bal = arma::zeros(X.n_cols);
  bal(balance.getL()).fill(1);
  bal(balance.getR()).fill(-1);
  return bal;

}


//' @export
// [[Rcpp::export]]
arma::vec find_predictive_balance(arma::mat Y, arma::vec x, int method = 1){

  Balance balance = Balance(default_nodes(Y.n_cols));
  if(method == 1){
    PredictiveBalance h = PredictiveBalance(&balance, Y, x);
    double score = h.setOptimal();
  }
  if(method == 2){
    PredictiveBalance2 h = PredictiveBalance2(&balance, Y, x);
    double score = h.setOptimal();
  }


  arma::vec bal = arma::zeros(Y.n_cols);
  bal(balance.getL()).fill(1);
  bal(balance.getR()).fill(-1);
  return bal;

}
