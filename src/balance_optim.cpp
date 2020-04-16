#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"


class PredictiveBalance: public EvaluateBalance {
  arma::mat lY;
  arma::mat x;
public:
  PredictiveBalance (arma::mat Y, arma::mat X): EvaluateBalance(X.n_cols){
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
  PredictiveBalance2 (arma::mat Y, arma::mat X): EvaluateBalance(X.n_cols){
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


// [[Rcpp::export]]
arma::vec find_predictive_balance(arma::mat Y, arma::vec x, int method = 1){
  EvaluateBalance h;
  if(method == 1){
    h = PredictiveBalance(Y, x);
    double score = h.setOptimal();
  }
  if(method == 2){
    h = PredictiveBalance2(Y, x);
    double score = h.setOptimal();
  }


  arma::vec bal = arma::zeros(Y.n_cols);
  bal(h.bal.getL()).fill(1);
  bal(h.bal.getR()).fill(-1);
  return bal;

}
