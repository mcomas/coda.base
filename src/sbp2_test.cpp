#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "sbp2.h"

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

    double sL = (nR/nL) * arma::accu(M(bal->getL(),bal->getL()));
    double sR = (nL/nR) * arma::accu(M(bal->getR(),bal->getR()));
    double sC = - 2*arma::accu(M(bal->getR(),bal->getL()));

    double variance = (sL + sR + sC) / (nL+nR);
    return variance;
    // arma::mat h = sqrt(nL*nR/(nL+nR)) * (mean(logX.cols(bal->getL()), 1) - mean(logX.cols(bal->getR()), 1));
    //
    // arma::mat v = var(h);
    // return v(0,0);
  }
  void print_state(){
    bal->print();
    // Rcpp::Rcout << logX.cols(bal->getL()) << std::endl;
    // Rcpp::Rcout << logX.cols(bal->getR());
    Rcpp::Rcout << eval() << std::endl;
  }
  void setOptimal(){
    bal->init();
    //print_state();

    double best_score = eval();
    arma::uvec L = arma::uvec(bal->getL());
    arma::uvec R = arma::uvec(bal->getR());


    unsigned int iter = 1;
    while(bal->hasNext()){
      if(iter % 10000 == 0){
        R_CheckUserInterrupt();
      }
      iter++;
      bal->nextBalance();
      //print_state();

      double score = eval();
      if(score > best_score){
        best_score = score;
        L = arma::uvec(bal->getL());
        R = arma::uvec(bal->getR());
      }
    }

    bal->init(L, R);
    //print_state();
  }
};

//' @export
// [[Rcpp::export]]
int sbp2_test_1(arma::mat X){

  Balance balance = Balance(default_nodes(X.n_cols));
  MaxVariance h = MaxVariance(&balance, X);
  h.setOptimal();

  // balance.init();
  // balance.print();
  // unsigned int iter = 1;
  // while(balance.hasNext()){
  //   if(iter % 10000 == 0){
  //     R_CheckUserInterrupt();
  //   }
  //   iter++;
  //   balance.nextBalance();
  //   balance.print();
  // }
  // return iter;


}
