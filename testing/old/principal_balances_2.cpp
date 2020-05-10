#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"

class PrincipalBalance2: public EvaluateBalance {
  arma::mat M;
public:
  PrincipalBalance2(arma::mat X): EvaluateBalance(X.n_cols){
    M = cov(clr_coordinates(X));
  }
  PrincipalBalance2 (std::map<int,arma::uvec> nodes, arma::mat X): EvaluateBalance(nodes, X.n_cols){
    M = cov(clr_coordinates(X));
  }
  double eval(Balance *bal){
    arma::vec b = getBalance(bal);
    arma::mat v = b.t() * M * b;
    return v(0);
  }
};


// [[Rcpp::export]]
arma::mat find_principal_balance2_01(arma::mat X){
  int D = X.n_cols;

  arma::mat PB = arma::zeros(D, D-1);
  std::vector<PrincipalBalance2> SOLS;

  SOLS.push_back(PrincipalBalance2(X));
  SOLS.back().setOptimal();

  for(int l=0;l<D-1;l++){
    double vBestSolution = 0;
    int iBestSolution = -1;
    for(unsigned int i =0; i< SOLS.size();i++){
      double v = SOLS[i].EvaluateBalance::eval();
      if(v > vBestSolution){
        vBestSolution = v;
        iBestSolution = i;
      }
    }
    PB.col(l) = SOLS[iBestSolution].getBalance();

    int n = SOLS[iBestSolution].bal.N;
    int nL = SOLS[iBestSolution].bal.L_length;
    int nR = SOLS[iBestSolution].bal.R_length;
    if(n > nL + nR){
      SOLS.push_back(PrincipalBalance2(SOLS[iBestSolution].bal.top(), X));
      SOLS.back().setOptimal();
    }
    if(nL > 1){
      SOLS.push_back(PrincipalBalance2(SOLS[iBestSolution].bal.left(), X));
      SOLS.back().setOptimal();
    }
    if(nR > 1){
      SOLS.push_back(PrincipalBalance2(SOLS[iBestSolution].bal.right(), X));
      SOLS.back().setOptimal();
    }
    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();
    Rcpp::checkUserInterrupt();
  }
  return PB;
}
