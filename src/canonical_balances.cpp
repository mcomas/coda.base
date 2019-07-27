#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"

class CanonicalBalance: public EvaluateBalance {
  arma::mat SX, SY, SXY;
public:
  CanonicalBalance(arma::mat X, arma::mat Y): EvaluateBalance(X.n_cols){
    SX = cov(clr_coordinates(X));
    SY = cov(Y);
    SXY = cov(clr_coordinates(X), Y);
  }
  CanonicalBalance (std::map<int,arma::uvec> nodes, arma::mat X, arma::mat Y): EvaluateBalance(nodes, X.n_cols){
    SX = cov(clr_coordinates(X));
    SY = cov(Y);
    SXY = cov(clr_coordinates(X), Y);
  }
  double eval(Balance *bal){
    //b = (chol2inv(chol(S22)) %*% t(S12) %*% a)
    //(t(a) %*% S12 %*% b)^2 / ((t(a) %*% S11 %*% a) * (t(b) %*% S22 %*% b))
    arma::vec a = getBalance(bal);
    arma::vec b = arma::pinv(SY) * SXY.t() * a;
    arma::mat numerator = a.t() * SXY * b;
    arma::mat denominator = (a.t() * SX * a) * (b.t() * SY * b);
    return numerator(0)*numerator(0)/denominator(0);
  }
};

//' @export
// [[Rcpp::export]]
arma::mat find_canonical_balance_01(arma::mat X, arma::mat Y){
  int D = X.n_cols;

  arma::mat PB = arma::zeros(D, D-1);
  std::vector<CanonicalBalance> SOLS;

  SOLS.push_back(CanonicalBalance(X, Y));
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
      SOLS.push_back(CanonicalBalance(SOLS[iBestSolution].bal.top(), X, Y));
      SOLS.back().setOptimal();
    }
    if(nL > 1){
      SOLS.push_back(CanonicalBalance(SOLS[iBestSolution].bal.left(), X, Y));
      SOLS.back().setOptimal();
    }
    if(nR > 1){
      SOLS.push_back(CanonicalBalance(SOLS[iBestSolution].bal.right(), X, Y));
      SOLS.back().setOptimal();
    }
    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();
    Rcpp::checkUserInterrupt();
  }
  return PB;
}

//' @export
// [[Rcpp::export]]
arma::mat find_canonical_balance_02(arma::mat X){
  int D = X.n_cols;

  arma::mat PB = arma::zeros(D, D-1);
  std::vector<CanonicalBalance> SOLS;

  SOLS.push_back(CanonicalBalance(X, X));
  SOLS.back().bal.rnd_init();
  SOLS.back().setLocalSearch();

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
      SOLS.push_back(CanonicalBalance(SOLS[iBestSolution].bal.top(), X, X));
      SOLS.back().bal.rnd_init();
      SOLS.back().setLocalSearch();
    }
    if(nL > 1){
      SOLS.push_back(CanonicalBalance(SOLS[iBestSolution].bal.left(), X, X));
      SOLS.back().bal.rnd_init();
      SOLS.back().setLocalSearch();
    }
    if(nR > 1){
      SOLS.push_back(CanonicalBalance(SOLS[iBestSolution].bal.right(), X, X));
      SOLS.back().bal.rnd_init();
      SOLS.back().setLocalSearch();
    }
    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();
    Rcpp::checkUserInterrupt();
  }
  return PB;
}
