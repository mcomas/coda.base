#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"

class PrincipalBalance: public EvaluateBalance {
  arma::mat M;
public:
  PrincipalBalance(arma::mat X): EvaluateBalance(X.n_cols){
    int n = bal.nodes.size();
    arma::mat M_all = cov(log(X));
    M = arma::zeros<arma::mat>(n, n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        M(i,j) = arma::accu(M_all(bal.nodes[i], bal.nodes[j]));
      }
    }
  }
  PrincipalBalance (std::map<int,arma::uvec> nodes, arma::mat X): EvaluateBalance(nodes, X.n_cols){
    int n = bal.nodes.size();
    arma::mat M_all = cov(log(X));
    M = arma::zeros<arma::mat>(n, n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        M(i,j) = arma::accu(M_all(bal.nodes[i], bal.nodes[j]));
      }
    }
  }
  double eval(Balance *bal){
    double nL = bal->get_nL(); //bal->L_length; //
    double nR = bal->get_nR(); //bal->R_length; //
    arma::uvec uL = bal->getL();
    arma::uvec uR = bal->getR();
    double sL = (nR/nL) * arma::accu(M(uL,uL));
    double sR = (nL/nR) * arma::accu(M(uR,uR));
    double sC = - 2*arma::accu(M(uR,uL));

    double variance = (sL + sR + sC) / (nL+nR);
    return variance;
  }
};

// [[Rcpp::export]]
arma::mat find_principal_balance_01(arma::mat X){
  int D = X.n_cols;

  arma::mat PB = arma::zeros(D, D-1);
  std::vector<PrincipalBalance> SOLS;

  SOLS.push_back(PrincipalBalance(X));
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
      SOLS.push_back(PrincipalBalance(SOLS[iBestSolution].bal.top(), X));
      SOLS.back().setOptimal();
    }
    if(nL > 1){
      SOLS.push_back(PrincipalBalance(SOLS[iBestSolution].bal.left(), X));
      SOLS.back().setOptimal();
    }
    if(nR > 1){
      SOLS.push_back(PrincipalBalance(SOLS[iBestSolution].bal.right(), X));
      SOLS.back().setOptimal();
    }
    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();
    Rcpp::checkUserInterrupt();
  }
  return PB;
}

// [[Rcpp::export]]
arma::mat find_principal_balance_02(arma::mat X){
  int D = X.n_cols;

  arma::mat PB = arma::zeros(D, D-1);
  std::vector<PrincipalBalance> SOLS;

  SOLS.push_back(PrincipalBalance(X));
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
      SOLS.push_back(PrincipalBalance(SOLS[iBestSolution].bal.top(), X));
      SOLS.back().bal.rnd_init();
      SOLS.back().setLocalSearch();
    }
    if(nL > 1){
      SOLS.push_back(PrincipalBalance(SOLS[iBestSolution].bal.left(), X));
      SOLS.back().bal.rnd_init();
      SOLS.back().setLocalSearch();
    }
    if(nR > 1){
      SOLS.push_back(PrincipalBalance(SOLS[iBestSolution].bal.right(), X));
      SOLS.back().bal.rnd_init();
      SOLS.back().setLocalSearch();
    }
    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();
    Rcpp::checkUserInterrupt();
  }
  return PB;
}
