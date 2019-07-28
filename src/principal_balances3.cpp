#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"

class PrincipalBalance3: public EvaluateBalance {
  arma::mat M;
  arma::mat PB;
  int npb = 0;
public:
  PrincipalBalance3(arma::mat X): EvaluateBalance(X.n_cols){
    M = cov(clr_coordinates(X));
  }
  PrincipalBalance3 (std::map<int,arma::uvec> nodes, arma::mat X, arma::mat PB_): EvaluateBalance(nodes, X.n_cols){
    M = cov(clr_coordinates(X));
    PB = PB_;
    npb = PB.n_cols;
    //Rcpp::Rcout << "Dim: " << npb << std::endl;
  }
  void setPB(arma::mat PB_){
    PB = PB_;
    npb = PB.n_cols;
  }
  double eval(Balance *bal){
    arma::vec b = getBalance(bal);
    int D = b.n_elem;
    arma::mat B = arma::mat(D, npb + 1);
    if(npb > 0){
      B.head_cols(npb) = PB;
    }
    B.col(npb) = b;

    arma::mat v = B.t() * M * B;
    // if(npb==2){
    //   Rcpp::Rcout << "B: "<< std::endl << B << std::endl;
    //   Rcpp::Rcout << "b: "<< std::endl << b << std::endl;
    //   // Rcpp::Rcout << "v: "<< std::endl << v << std::endl;
    //   // Rcpp::Rcout << "M: "<< std::endl << M << std::endl;
    //   Rcpp::Rcout << "det: " << arma::det(v) << std::endl;
    // }
    //
    return arma::det(v);
  }
};

//' @export
// [[Rcpp::export]]
arma::mat find_principal_balance3_01(arma::mat X){
  int D = X.n_cols;

  arma::mat PB = arma::zeros(D, D-1);
  std::vector<PrincipalBalance3> SOLS;

  SOLS.push_back(PrincipalBalance3(X));
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
      SOLS.push_back(PrincipalBalance3(SOLS[iBestSolution].bal.top(), X, PB.head_cols(l+1)));
      //SOLS.back().setOptimal();
    }
    if(nL > 1){
      SOLS.push_back(PrincipalBalance3(SOLS[iBestSolution].bal.left(), X, PB.head_cols(l+1)));
      //SOLS.back().setOptimal();
    }
    if(nR > 1){
      SOLS.push_back(PrincipalBalance3(SOLS[iBestSolution].bal.right(), X, PB.head_cols(l+1)));
      //SOLS.back().setOptimal();
    }
    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();
    for(int i=0; i < SOLS.size(); i++){
      SOLS[i].setPB(PB.head_cols(l+1));
      SOLS[i].setOptimal();
    }
    Rcpp::checkUserInterrupt();
  }
  return PB;
}
