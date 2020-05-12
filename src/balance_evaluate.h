#ifndef balance_evaluate_H
#define balance_evaluate_H

#include <RcppArmadillo.h>

class EvaluateBalance {

public:
  EvaluateBalance(){

  }

  virtual double eval(arma::uvec& L, arma::uvec& R, int l, int r){
    return(-1);
  }

};

class MaximumVariance: public EvaluateBalance {
  arma::mat M;
  arma::vec N;
  std::map<int,arma::uvec> nodes;

  double bestScore = -1;

public:
  arma::uvec bestL, bestR;

  MaximumVariance(){
  }
  MaximumVariance(std::map<int,arma::uvec>& nodes0, arma::mat& X){

    nodes = nodes0;
    arma::mat S = cov(log(X));
    M = arma::mat(nodes0.size(), nodes0.size());

    N = arma::vec(nodes0.size());
    for(int i = 0; i < nodes0.size(); i++){
      N(i) = nodes0[i].n_elem;
      for(int j = 0; j < nodes0.size(); j++){
        M(i,j) = arma::accu(S(nodes0[i], nodes0[j]));
      }
    }
  }
  void init(){
    bestScore = -1;
  }
  double eval(arma::uvec& L, arma::uvec& R, int l, int r){
    double nL = 0;
    for(unsigned int i=0; i<l;nL+=N[L[i++]]);

    double nR = 0;
    for(unsigned int i=0; i<r;nR+=N[R[i++]]);

    double variance = 0;
    for(int i=0;i<l;i++){
      for(int j=0;j<l;j++){
        variance += (nR/nL) * M(L[i],L[j]);
      }
      for(int j=0;j<r;j++){
        variance += - 2 * M(L[i],R[j]);
      }
    }
    for(int i=0;i<r;i++){
      for(int j=0;j<r;j++){
        variance += (nL/nR) * M(R[i],R[j]);
      }
    }
    variance = variance / (nL+nR);
    if(variance > bestScore){
      bestScore = variance;
      bestL = arma::uvec(L.head(l));
      bestR = arma::uvec(R.head(r));
    }
    return variance;
  }
//
//   double eval(arma::uvec& L, arma::uvec& R){
//     return(eval(L, R, L.size(), R.size()));
//   }

};

// class MaximumDotProduct: public EvaluateBalance {
//
//   std::map<int,arma::uvec>& nodes;
//   arma::vec V;
// public:
//   MaximumDotProduct(std::map<int,arma::uvec>& nodes0, arma::vec& V0){
//     nodes = nodes0;
//     V = V0;
//   }
//   double eval(arma::uvec& L, arma::uvec& R, int l, int r){
//
//     double nL = 0, nR = 0;
//     for(unsigned int i = 0; i< l; i++) nL+=nodes[L[i]].size();
//     for(unsigned int i = 0; i< r; i++) nR+=nodes[R[i]].size();
//
//     arma::vec b = arma::zeros(D);
//     for(unsigned int i = 0; i< l; i++) b(nodes[L[i]]).fill(-1/nL * sqrt(nL*nR/(nL+nR)));
//     for(unsigned int i = 0; i< r; i++) b(nodes[R[i]]).fill(+1/nR * sqrt(nL*nR/(nL+nR)));
//
//     return(fabs(dot(b, V)));
//   }
// };

#endif
