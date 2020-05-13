#ifndef Balance_H
#define Balance_H

#include <RcppArmadillo.h>
#include "balance_evaluate.h"
#include "balance_optimal.h"

template <class EB>
class Balance {
public:
  int D;

  arma::uvec L, R;
  unsigned int L_length, R_length;

  std::map<int,arma::uvec> nodes;
  int n_nodes;

  EB ebalance = EB();
  double ebalance_value = 0;

  Balance(int D0){
    D = D0;
    for(int i=0;i<D;i++){
      nodes[i] = arma::uvec(1);
      nodes[i][0] = i;
    }
    n_nodes = D;

    L = arma::uvec(n_nodes);
    L_length = 0;
    R = arma::uvec(n_nodes);
    R_length = 0;
  }
  Balance(int D0, std::map<int,arma::uvec> nodes0){
    D = D0;

    nodes = nodes0;
    n_nodes = nodes0.size();

    L = arma::uvec(n_nodes);
    L_length = 0;
    R = arma::uvec(n_nodes);
    R_length = 0;
  }
  void setEvaluator(EB ebalance0){
    ebalance = ebalance0;
  }
  void set(arma::uvec uL, arma::uvec uR){
    L_length = uL.size();
    L.head(L_length) = uL;

    R_length = uR.size();
    R.head(R_length) = uR;

    ebalance_value = ebalance.eval(L, R, L_length, R_length);
  }
  void setWithExhaustiveSearch(){
    int N = nodes.size();
    std::vector<arma::uvec> P = std::vector<arma::uvec>(3);
    P[0] = arma::uvec(N);
    P[1] = arma::uvec(N);
    P[2] = arma::uvec(N);
    int p[3];
    p[0] = N-2;
    p[1] = 1;
    p[2] = 1;
    for(int i = 0;i < N-2; i++) P[0][i] = i;
    P[1][0] = N-2;
    P[2][0] = N-1;
    arma::uvec I = arma::uvec(N);
    for(int i = 0;i < N-2; i++) I(i) = i;
    I(N-2) = 0;
    I(N-1) = 0;

    arma::uvec A = arma::uvec(N+1);
    A.fill(0);
    A[N-1] = 1;
    A[N] = 2;
    f<EB>(3, N, 0, I, A, P, p, ebalance);

    set(ebalance.bestL, ebalance.bestR);

  }
  void setWithLogContrast(arma::vec V){

    int imin = index_min(V);
    int imax = index_max(V);

    V(imin) = 0;
    V(imax) = 0;
    arma::uvec ord = sort_index(abs(V), "descend");
    arma::uvec uL(ord.size()), uR(ord.size());
    uL[0] = imin; uR[0] = imax;
    unsigned l = 1, r = 1;

    ebalance.eval(uL, uR, l, r);
    for(int i = 0; i < n_nodes-2; i++){
      if(V(ord[i]) < 0) uL(l++) = ord[i];
      else uR(r++) = ord[i];

      ebalance.eval(uL, uR, l, r);
    }
    set(ebalance.bestL, ebalance.bestR);
  }

  void setWithLogContrastForceBranch(arma::vec V, arma::uvec forced){

    double vforced = V(forced[0]);
    V(forced[0]) = 0;
    for(int i=1; i < forced.size(); i++){
      vforced += V(forced[i]);
      V(forced[i]) = 0;
    }
    int imax = index_max(abs(V));
    arma::uvec ord;
    if(V(imax) > 0){
       ord = sort_index(V, "descend");
    }else{
      ord = sort_index(V, "ascend");
    }
    unsigned pos = 0;
    arma::uvec uR(ord.size());
    while(V(pos) != 0){
      uR[pos] = ord[pos];
      pos++;
      ebalance.eval(forced, uR, forced.size(), pos);
    }
    set(ebalance.bestL, ebalance.bestR);

  }


  arma::vec getBalance(){
    double nL = 0, nR = 0;
    for(unsigned int i = 0; i< L_length; i++) nL+=nodes[L[i]].size();
    for(unsigned int i = 0; i< R_length; i++) nR+=nodes[R[i]].size();

    arma::vec b = arma::zeros(D);
    double l_v = -1/nL * sqrt(nL*nR/(nL+nR));
    double r_v = +1/nR * sqrt(nL*nR/(nL+nR));
    for(unsigned int i = 0; i< L_length; i++) b(nodes[L[i]]).fill(l_v);
    for(unsigned int i = 0; i< R_length; i++) b(nodes[R[i]]).fill(r_v);
    return(b);
  }
  double eval(){
    return(ebalance_value);
  }


  Balance<EB> top(){
    arma::uvec O = arma::zeros<arma::uvec>(nodes.size());
    O(L.head(L_length)).fill(1);
    O(R.head(R_length)).fill(1);

    arma::uvec uV = find(O == 1);
    arma::uvec uI = find(O == 0);

    int nI = uI.n_elem;
    std::map<int,arma::uvec> node0;
    for(int i = 0; i<nI; i++){
      node0[i] = nodes[uI[i]];
    }

    int nV = 0;
    for(unsigned int i = 0; i< L_length; i++) nV+=nodes[L[i]].size();
    for(unsigned int i = 0; i< R_length; i++) nV+=nodes[R[i]].size();
    arma::uvec V(nV);
    int k = 0;
    for(unsigned int i= 0; i < L_length; i++)
      for(unsigned int j= 0; j < nodes[L[i]].n_elem; j++, k++)
        V(k) = nodes[L[i]][j];
    for(unsigned int i=0; i < R_length; i++)
      for(unsigned int j=0; j < nodes[R[i]].n_elem; j++, k++)
        V(k) = nodes[R[i]][j];
    node0[nI] = V;

    return(Balance<EB>(D, node0));
  }

  Balance<EB> left(){
    arma::uvec uL = L.head(L_length);
    arma::uvec uR = R.head(R_length);
    std::map<int,arma::uvec> node0;
    for(unsigned int i = 0; i<uL.n_elem; i++){
      node0[i] = arma::uvec(nodes[uL[i]]);
    }
    return(Balance<EB>(D, node0));
  }

  Balance<EB> right(){
    arma::uvec uL = L.head(L_length);
    arma::uvec uR = R.head(R_length);
    std::map<int,arma::uvec> node0;
    for(unsigned int i = 0; i<uR.n_elem; i++){
      node0[i] = arma::uvec(nodes[uR[i]]);
    }
    return(Balance<EB>(D, node0));
  }

  void print(){
    Rcpp::Rcout << "Elements: ";
    for(unsigned int i=0; i<n_nodes;i++){
      Rcpp::Rcout << "{";
      for(int j = 0; j<nodes[i].n_elem;j++) Rcpp::Rcout << " " << nodes[i][j];
      Rcpp::Rcout << " } ";
    }
    Rcpp::Rcout << "L:";
    for(unsigned int i=0; i < L_length; i++){
      Rcpp::Rcout << " " << L[i];
    }
    Rcpp::Rcout << ", R:";
    for(unsigned int i=0; i < R_length; i++){
      Rcpp::Rcout << " " << R[i];
    }
    Rcpp::Rcout << std::endl;
    if(L_length > 0 & R_length > 0){
      Rcpp::Rcout << getBalance().t() << "Value = " << ebalance.eval(L, R, L_length, R_length) << std::endl;
    }else{
      Rcpp::Rcout << "Balance is not defined" << std::endl;
    }
  }
};

#endif
