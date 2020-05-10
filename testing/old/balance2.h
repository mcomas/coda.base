#ifndef balance2_H
#define balance2_H

#include <RcppArmadillo.h>
#include "evaluate_balance2.h"
#include "exhaustive.h"

template <class EB>
class Balance2 {
public:
  int D;

  arma::uvec L, R;
  unsigned int L_length, R_length;

  std::map<int,arma::uvec> nodes;
  int n_nodes;

  EB ebalance;

  Balance2(int D0){
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
  Balance2(int D0, std::map<int,arma::uvec> nodes0){
    D = D0;

    nodes = nodes0;
    n_nodes = nodes0.size();

    L = arma::uvec(n_nodes);
    L_length = 0;
    R = arma::uvec(n_nodes);
    R_length = 0;
  }
  // int size(){
  //   return(nodes.size());
  // }
  void addL(unsigned I){
    L(L_length) = I;
    L_length++;
  }
  void addR(unsigned I){
    R(R_length) = I;
    R_length++;
  }
  void setL(arma::uvec uL){
    L_length = uL.size();
    for(int i=0;i< L_length;i++) L(i) = uL(i);
  }
  void setR(arma::uvec uR){
    R_length = uR.size();
    for(int i=0;i< R_length;i++) R(i) = uR(i);
  }
  void setL(int iL){
    L_length = 1;
    L(0) = iL;
  }
  void setR(int iR){
    R_length = 1;
    R(0) = iR;
  }
  // Number of parts in left side
  int get_nL(){
    double nL = 0;
    for(unsigned int i=0; i<L_length;nL+=nodes[L[i++]].n_elem);
    return(nL);
  }
  // Number of parts in right side
  int get_nR(){
    double nR = 0;
    for(unsigned int i=0; i<R_length;nR+=nodes[R[i++]].n_elem);
    return(nR);
  }
  Balance2 top(){
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

    int nV = get_nL() + get_nR();
    arma::uvec V(nV);
    int k = 0;
    for(unsigned int i= 0; i < L_length; i++)
      for(unsigned int j= 0; j < nodes[L[i]].n_elem; j++, k++)
        V(k) = nodes[L[i]][j];
    for(unsigned int i=0; i < R_length; i++)
      for(unsigned int j=0; j < nodes[R[i]].n_elem; j++, k++)
        V(k) = nodes[R[i]][j];
    node0[nI] = V;
    Balance2 balance(D, node0);
    return(balance);
  }
  Balance2 left(){
    arma::uvec uL = L.head(L_length);
    arma::uvec uR = R.head(R_length);
    std::map<int,arma::uvec> node0;
    for(unsigned int i = 0; i<uL.n_elem; i++){
      node0[i] = arma::uvec(nodes[uL[i]]);
    }
    Balance2 balance(D, node0);
    return(balance);
  }
  Balance2 right(){
    arma::uvec uL = L.head(L_length);
    arma::uvec uR = R.head(R_length);
    std::map<int,arma::uvec> node0;
    for(unsigned int i = 0; i<uR.n_elem; i++){
      node0[i] = arma::uvec(nodes[uR[i]]);
    }
    Balance2 balance(D, node0);
    return(balance);
  }
  arma::vec getBalance(){
    double nL = 0, nR = 0;
    for(unsigned int i = 0; i< L_length; i++) nL+=nodes[L[i]].size();
    for(unsigned int i = 0; i< R_length; i++) nR+=nodes[R[i]].size();

    arma::vec b = arma::zeros(D);
    for(unsigned int i = 0; i< L_length; i++) b(nodes[L[i]]).fill(-1/nL * sqrt(nL*nR/(nL+nR)));
    for(unsigned int i = 0; i< R_length; i++) b(nodes[R[i]]).fill(+1/nR * sqrt(nL*nR/(nL+nR)));
    return(b);
  }
  arma::vec getBalanceIfAdd(arma::uvec uL, arma::uvec uR){
    double nL = 0, nR = 0;
    for(unsigned int i = 0; i< L_length; i++) nL+=nodes[L[i]].size();
    for(unsigned int i = 0; i< uL.size(); i++) nL+=nodes[uL[i]].size();
    for(unsigned int i = 0; i< R_length; i++) nR+=nodes[R[i]].size();
    for(unsigned int i = 0; i< uR.size(); i++) nR+=nodes[uR[i]].size();

    arma::vec b = arma::zeros(D);
    for(unsigned int i = 0; i< L_length; i++) b(nodes[L[i]]).fill(-1/nL * sqrt(nL*nR/(nL+nR)));
    for(unsigned int i = 0; i< uL.size(); i++) b(nodes[uL[i]]).fill(-1/nL * sqrt(nL*nR/(nL+nR)));
    for(unsigned int i = 0; i< R_length; i++) b(nodes[R[i]]).fill(+1/nR * sqrt(nL*nR/(nL+nR)));
    for(unsigned int i = 0; i< uR.size(); i++) b(nodes[uR[i]]).fill(+1/nR * sqrt(nL*nR/(nL+nR)));
    return(b);
  }
  arma::vec getBalanceIfAddL(unsigned iL){
    double nL = nodes[iL].size(), nR = 0;
    for(unsigned int i = 0; i< L_length; i++) nL+=nodes[L[i]].size();
    for(unsigned int i = 0; i< R_length; i++) nR+=nodes[R[i]].size();

    arma::vec b = arma::zeros(D);
    b(nodes[iL]).fill(-1/nL * sqrt(nL*nR/(nL+nR)));
    for(unsigned int i = 0; i< L_length; i++) b(nodes[L[i]]).fill(-1/nL * sqrt(nL*nR/(nL+nR)));
    for(unsigned int i = 0; i< R_length; i++) b(nodes[R[i]]).fill(+1/nR * sqrt(nL*nR/(nL+nR)));
    return(b);
  }
  arma::vec getBalanceIfAddR(unsigned iR){
    double nL = 0, nR = nodes[iR].size();
    for(unsigned int i = 0; i< L_length; i++) nL+=nodes[L[i]].size();
    for(unsigned int i = 0; i< R_length; i++) nR+=nodes[R[i]].size();

    arma::vec b = arma::zeros(D);
    b(nodes[iR]).fill(+1/nR * sqrt(nL*nR/(nL+nR)));
    for(unsigned int i = 0; i< L_length; i++) b(nodes[L[i]]).fill(-1/nL * sqrt(nL*nR/(nL+nR)));
    for(unsigned int i = 0; i< R_length; i++) b(nodes[R[i]]).fill(+1/nR * sqrt(nL*nR/(nL+nR)));
    return(b);
  }
  double approximateLogContrast(arma::vec LC, arma::mat *lX = NULL){
    arma::vec V = arma::zeros(n_nodes);
    for(int i=0; i < n_nodes; i++){
      for(int j: nodes[i]){
        V(i) += LC(j);
      }
    }

    int imin = index_min(V);
    int imax = index_max(V);

    setL(imin);
    setR(imax);

    arma::vec bal = getBalance();

    V(imin) = 0;
    V(imax) = 0;
    arma::uvec ord = sort_index(abs(V), "descend");
    arma::uvec uL(ord.size()), uR(ord.size());
    int nL = 0, nR = 0;
    int mL = 0, mR = 0;
    double score_max = dot(bal,LC);
    if(lX){
      score_max = var((*lX) * bal);
    }
    for(int i = 0; i < n_nodes-2; i++){

      if(V(ord[i]) < 0) uL(nL++) = ord[i];
      else uR(nR++) = ord[i];

      arma::vec bal = getBalanceIfAdd(uL.head(nL), uR.head(nR));
      double score = fabs(dot(bal,LC));
      if(lX){
        score = var((*lX) * bal);
      }
      if(score > score_max){
        score_max = score;
        mL = nL;
        mR = nR;
      }
    }

    for(int i=0; i < mL; i++) addL(uL(i));
    for(int i=0; i < mR; i++) addR(uR(i));

    return(score_max);
  }
  double iterateLogContrast(arma::vec LC){
    arma::vec V = arma::zeros(n_nodes);
    for(int i=0; i < n_nodes; i++){
      for(int j: nodes[i]){
        V(i) += LC(j);
      }
    }

    int imin = index_min(V);
    int imax = index_max(V);

    setL(imin);
    setR(imax);

    V(imin) = 0;
    V(imax) = 0;
    arma::uvec ord = sort_index(abs(V), "descend");
    arma::uvec uL(ord.size()), uR(ord.size());
    int nL = 0, nR = 0;
    int mL = 0, mR = 0;
    double score_max = ebalance->eval(this);
    for(int i = 0; i < n_nodes-2; i++){
      if(V(ord[i]) < 0) uL(nL++) = ord[i];
      else uR(nR++) = ord[i];

      double score = ebalance->evalIfAdd(this, uL.head(nL), uR.head(nR));
      if(score > score_max){
        score_max = score;
        mL = nL;
        mR = nR;
      }
    }
    for(int i=0; i < mL; i++) addL(uL(i));
    for(int i=0; i < mR; i++) addR(uR(i));

    return(score_max);
  }
  double exhaustiveIteration(){
    std::vector<arma::uvec> P = std::vector<arma::uvec>(3);
    P[0] = arma::uvec(D);
    P[1] = arma::uvec(D);
    P[2] = arma::uvec(D);
    int p[3];
    p[0] = D-2;
    p[1] = 1;
    p[2] = 1;
    for(int i = 0;i < D-2; i++) P[0][i] = i;
    P[1][0] = D-2;
    P[2][0] = D-1;
    arma::uvec I = arma::uvec(D);
    for(int i = 0;i < D-2; i++) I(i) = i;
    I(D-2) = 0;
    I(D-1) = 0;

    arma::uvec A = arma::uvec(D+1);
    A.fill(0);
    A[D-1] = 1;
    A[D] = 2;
    f(3, D, 0, I, A, P, p, ebalance);
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
      Rcpp::Rcout << getBalance().t() << std::endl;
    }
  }

};

#endif
