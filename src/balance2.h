#ifndef balance2_H
#define balance2_H

#include <RcppArmadillo.h>

using namespace Rcpp;

class Balance2 {
public:
  int D;

  arma::uvec L;
  unsigned int L_length;
  arma::uvec R;
  unsigned int R_length;

  std::map<int,arma::uvec> nodes;
  int n_nodes;

  Balance2(int D0, std::map<int,arma::uvec> nodes0){
    D = D0;

    nodes = nodes0;
    n_nodes = nodes0.size();

    L = arma::uvec(n_nodes);
    L_length = 0;
    R = arma::uvec(n_nodes);
    R_length = 0;
  }
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
  void approximateLogContrast(arma::vec LC){
    arma::vec V = arma::zeros(n_nodes);
    for(int i=0; i < n_nodes; i++){
      for(int j: nodes[i]){
        V(i) += LC(j);
      }
    }
    int imin = index_min(V);
    int imax = index_max(V);
    V(imin) = 0;
    V(imax) = 0;
    Rcout << "Min:" << imin << std::endl;
    Rcout << "Max:" << imax << std::endl;
    addL(imin);
    addR(imax);
    print();
    Rcout << "Min/Max:" << std::endl << getBalance() << std::endl;
    arma::uvec ord = sort_index(abs(V), "descend");
    for(int i = 0; i < n_nodes-2; i++){
      if(V(ord[i]) < 0){
        addL(ord[i]);
      }else{
        addR(ord[i]);
      }
      print();
      Rcout << "Node included:" << std::endl << getBalance() << std::endl;
    }
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
  void print(){
    Rcout << "Elements: ";
    for(unsigned int i=0; i<n_nodes;i++){
      Rcout << "{";
      for(int j = 0; j<nodes[i].n_elem;j++) Rcout << " " << nodes[i][j];
      Rcout << " } ";
    }
    Rcout << "L:";
    for(unsigned int i=0; i < L_length; i++){
      Rcout << " " << L[i];
    }
    Rcout << ", R:";
    for(unsigned int i=0; i < R_length; i++){
      Rcout << " " << R[i];
    }
    Rcout << std::endl;
  }

};

#endif
