#ifndef balance2_H
#define balance2_H

#include <RcppArmadillo.h>

using namespace Rcpp;

class Balance2 {
public:
  int N;

  arma::uvec L;
  unsigned int L_length;
  arma::uvec R;
  unsigned int R_length;

  std::map<int,arma::uvec> nodes;

  Balance2(std::map<int,arma::uvec> nodes0){

    nodes = nodes0;
    N = nodes.size();
    L = arma::uvec(N);
    L_length = 0;
    R = arma::uvec(N);
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
  void approximateLogContrast(arma::vec LC){
    arma::vec V = arma::zeros(N);
    for(int i=0; i <N; i++){
      for(int j: nodes[i]){
        V(i) += LC(j);
      }
    }
    Rcout << LC << std::endl<<  V << std::endl;
    int imin = index_min(V);
    int imax = index_max(V);
    Rcout << "Min:" << imin << std::endl;
    Rcout << "Max:" << imax << std::endl;
    V(imin) = 0;
    V(imax) = 0;
    for(int i: nodes[imin]) addL(i);
    for(int i: nodes[imax]) addR(i);
    Rcout << sort_index(abs(V), "descend");
  }
  void print(){
    Rcout << "Elements: ";
    for(unsigned int i=0; i<nodes.size();i++){
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

std::map<int,arma::uvec> default_node(int size){
  std::map<int,arma::uvec> node;
  for(int i=0;i<size;i++){
    node[i] = arma::uvec(1);
    node[i][0] = i;
  }
  return(node);
}

#endif
