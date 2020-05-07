#ifndef balance2_H
#define balance2_H

#include <RcppArmadillo.h>

using namespace Rcpp;
class Balance2;
class EvaluateBalance2;

class Balance2 {
public:
  int D;

  arma::uvec L;
  unsigned int L_length;
  arma::uvec R;
  unsigned int R_length;

  std::map<int,arma::uvec> nodes;
  int n_nodes;
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
  int size(){
    return(nodes.size());
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
    // Rcpp::Rcout << "uV" << uV;
    // Rcpp::Rcout << "uI" << uI;
    // New nodes pointers
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
  template <class EB>
  double iterateLogContrast(arma::vec LC, EB *ebalance){
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
    if(L_length > 0 & R_length > 0){
      Rcout << getBalance().t() << std::endl;
    }
  }

};

class EvaluateBalance2 {

protected:
  virtual double evalBalance(arma::vec balance){
    return(-1);
  }

public:
  EvaluateBalance2(){ }

  double eval(Balance2 *bal){
    return(evalBalance(bal->getBalance()));
  }
  double evalIfAddL(Balance2 *bal, unsigned iL){
    return(evalBalance(bal->getBalanceIfAddL(iL)));
  }
  double evalIfAddR(Balance2 *bal, unsigned iR){
    return(evalBalance(bal->getBalanceIfAddR(iR)));
  }
  double evalIfAdd(Balance2 *bal, arma::uvec uL, arma::uvec uR){
    return(evalBalance(bal->getBalanceIfAdd(uL, uR)));
  }
};

class MaximumDotProduct: public EvaluateBalance2 {

  arma::vec V;
public:
  MaximumDotProduct(arma::vec V0){
    V = V0;
  }
  double evalBalance(arma::vec balance){
    return(fabs(dot(balance, V)));
  }
};

// Better when N small
class MaximumVariance1: public EvaluateBalance2 {
  arma::mat lX;
public:
  MaximumVariance1(arma::mat X){
    lX = log(X);
  }
  double evalBalance(arma::vec balance){
    return(var(lX * balance));
  }
};

// Better when N>>0, D small
class MaximumVariance2: public EvaluateBalance2 {

  arma::mat M;
public:
  MaximumVariance2(arma::mat X){
    M = cov(log(X));
  }
  double evalBalance(arma::vec balance){
    return( arma::accu((balance.t() * M * balance)) );
  }
};

// Better when N>>0 or D>>0, small number of parts.
class MaximumVariance3: public EvaluateBalance2 {

  arma::mat M;
public:
  MaximumVariance3(arma::mat X){
    M = cov(log(X));
  }

  double eval(Balance2 *bal){
    std::map<int,arma::uvec> nodes = bal->nodes;
    arma::uvec L = bal->L;
    arma::uvec R = bal->R;

    double nL = bal->get_nL();
    double nR = bal->get_nR();
    double variance = 0;

    for(int i=0;i<bal->L_length;i++){
      for(int j=0;j<bal->L_length;j++){
        variance += (nR/nL) * arma::accu(M(nodes[L[i]],nodes[L[j]]));
      }
      for(int j=0;j<bal->R_length;j++){
        variance += - 2 * arma::accu(M(nodes[L[i]],nodes[R[j]]));
      }
    }
    for(int i=0;i<bal->R_length;i++){
      for(int j=0;j<bal->R_length;j++){
        variance += (nL/nR) * arma::accu(M(nodes[R[i]],nodes[R[j]]));
      }
    }
    return variance / (nL+nR);
  }
  double evalIfAddL(Balance2 *bal, unsigned iL){
    std::map<int,arma::uvec> nodes = bal->nodes;
    arma::uvec L = bal->L;
    arma::uvec R = bal->R;

    double nL = bal->get_nL() + nodes[iL].size();
    double nR = bal->get_nR();
    double variance = 0;

    // Node iL
    variance += (nR/nL) * arma::accu(M(nodes[iL],nodes[iL]));
    for(int j=0;j<bal->L_length;j++){
      variance += 2 * (nR/nL) * arma::accu(M(nodes[iL],nodes[L[j]]));
    }
    for(int j=0;j<bal->R_length;j++){
      variance += - 2 * arma::accu(M(nodes[iL],nodes[R[j]]));
    }

    for(int i=0;i<bal->L_length;i++){
      variance += (nR/nL) * arma::accu(M(nodes[L[i]],nodes[L[i]]));
      for(int j=0;j<bal->R_length;j++){
        variance += - 2 * arma::accu(M(nodes[L[i]],nodes[R[j]]));
      }
    }
    for(int i=0;i<bal->R_length;i++){
      variance += (nL/nR) * arma::accu(M(nodes[R[i]],nodes[R[i]]));
    }
    return variance / (nL+nR);
  }
  double evalIfAddR(Balance2 *bal, unsigned iR){
    std::map<int,arma::uvec> nodes = bal->nodes;
    arma::uvec L = bal->L;
    arma::uvec R = bal->R;

    double nL = bal->get_nL();
    double nR = bal->get_nR() + nodes[iR].size();
    double variance = 0;

    // Node iR
    variance += (nL/nR) * arma::accu(M(nodes[iR],nodes[iR]));
    for(int j=0;j<bal->L_length;j++){
      variance += - 2 * arma::accu(M(nodes[iR],nodes[L[j]]));
    }
    for(int j=0;j<bal->R_length;j++){
      variance += 2 * (nL/nR) * arma::accu(M(nodes[iR],nodes[R[j]]));
    }

    for(int i=0;i<bal->L_length;i++){
      variance += (nR/nL) * arma::accu(M(nodes[L[i]],nodes[L[i]]));
      for(int j=0;j<bal->R_length;j++){
        variance += - 2 * arma::accu(M(nodes[L[i]],nodes[R[j]]));
      }
    }
    for(int i=0;i<bal->R_length;i++){
      variance += (nL/nR) * arma::accu(M(nodes[R[i]],nodes[R[i]]));
    }
    return variance / (nL+nR);
  }
  double evalIfAdd(Balance2 *bal, arma::uvec uL, arma::uvec uR){
    std::map<int,arma::uvec> nodes = bal->nodes;

    double nL = bal->get_nL();
    for(unsigned int i = 0; i< uL.size(); i++) nL+=nodes[uL[i]].size();
    double nR = bal->get_nR();
    for(unsigned int i = 0; i< uR.size(); i++) nR+=nodes[uR[i]].size();
    double variance = 0;

    arma::uvec joinL = join_cols(uL, bal->L.head(bal->L_length));
    arma::uvec joinR = join_cols(uR, bal->R.head(bal->R_length));

    for(int i=0;i<joinL.size();i++){
      for(int j=0;j<joinL.size();j++){
        variance += (nR/nL) * arma::accu(M(nodes[joinL[i]],nodes[joinL[j]]));
      }
      for(int j=0;j<joinR.size();j++){
        variance += - 2 * arma::accu(M(nodes[joinL[i]],nodes[joinR[j]]));
      }
    }
    for(int i=0;i<joinR.size();i++){
      for(int j=0;j<joinR.size();j++){
        variance += (nL/nR) * arma::accu(M(nodes[joinR[i]],nodes[joinR[j]]));
      }
    }
    return variance / (nL+nR);
  }
};


#endif
