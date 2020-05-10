#ifndef evaluate_balance2_H
#define evaluate_balance2_H

#include <RcppArmadillo.h>
#include "balance2.h"

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
