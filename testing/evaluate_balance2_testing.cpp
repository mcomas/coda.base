#include <RcppArmadillo.h>
#include "../src/balance2.h"
#include "../src/coda.h"

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

class MaximumVariance: public EvaluateBalance2 {

  arma::mat lX;
public:
  MaximumVariance(arma::mat X){
    lX = log(X);
  }
  double evalBalance(arma::vec balance){
    return(var(lX * balance));
  }
};

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

class MaximumVariance3: public EvaluateBalance2 {

  arma::mat M;
public:
  MaximumVariance3(arma::mat X){
    M = cov(log(X));
    Rcout << M << std::endl;
  }
  // Need to be protected
  // double evalBalance(arma::vec balance){
  //   return( arma::accu((balance.t() * M * balance)) );
  // }
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
    Rcout << "Left:" << joinL.t() << ", size: " << joinL.size() << std::endl;
    Rcout << "Right:" << joinR.t() << ", size: " << joinR.size() << std::endl;
    for(int i=0;i<joinL.size();i++){
      variance += (nR/nL) * arma::accu(M(nodes[joinL[i]],nodes[joinL[i]]));
      for(int j=0;j<joinR.size();j++){
        variance += - 2 * arma::accu(M(nodes[joinL[i]],nodes[joinR[j]]));
      }
    }
    for(int i=0;i<joinR.size();i++){
      variance += (nL/nR) * arma::accu(M(nodes[joinR[i]],nodes[joinR[i]]));
    }
    Rcout << "Variance:" << variance / (nL+nR) << std::endl << std::endl;
    return variance / (nL+nR);
  }
};

// [[Rcpp::export]]
void eb_test_01a(arma::mat X, arma::vec V) {
  Balance2 balance = Balance2(X.n_cols);
  balance.addL(0);
  balance.addR(1);
  MaximumDotProduct ebalance1 = MaximumDotProduct(V);
  Rcpp::Rcout << "Current value:" << ebalance1.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance1.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance1.evalIfAddR(&balance, 2) << std::endl;
  Rcpp::Rcout << "Optimised minimum angle:" << balance.iterateLogContrast<MaximumDotProduct>(V, &ebalance1) << std::endl;
  balance.print();
}

// [[Rcpp::export]]
void eb_test_01b(arma::mat X, arma::vec V) {
  Balance2 balance = Balance2(X.n_cols);
  balance.addL(0);
  balance.addR(1);
  MaximumVariance ebalance2 = MaximumVariance(X);
  Rcpp::Rcout << "Current value:" << ebalance2.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance2.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance2.evalIfAddR(&balance, 2) << std::endl;
  Rcpp::Rcout << "Optimised maximum variance:" << balance.iterateLogContrast<MaximumVariance>(V, &ebalance2) << std::endl;
  balance.print();
}

// [[Rcpp::export]]
void eb_test_01c(arma::mat X, arma::vec V) {
  Balance2 balance = Balance2(X.n_cols);
  balance.addL(0);
  balance.addR(1);
  MaximumVariance3 ebalance2 = MaximumVariance3(X);
  Rcpp::Rcout << "Current value:" << ebalance2.eval(&balance) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to L:" << ebalance2.evalIfAddL(&balance, 2) << std::endl;
  Rcpp::Rcout << "Value if node 2 added to R:" << ebalance2.evalIfAddR(&balance, 2) << std::endl;
  Rcpp::Rcout << "Optimised maximum variance:" << balance.iterateLogContrast<MaximumVariance3>(V, &ebalance2) << std::endl;
  // balance.print();
}

// [[Rcpp::export]]
void eb_test_02(arma::mat X, arma::vec V) {
  Balance2 balance = Balance2(X.n_cols);
  balance.addL(0);
  balance.addR(1);
  MaximumDotProduct ebalance1 = MaximumDotProduct(V);
  MaximumVariance ebalance2 = MaximumVariance(X);
  MaximumVariance2 ebalance3 = MaximumVariance2(X);
  MaximumVariance3 ebalance4 = MaximumVariance3(X);
  clock_t t0 = clock();
  Rcpp::Rcout << "Optimised minimum angle:" << balance.iterateLogContrast<MaximumDotProduct>(V, &ebalance1) << std::endl;
  clock_t t1 = clock();
  Rcpp::Rcout << "Optimised maximum variance:" << balance.iterateLogContrast<MaximumVariance>(V, &ebalance2) << std::endl;
  clock_t t2 = clock();
  Rcpp::Rcout << "Optimised maximum variance2:" << balance.iterateLogContrast<MaximumVariance2>(V, &ebalance3) << std::endl;
  clock_t t3 = clock();
  Rcpp::Rcout << "Optimised maximum variance3:" << balance.iterateLogContrast<MaximumVariance3>(V, &ebalance4) << std::endl;
  clock_t t4 = clock();
  Rcout << "Timing angle:" << (double) (t1-t0) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance:" << (double) (t2-t1) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance2:" << (double) (t3-t2) / CLOCKS_PER_SEC * 1000.0 << std::endl;
  Rcout << "Timing variance3:" << (double) (t4-t3) / CLOCKS_PER_SEC * 1000.0 << std::endl;
}

/*** R
set.seed(2)
X = matrix(rlnorm(10*7), ncol = 7)
eb_test_01a(X, eigen(cov(X))$vector[,1,drop=FALSE])
eb_test_01b(X, eigen(cov(X))$vector[,1,drop=FALSE])
eb_test_01c(X, eigen(cov(X))$vector[,1,drop=FALSE])
*/
