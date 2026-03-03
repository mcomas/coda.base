#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"

const double SQR2DIV2 = 0.70710678118;

void optimise(Balance<MaximumVariance>& balance, const arma::mat& X) {
  MaximumVariance ebalance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
}

// [[Rcpp::export]]
arma::mat find_PB(const arma::mat& X) {
  const int K = X.n_cols;

  arma::mat pb_mat(K, K - 1);
  std::vector<Balance<MaximumVariance>> sols;

  Balance<MaximumVariance> balance(X.n_cols);
  optimise(balance, X);
  sols.push_back(balance);

  for (int l = 0; l < K - 1; ++l) {
    double vBestSolution = -std::numeric_limits<double>::infinity();
    int iBestSolution = 0;

    for (unsigned int i = 0; i < sols.size(); ++i) {
      const double v = sols[i].eval();
      if (v > vBestSolution) {
        vBestSolution = v;
        iBestSolution = i;
      }
    }

    pb_mat.col(l) = sols[iBestSolution].getBalance();

    Balance<MaximumVariance> top = sols[iBestSolution].top();
    Balance<MaximumVariance> left = sols[iBestSolution].left();
    Balance<MaximumVariance> right = sols[iBestSolution].right();

    if (top.nodes.size() > 1) {
      optimise(top, X);
      sols.push_back(top);
    }
    if (left.nodes.size() > 1) {
      optimise(left, X);
      sols.push_back(left);
    }
    if (right.nodes.size() > 1) {
      optimise(right, X);
      sols.push_back(right);
    }

    sols[iBestSolution] = sols.back();
    sols.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return pb_mat;
}

// [[Rcpp::export]]
arma::vec get_balance_using_pc(arma::mat& X, bool angle = false){
  unsigned D = X.n_cols;
  if(D == 2){
    arma::vec balance = {+SQR2DIV2, -SQR2DIV2};
    return(balance);
  }else{
    arma::vec eigval;
    arma::mat eigvec;

    arma::mat S = cov(clr_coordinates(X));
    eig_sym( eigval, eigvec, S);

    arma::vec V = eigvec.tail_cols(1);
    //Rcpp::Rcout << V.t();

    unsigned imin = index_min(V);
    unsigned imax = index_max(V);

    V(imin) = 0;
    V(imax) = 0;
    arma::uvec ord = sort_index(abs(V), "descend");
    //Rcpp::Rcout << ord;
    arma::uvec uL(D), uR(D);
    uL[0] = imin; uR[0] = imax;
    unsigned l = 1, r = 1;
    arma::vec balance = arma::zeros(D);
    balance[imin] = -SQR2DIV2;
    balance[imax] = +SQR2DIV2;
    double bestScore = as_scalar(balance.t() * S * balance);
    if(angle)  bestScore = fabs(dot(eigvec.tail_cols(1), balance));

    double bestR = 1, bestL = 1;
    for(unsigned i = 0; i < D-2; i++){
      if(V(ord[i]) < 0) uL(l++) = ord[i];
      else uR(r++) = ord[i];

      balance(uL.head(l)).fill(-1.0/l * sqrt((double)l*r/(l+r)));
      balance(uR.head(r)).fill(+1.0/r * sqrt((double)l*r/(l+r)));

      double score = as_scalar(balance.t() * S * balance);
      if(angle) score = fabs(dot(eigvec.tail_cols(1), balance));

      if(score > bestScore){
        bestScore = score;
        bestL = l;
        bestR = r;
      }
    }

    balance.fill(0);
    balance(uL.head(bestL)).fill(-1.0/bestL * sqrt((double)bestL*bestR/(bestL+bestR)));
    balance(uR.head(bestR)).fill(+1.0/bestR * sqrt((double)bestL*bestR/(bestL+bestR)));
    return(balance);

  }
}
