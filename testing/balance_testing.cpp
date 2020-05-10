#include <RcppArmadillo.h>
#include "../src/balance.h"
#include "../src/balance_evaluate.h"
#include "../src/coda.h"

void optimise(Balance<MaximumVariance>& balance, arma::mat &X){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
}

// [[Rcpp::export]]
void testing_01(arma::mat X) {

  Balance<MaximumVariance> balance = Balance<MaximumVariance>(X.n_cols);
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
  balance.print();

  Balance<MaximumVariance> top = balance.top();
  Balance<MaximumVariance> left = balance.left();
  Balance<MaximumVariance> right = balance.right();

  if(top.nodes.size() > 1){
    optimise(top, X);
    top.print();
  }
  if(left.nodes.size() > 1){
    optimise(left, X);
    left.print();
  }
  if(right.nodes.size() > 1){
    optimise(right, X);
    right.print();
  }

}



/*** R
set.seed(2)
X = matrix(rlnorm(10*7), ncol = 7)

testing_01(X)
*/
