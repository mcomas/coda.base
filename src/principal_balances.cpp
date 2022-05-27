// #ifndef balance_optimal_H
// #define balance_optimal_H

#include <RcppArmadillo.h>
#include "balance.h"
#include "coda.h"

const double SQR2DIV2 = 0.70710678118;

void optimise(Balance<MaximumVariance>& balance, arma::mat& X){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  balance.setWithExhaustiveSearch();
}

void optimise_balance_using_pc(Balance<MaximumVariance>& balance, arma::mat& X){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  if(balance.nodes.size() == 2){
    arma::uvec uL(1), uR(1);
    uL[0] = 0; uR[0] = 1;
    ebalance.eval(uL, uR, 1, 1);
    balance.set(uL, uR);
  }else{
    arma::mat Xsub = arma::mat(X.n_rows, balance.nodes.size());
    for(int i = 0; i < Xsub.n_cols; i++){
      Xsub.col(i) = X.col(balance.nodes[i][0]);
      for(int j = 1; j < balance.nodes[i].n_elem; j++){
        Xsub.col(i) += X.col(balance.nodes[i][j]);
      }
    }
    arma::mat U, V;
    arma::vec s;

    arma::svd_econ(U, s, V, clr_coordinates(Xsub));
    // Rcpp::Rcout << V.col(0).t();
    balance.setWithLogContrast(V.col(0));

  }
}

// [[Rcpp::export]]
arma::vec get_balance_using_pc(arma::mat& X){
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
    // Rcpp::Rcout << bestScore << std::endl;
    //double bestScore = fabs(dot(eigvec.tail_cols(1), balance));
    double bestR = 1, bestL = 1;
    for(unsigned i = 0; i < D-2; i++){
      if(V(ord[i]) < 0) uL(l++) = ord[i];
      else uR(r++) = ord[i];

      balance(uL.head(l)).fill(-1.0/l * sqrt((double)l*r/(l+r)));
      balance(uR.head(r)).fill(+1.0/r * sqrt((double)l*r/(l+r)));

      double score = as_scalar(balance.t() * S * balance);
      //double score = fabs(dot(eigvec.tail_cols(1), balance));
      // Rcpp::Rcout << score << std::endl;
      //Rcpp::Rcout << balance.t();
      //Rcpp::Rcout << "Value:" << score <<std::endl;
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

void optimise_using_pc_forcing_branch(Balance<MaximumVariance>& balance, arma::mat& X, unsigned forced){
  MaximumVariance ebalance = MaximumVariance(balance.nodes, X);
  balance.setEvaluator(ebalance);
  if(balance.nodes.size() == 2){
    arma::uvec uL(1), uR(1);
    uL[0] = 0; uR[0] = 1;
    ebalance.eval(uL, uR, 1, 1);
    balance.set(uL, uR);
  }else{
    arma::mat Xsub = arma::mat(X.n_rows, balance.nodes.size());
    for(int i = 0; i < Xsub.n_cols; i++){
      Xsub.col(i) = X.col(balance.nodes[i][0]);
      for(int j = 1; j < balance.nodes[i].n_elem; j++){
        Xsub.col(i) += X.col(balance.nodes[i][j]);
      }
    }
    arma::mat U, V;
    arma::vec s;

    arma::svd_econ(U, s, V, clr_coordinates(Xsub));
    arma::uvec uF = arma::uvec(1);
    uF[0] = forced;
    balance.setWithLogContrastForceBranch(V.col(0), uF);

  }
}

void optimise_recursively(Balance<MaximumVariance>& balance, arma::mat& X, arma::mat& pb_mat, int *pb_size){

  optimise_balance_using_pc(balance, X);
  pb_mat.col(*pb_size) = balance.getBalance();
  (*pb_size)++;

  Balance<MaximumVariance> top = balance.top();
  if(top.nodes.size() > 1){
    optimise_recursively(top, X, pb_mat, pb_size);
  }
  Balance<MaximumVariance> left = balance.left();
  if(left.nodes.size() > 1){
    optimise_recursively(left, X, pb_mat, pb_size);
  }
  Balance<MaximumVariance> right  = balance.right();
  if(right.nodes.size() > 1){
    optimise_recursively(right, X, pb_mat, pb_size);
  }

}

void optimise_recursively_forcing_parents(Balance<MaximumVariance>& balance, arma::mat& X, arma::mat& pb_mat, int *pb_size){

  optimise_balance_using_pc(balance, X);
  pb_mat.col(*pb_size) = balance.getBalance();
  (*pb_size)++;

  Balance<MaximumVariance> top = balance.top();
  while(top.nodes.size() > 1){
    optimise_using_pc_forcing_branch(top, X, top.nodes.size()-1);
    pb_mat.col(*pb_size) = top.getBalance();
    (*pb_size)++;
    Balance<MaximumVariance> left = top.left();
    if(left.nodes.size() > 1){
      optimise_recursively_forcing_parents(left, X, pb_mat, pb_size);
    }
    Balance<MaximumVariance> right  = top.right();
    if(right.nodes.size() > 1){
      optimise_recursively_forcing_parents(right, X, pb_mat, pb_size);
    }
    top = top.top();
  }
  Balance<MaximumVariance> left = balance.left();
  if(left.nodes.size() > 1){
    optimise_recursively_forcing_parents(left, X, pb_mat, pb_size);
  }
  Balance<MaximumVariance> right  = balance.right();
  if(right.nodes.size() > 1){
    optimise_recursively_forcing_parents(right, X, pb_mat, pb_size);
  }

}


// [[Rcpp::export]]
arma::mat find_PB(arma::mat& X){
  int K = X.n_cols;

  arma::mat pb_mat = arma::mat(K,K-1);

  std::vector<Balance<MaximumVariance>> SOLS;

  Balance<MaximumVariance> balance = Balance<MaximumVariance>(X.n_cols);

  optimise(balance, X);

  SOLS.push_back(balance);

  for(int l=0;l<K-1;l++){
    //Rcpp::Rcout << "Starting step " << l + 1 << " of " << K-1 << std::endl;
    //print_list(SOLS);
    double vBestSolution = 0;
    int iBestSolution = -1;
    for(unsigned int i =0; i< SOLS.size();i++){
      double v = SOLS[i].eval();
      if(v > vBestSolution){
        vBestSolution = v;
        iBestSolution = i;
      }
    }

    pb_mat.col(l) = SOLS[iBestSolution].getBalance();

    Balance<MaximumVariance> top = SOLS[iBestSolution].top();
    Balance<MaximumVariance> left = SOLS[iBestSolution].left();
    Balance<MaximumVariance> right  = SOLS[iBestSolution].right();

    if(top.nodes.size() > 1){
      optimise(top, X);
      SOLS.push_back(top);
    }
    if(left.nodes.size() > 1){
      optimise(left, X);
      SOLS.push_back(left);
    }
    if(right.nodes.size() > 1){
      optimise(right, X);
      SOLS.push_back(right);
    }

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat find_PB_using_pc(arma::mat& X){
  int K = X.n_cols;

  arma::mat pb_mat = arma::mat(K,K-1);

  std::vector<Balance<MaximumVariance>> SOLS;

  Balance<MaximumVariance> balance = Balance<MaximumVariance>(X.n_cols);

  optimise_balance_using_pc(balance, X);

  SOLS.push_back(balance);

  for(int l=0;l<K-1;l++){
    //Rcpp::Rcout << "Starting step " << l + 1 << " of " << K-1 << std::endl;
    //print_list(SOLS);
    double vBestSolution = 0;
    int iBestSolution = -1;
    for(unsigned int i =0; i< SOLS.size();i++){
      double v = SOLS[i].eval();
      if(v > vBestSolution){
        vBestSolution = v;
        iBestSolution = i;
      }
    }

    pb_mat.col(l) = SOLS[iBestSolution].getBalance();

    Balance<MaximumVariance> top = SOLS[iBestSolution].top();
    Balance<MaximumVariance> left = SOLS[iBestSolution].left();
    Balance<MaximumVariance> right  = SOLS[iBestSolution].right();

    if(top.nodes.size() > 1){
      optimise_balance_using_pc(top, X);
      SOLS.push_back(top);
    }
    if(left.nodes.size() > 1){
      optimise_balance_using_pc(left, X);
      SOLS.push_back(left);
    }
    if(right.nodes.size() > 1){
      optimise_balance_using_pc(right, X);
      SOLS.push_back(right);
    }

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat find_PB_using_pc_recursively(arma::mat& X){
  int K = X.n_cols;

  arma::mat pb_mat = arma::zeros(K,K-1);

  std::vector<Balance<MaximumVariance>> SOLS;
  int pb_size = 0;
  Balance<MaximumVariance> root = Balance<MaximumVariance>(X.n_cols);
  optimise_recursively(root, X, pb_mat, &pb_size);


  // Rcpp::Rcout << pb_size << std::endl;
  return(pb_mat);
}

/*
 *
 *
 *
 */



// [[Rcpp::export]]
arma::mat find_PB_using_pc_recursively_forcing_parents(arma::mat& X){
  arma::mat pb_mat = arma::zeros(X.n_cols,X.n_cols-1);

  unsigned pb_i = 0;
  int K = X.n_cols;
  arma::mat S0 = cov(clr_coordinates(X));

  pb_mat.col(0) = get_balance_using_pc(X);
  // Rcpp::Rcout << pb_mat << std::endl;
  arma::uvec left = find(pb_mat.col(0) < 0);
  arma::uvec right = find(pb_mat.col(0) > 0);
  if(left.n_elem > 1){
    if(left.n_elem == 2){
      arma::vec balance = arma::zeros(K);
      balance[left(0)] = -SQR2DIV2;
      balance[left(1)] = +SQR2DIV2;
      pb_mat.col(++pb_i) = balance;
    }else{
      arma::mat Xsub = arma::mat(X.n_rows, left.n_elem);
      for(unsigned i = 0; i < left.n_elem; i++) Xsub.col(i) = X.col(left(i));

      unsigned new_cols = left.n_elem - 1;
      arma::ucolvec inew_cols = arma::ucolvec(new_cols);
      for(unsigned i = 0; i < new_cols; i++) inew_cols(i) = 1+pb_i + i;
      pb_mat.submat(left, inew_cols) = find_PB_using_pc_recursively_forcing_parents(Xsub);
      pb_i+=new_cols;
    }
  }

  if(right.n_elem > 1){
    if(right.n_elem == 2){
      arma::vec balance = arma::zeros(K);
      balance[right(0)] = -SQR2DIV2;
      balance[right(1)] = +SQR2DIV2;
      pb_mat.col(++pb_i) = balance;
    }else{
      arma::mat Xsub = arma::mat(X.n_rows, right.n_elem);
      for(unsigned i = 0; i < right.n_elem; i++) Xsub.col(i) = X.col(right(i));

      unsigned new_cols = right.n_elem - 1;
      arma::ucolvec inew_cols = arma::ucolvec(new_cols);
      for(unsigned i = 0; i < new_cols; i++) inew_cols(i) = 1+pb_i + i;
      pb_mat.submat(right, inew_cols) = find_PB_using_pc_recursively_forcing_parents(Xsub);
      pb_i+=new_cols;
    }
  }

  arma::uvec zeros = find( pb_mat.col(0) == 0 );
  arma::uvec no_zeros = find( pb_mat.col(0) != 0 );

  arma::vec bal = arma::zeros(K);


  while( zeros.n_elem > 0 ){
    double l = no_zeros.n_elem;
    unsigned r = 0, bestR;
    arma::uvec uR(K);
    if(zeros.n_elem <= 2){
      bestR = 1;
      bal(no_zeros).fill(-1.0/l * sqrt(l/(l+1)));
      bal(zeros.head(1)).fill(sqrt(l/(l+1)));  // When zeros.n_elem == 2 the decision is arbitrary.
    }else{
      arma::vec eigval;
      arma::mat eigvec;
      arma::mat Xsub = arma::mat(X.n_rows, zeros.n_elem);
      for(int i = 0; i < zeros.n_elem; i++) Xsub.col(i) = X.col(zeros[i]);
      Rcpp::Rcout << Xsub;
      arma::mat S = arma::cov(clr_coordinates(Xsub));
      eig_sym( eigval, eigvec, S);

      arma::vec V = eigvec.tail_cols(1);
      arma::uvec ord;
      if( fabs(min(V)) >= max(V) ){
        ord = sort_index(V, "ascend");
      }else{
        ord = sort_index(V, "descend");
      }

      double bestScore = -1;
      for(unsigned i = 0; arma::as_scalar(V(ord[0])) * arma::as_scalar(V(ord[i])) > 0; i++){
        uR(r++) = zeros[ord[i]];

        bal(no_zeros).fill(-1.0/l * sqrt((double)l*r/(l+r)));
        bal(uR.head(r)).fill(+1.0/r * sqrt((double)l*r/(l+r)));
        double score = as_scalar(bal.t() * S0 * bal);
        // double score = dot(V, bal);
        if(score > bestScore){
          bestScore = score;
          bestR = r;
        }
      }
      bal.fill(0);
      bal(no_zeros).fill(-1.0/l * sqrt((double)l*bestR/(l+bestR)));
      for(unsigned i = 0; i < bestR; i++){
        bal(uR(i)) = +1.0/bestR * sqrt((double)l*bestR/(l+bestR));
      }
    }
    zeros = find( bal == 0 );
    no_zeros = find( bal != 0 );

    pb_i++;
    pb_mat.col( pb_i ) = bal;

    right = uR.head(bestR);
    if(bestR > 1){
      if(bestR == 2){
        arma::vec balance = arma::zeros(K);
        balance[uR(0)] = -SQR2DIV2;
        balance[uR(1)] = +SQR2DIV2;
        pb_mat.col(++pb_i) = balance;
      }else{
        arma::mat Xsub = arma::mat(X.n_rows, right.n_elem);
        for(unsigned i = 0; i < right.n_elem; i++) Xsub.col(i) = X.col(right(i));

        unsigned new_cols = right.n_elem - 1;
        arma::ucolvec inew_cols = arma::ucolvec(new_cols);
        for(unsigned i = 0; i < new_cols; i++) inew_cols(i) = 1+pb_i + i;
        pb_mat.submat(right, inew_cols) = find_PB_using_pc_recursively_forcing_parents(Xsub);
        pb_i+=new_cols;
      }
    }

  }

  return(pb_mat);
}


/*** R
SEED = round(1000*runif(1))
set.seed(363)
SEED
D = 10
X = matrix(rlnorm(100*D), ncol = D)
sort(diag(cov(coordinates(X, find_PB(X)))), decreasing = TRUE)
sort(diag(cov(coordinates(X, find_PB_using_pc(X)))), decreasing = TRUE)
sort(diag(cov(coordinates(X, find_PB_using_pc_recursively(X)))), decreasing = TRUE)
sort(diag(cov(coordinates(X, find_PB_using_pc_recursively_forcing_parents(X)))), decreasing = TRUE)
*/

// #endif
