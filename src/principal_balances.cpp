// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "sbp.h"
#include "principal_balances.h"

void print_list(std::vector<SBP> SOLS){
  Rcpp::Rcout << ":Begin list:" << std::endl;
  for(unsigned int i=0; i<SOLS.size(); i++){
    Rcpp::Rcout << "Element " << i << std::endl;
    SOLS[i].print_status(true,true,true);
  }
  Rcpp::Rcout << ":End list:" << std::endl;
}
arma::vec balance(int K, SBP sbp){
  arma::uvec L = sbp.getL();
  arma::uvec R = sbp.getR();
  arma::uvec partition = arma::zeros<arma::uvec>(K);
  for(unsigned int i = 0; i< L.n_elem; i++){
    partition(sbp.get_indices(L[i])).fill(1);
  }
  for(unsigned int i = 0; i< R.n_elem; i++){
    partition(sbp.get_indices(R[i])).fill(2);
  }
  arma::uvec O = arma::zeros<arma::uvec>(K);

  arma::uvec iL = find(partition == 1);
  arma::uvec iR = find(partition == 2);

  double nL = (double)iL.n_elem;
  double nR = (double)iR.n_elem;

  arma::vec bal = arma::zeros(K);
  bal(iL).fill(1/nL * sqrt(nL*nR/(nL+nR)));
  bal(iR).fill(-1/nR * sqrt(nL*nR/(nL+nR)));
  return(bal);
}

// [[Rcpp::export]]
arma::mat find_PB(arma::mat M, int rep = 1){
  int K = M.n_cols;
  std::vector<SBP> PB;
  std::vector<SBP> SOLS;

  SOLS.push_back(SBP(M, default_node(K)));
  SOLS.back().local_search(rep);
  for(int l=0;l<K-1;l++){
    //Rcpp::Rcout << "Starting step " << l + 1 << " of " << K-1 << std::endl;
    //print_list(SOLS);
    double vBestSolution = 0;
    int iBestSolution = -1;
    for(unsigned int i =0; i< SOLS.size();i++){
      double v = SOLS[i].var();
      if(v > vBestSolution){
        vBestSolution = v;
        iBestSolution = i;
      }
    }
    PB.push_back(SOLS[iBestSolution]);
    //Rcpp::Rcout << "Best solution found" << std::endl;
    //SOLS[iBestSolution].print_status(true,true,true);

    //Rcpp::Rcout << "Status:" << SOLS.size() << std::endl;
    //SOLS.back().print_status(true,true,true);
    int n = SOLS[iBestSolution].get_n();
    int nL = SOLS[iBestSolution].getL().n_elem;
    int nR = SOLS[iBestSolution].getR().n_elem;
    //Rcpp::Rcout << n << " " << nL << " " << nR << std::endl;
    if(n > nL + nR){
      //Rcpp:: Rcout << "Start top...";
      //Rcpp::Rcout << "Including top... ";
      SOLS.push_back(SOLS[iBestSolution].top());
      //Rcpp::Rcout << "Top included!" << std::endl;
      //SOLS.back().print_status(true,true,true);
      SOLS.back().local_search(rep);
      //SOLS.back().print_status(true,true,true);
      //Rcpp::Rcout << "End top" << std::endl;
    }
    if(nL > 1){
      //Rcpp:: Rcout << "Start left..." << std::endl;
      //SOLS[1].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].left());
      SOLS.back().local_search(rep);
      //Rcpp::Rcout << "End left" << std::endl;
    }
    if(nR > 1){
      //Rcpp:: Rcout << "Start right...";
      //SOLS[2].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].right());
      SOLS.back().local_search(rep);
      //Rcpp::Rcout << "End right" << std::endl;
    }
    //SOLS.back().print_status(true,true,true);

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();
  }
  //Rcpp::Rcout << "End" << std::endl;
  arma::mat pb_mat = arma::mat(K,PB.size());
  for(unsigned int i=0; i<PB.size(); i++){
    pb_mat.col(i) = balance(K,PB[i]);
  }
  return(pb_mat);
  //Rcpp::Rcout << "Principal balances: "<< std::endl;
  //print_list(PB);
}

// [[Rcpp::export]]
arma::mat find_PB2(arma::mat M, int steps = 100){
  int K = M.n_cols;
  std::vector<SBP> PB;
  std::vector<SBP> SOLS;

  SOLS.push_back(SBP(M, default_node(K)));
  SOLS.back().simulated_annealing(steps);
  for(int l=0;l<K-1;l++){
    //Rcpp::Rcout << "Starting step " << l + 1 << " of " << K-1 << std::endl;
    //print_list(SOLS);
    double vBestSolution = 0;
    int iBestSolution = -1;
    for(unsigned int i =0; i< SOLS.size();i++){
      double v = SOLS[i].var();
      if(v > vBestSolution){
        vBestSolution = v;
        iBestSolution = i;
      }
    }
    PB.push_back(SOLS[iBestSolution]);
    //Rcpp::Rcout << "Best solution found" << std::endl;
    //SOLS[iBestSolution].print_status(true,true,true);

    //Rcpp::Rcout << "Status:" << SOLS.size() << std::endl;
    //SOLS.back().print_status(true,true,true);
    int n = SOLS[iBestSolution].get_n();
    int nL = SOLS[iBestSolution].getL().n_elem;
    int nR = SOLS[iBestSolution].getR().n_elem;
    //Rcpp::Rcout << n << " " << nL << " " << nR << std::endl;
    if(n > nL + nR){
      //Rcpp:: Rcout << "Start top...";
      //Rcpp::Rcout << "Including top... ";
      SOLS.push_back(SOLS[iBestSolution].top());
      //Rcpp::Rcout << "Top included!" << std::endl;
      //SOLS.back().print_status(true,true,true);
      SOLS.back().simulated_annealing(steps);
      //SOLS.back().print_status(true,true,true);
      //Rcpp::Rcout << "End top" << std::endl;
    }
    if(nL > 1){
      //Rcpp:: Rcout << "Start left..." << std::endl;
      //SOLS[1].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].left());
      SOLS.back().simulated_annealing(steps);
      //Rcpp::Rcout << "End left" << std::endl;
    }
    if(nR > 1){
      //Rcpp:: Rcout << "Start right...";
      //SOLS[2].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].right());
      SOLS.back().simulated_annealing(steps);
      //Rcpp::Rcout << "End right" << std::endl;
    }
    //SOLS.back().print_status(true,true,true);

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();
  }
  //Rcpp::Rcout << "End" << std::endl;
  arma::mat pb_mat = arma::mat(K,PB.size());
  for(unsigned int i=0; i<PB.size(); i++){
    pb_mat.col(i) = balance(K,PB[i]);
  }
  return(pb_mat);
  //Rcpp::Rcout << "Principal balances: "<< std::endl;
  //print_list(PB);
}

// Helpers
std::map<int,arma::uvec> default_node(int size){
  std::map<int,arma::uvec> node;
  for(int i=0;i<size;i++){
    node[i] = arma::uvec(1);
    node[i][0] = i;
  }
  return(node);
}
