#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "sbp.h"
#include "principal_balances.h"
#include "coda.h"

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
arma::mat find_PB_rnd_local_search(arma::mat M, int rep = 1){
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

arma::vec get_pc1(arma::mat M){
  arma::vec b_old = arma::zeros<arma::vec>(M.n_cols);
  arma::vec b_new = arma::ones<arma::vec>(M.n_cols);
  int iter = 0;
  while( (max(abs(b_old-b_new)) > 10e-5) & (iter < 1000) ){
    iter++;
    b_old = b_new;
    b_new = normalise(M * b_old);
  }
  return(b_new);
}
arma::vec get_pc1_error(arma::mat M){
  int k = M.n_rows;

  arma::vec b = arma::ones<arma::vec>(k);
  arma::vec b_old = arma::zeros<arma::vec>(k);
  double mu_prev = 0, mu = 1;

  while( std::abs(mu - mu_prev) > 0.001 ){
    b = normalise(M * b);
    mu_prev = mu;
    mu = ((arma::mat)(b.t() * M * b))[0];
  }

  arma::mat inv_M = arma::mat(k, k);
  arma::mat MuId = arma::zeros(k, k);
  MuId.diag().fill(mu);

  while( max(abs(b_old-b)) > 10e-10){
    b_old = b;
    b = normalise( inv(M - MuId) * b );
    mu = ((arma::mat)(b.t() * M * b))(0,0);
    MuId.diag().fill( mu );
  }
  return(b);
}

// [[Rcpp::export]]
arma::mat find_PB_pc_local_search(arma::mat X){
  arma::mat M = cov(log(X));

  arma::mat Xclr = clr_coordinates(X);
  arma::mat Mclr = cov(Xclr);
  arma::vec PC1 = get_pc1(Mclr);

  int K = M.n_cols;

  std::vector<SBP> PB;
  arma::mat pb_mat = arma::mat(K,K-1);
  std::vector<SBP> SOLS;

  SOLS.push_back(SBP(M, default_node(K)));

  SOLS.back().first_pc_local_search(PC1);
  //print_list(SOLS);
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
    // Rcpp::Rcout << "Best solution found" << std::endl;
    // SOLS[iBestSolution].print_status(true,true,true);

    pb_mat.col(l) = balance(K,PB[l]);
    Xclr = Xclr - Xclr * pb_mat.col(l) * pb_mat.col(l).t();
    Mclr = cov(Xclr);

    arma::vec PC1 = get_pc1(Mclr);


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
      SOLS.back().first_pc_local_search(PC1);
      //SOLS.back().print_status(true,true,true);
      //Rcpp::Rcout << "End top" << std::endl;
    }
    if(nL > 1){
      //Rcpp:: Rcout << "Start left..." << std::endl;
      //SOLS[1].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].left());
      SOLS.back().first_pc_local_search(PC1);
      //Rcpp::Rcout << "End left" << std::endl;
    }
    if(nR > 1){
      //Rcpp:: Rcout << "Start right...";
      //SOLS[2].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].right());
      SOLS.back().first_pc_local_search(PC1);
      //Rcpp::Rcout << "End right" << std::endl;
    }
    //SOLS.back().print_status(true,true,true);

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat find_PB(arma::mat X){
  arma::mat M = cov(log(X));

  int K = M.n_cols;

  std::vector<SBP> PB;
  arma::mat pb_mat = arma::mat(K,K-1);
  std::vector<SBP> SOLS;

  SOLS.push_back(SBP(M, default_node(K)));

  SOLS.back().optimal();
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
    // Rcpp::Rcout << "Best solution found" << std::endl;
    // SOLS[iBestSolution].print_status(true,true,true);
    pb_mat.col(l) = balance(K,PB[l]);

    int n = SOLS[iBestSolution].get_n();
    int nL = SOLS[iBestSolution].getL().n_elem;
    int nR = SOLS[iBestSolution].getR().n_elem;
    if(n > nL + nR){
      //Rcpp:: Rcout << "Start top...";
      SOLS.push_back(SOLS[iBestSolution].top());
      SOLS.back().optimal();
      //SOLS.back().print_status(true,true,false);
      //Rcpp::Rcout << "End top" << std::endl;
    }
    if(nL > 1){
      //Rcpp:: Rcout << "Start left..." << std::endl;
      SOLS.push_back(SOLS[iBestSolution].left());
      SOLS.back().optimal();
      //SOLS.back().print_status(true,true,false);
      //Rcpp::Rcout << "End left" << std::endl;
    }
    if(nR > 1){
      //Rcpp:: Rcout << "Start right...";
      SOLS.push_back(SOLS[iBestSolution].right());
      SOLS.back().optimal();
      //SOLS.back().print_status(true,true,false);
      //Rcpp::Rcout << "End right" << std::endl;
    }

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return(pb_mat);
}

// [[Rcpp::export]]
arma::mat find_PB2(arma::mat M, int random = 100, int optim = 0){
  int K = M.n_cols;
  std::vector<SBP> PB;
  std::vector<SBP> SOLS;

  SOLS.push_back(SBP(M, default_node(K)));
  SOLS.back().simulated_annealing(random, optim);
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
      int random_scaled = (int)ceil(random * (double)(nL + nR)/K);
      int optim_scaled =  (int)ceil(optim * (double)(nL + nR)/K);
      SOLS.back().simulated_annealing(random_scaled , optim_scaled);
      //SOLS.back().print_status(true,true,true);
      //Rcpp::Rcout << "End top" << std::endl;
    }
    if(nL > 1){
      //Rcpp:: Rcout << "Start left..." << std::endl;
      //SOLS[1].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].left());
      int random_scaled = (int)ceil(random * (double)(nL)/K);
      int optim_scaled =  (int)ceil(optim * (double)(nL)/K);
      SOLS.back().simulated_annealing(random_scaled , optim_scaled);
      //Rcpp::Rcout << "End left" << std::endl;
    }
    if(nR > 1){
      //Rcpp:: Rcout << "Start right...";
      //SOLS[2].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].right());
      int random_scaled = (int)ceil(random * (double)(nR)/K);
      int optim_scaled =  (int)ceil(optim * (double)(nR)/K);
      SOLS.back().simulated_annealing(random_scaled , optim_scaled);
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
arma::mat find_PB3(arma::mat M, int steps, int random = 100, int optim = 0, int k = 0){
  int K = M.n_cols;
  if(k == 0) k = K;
  std::vector<SBP> PB;
  std::vector<SBP> SOLS;

  SOLS.push_back(SBP(M, default_node(K)));
  SOLS.back().simulated_annealing2(steps, random, optim, k);
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
      int random_scaled = (int)ceil(random * (double)(nL + nR)/K);
      int optim_scaled =  (int)ceil(optim * (double)(nL + nR)/K);
      SOLS.back().simulated_annealing2(steps, random_scaled , optim_scaled, k);
      //SOLS.back().print_status(true,true,true);
      //Rcpp::Rcout << "End top" << std::endl;
    }
    if(nL > 1){
      //Rcpp:: Rcout << "Start left..." << std::endl;
      //SOLS[1].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].left());
      int random_scaled = (int)ceil(random * (double)(nL)/K);
      int optim_scaled =  (int)ceil(optim * (double)(nL)/K);
      SOLS.back().simulated_annealing2(steps, random_scaled , optim_scaled, k);
      //Rcpp::Rcout << "End left" << std::endl;
    }
    if(nR > 1){
      //Rcpp:: Rcout << "Start right...";
      //SOLS[2].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].right());
      int random_scaled = (int)ceil(random * (double)(nR)/K);
      int optim_scaled =  (int)ceil(optim * (double)(nR)/K);
      SOLS.back().simulated_annealing2(steps, random_scaled , optim_scaled, k);
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
arma::mat find_PB4(arma::mat M){
  int K = M.n_cols;

  std::vector<SBP> PB;
  std::vector<SBP> SOLS;

  SOLS.push_back(SBP(M, default_node(K)));
  SOLS.back().first_component_approximation();
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
      SOLS.back().first_component_approximation();
      //SOLS.back().print_status(true,true,true);
      //Rcpp::Rcout << "End top" << std::endl;
    }
    if(nL > 1){
      //Rcpp:: Rcout << "Start left..." << std::endl;
      //SOLS[1].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].left());
      SOLS.back().first_component_approximation();
      //Rcpp::Rcout << "End left" << std::endl;
    }
    if(nR > 1){
      //Rcpp:: Rcout << "Start right...";
      //SOLS[2].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].right());
      SOLS.back().first_component_approximation();
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
arma::mat find_PB5(arma::mat X){
  arma::mat M = cov(log(X));

  arma::mat Xclr = clr_coordinates(X);
  arma::mat Mclr = cov(Xclr);

  int K = M.n_cols;

  std::vector<SBP> PB;
  arma::mat pb_mat = arma::mat(K,K-1);
  std::vector<SBP> SOLS;

  SOLS.push_back(SBP(M, default_node(K)));

  SOLS.back().first_component_approximation2(Mclr);
  //print_list(SOLS);
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
    // Rcpp::Rcout << "Best solution found" << std::endl;
    // SOLS[iBestSolution].print_status(true,true,true);
    pb_mat.col(l) = balance(K,PB[l]);
    Xclr = Xclr - Xclr * pb_mat.col(l) * pb_mat.col(l).t();
    Mclr = cov(Xclr);


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
      SOLS.back().first_component_approximation2(Mclr);
      //SOLS.back().print_status(true,true,true);
      //Rcpp::Rcout << "End top" << std::endl;
    }
    if(nL > 1){
      //Rcpp:: Rcout << "Start left..." << std::endl;
      //SOLS[1].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].left());
      SOLS.back().first_component_approximation2(Mclr);
      //Rcpp::Rcout << "End left" << std::endl;
    }
    if(nR > 1){
      //Rcpp:: Rcout << "Start right...";
      //SOLS[2].print_status(true,true,true);
      SOLS.push_back(SOLS[iBestSolution].right());
      SOLS.back().first_component_approximation2(Mclr);
      //Rcpp::Rcout << "End right" << std::endl;
    }
    //SOLS.back().print_status(true,true,true);

    SOLS[iBestSolution] = SOLS.back();
    SOLS.pop_back();

    Rcpp::checkUserInterrupt();
  }

  return(pb_mat);
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
