#include <RcppArmadillo.h>

class SBP {
  arma::mat M;

  arma::uvec L;
  arma::uvec R;

  unsigned int L_length;
  unsigned int R_length;

  std::map<int,arma::uvec> node;
  double sL, sR, sC, variance;

  bool initialized;

  void print_nodes();
  void print_LR();
  void print_M(){ Rcpp::Rcout << M; }

  void best_improve();
  void k_best_improve(int k);
public:
  double var(){ return(variance); }
  SBP (arma::mat, std::map<int,arma::uvec>);

  void init();
  void init(arma::uvec,arma::uvec);

  int get_nL();
  int get_nR();
  int get_n(){ return(M.n_cols); }
  int get_ni(int I){ return(node[I].n_elem); }

  arma::uvec get_indices(int I){ return(node[I]); }

  arma::uvec getL(){ return(L.head(L_length)); }
  arma::uvec getR(){ return(R.head(R_length)); }

  SBP top();
  SBP left();
  SBP right();


  void addL(int);
  double v_addL(int);
  void addR(int);
  double v_addR(int);
  void removeL(int);
  double v_removeL(int);
  void removeR(int);
  double v_removeR(int);

  void local_search(int rep);
  void simulated_annealing(int steps, int optim);
  void simulated_annealing2(int steps, int random, int optim, int k);
  void first_component_approximation();
  void first_component_approximation2(arma::mat Mclr);
  void first_component_approximation2(arma::vec Mclr);
  void print_status(bool, bool, bool);
};
