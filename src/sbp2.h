#include <RcppArmadillo.h>

class Balance {
  double X;

  double score = 0;
  std::map<int,arma::uvec> nodes;

  arma::uvec L;
  arma::uvec R;
  unsigned int L_length;
  unsigned int R_length;

  bool initialized = false;

public:
  Balance (std::map<int,arma::uvec> nodes0){
    nodes = nodes0;

    L = arma::uvec(nodes.size());
    L_length = 0;
    R = arma::uvec(nodes.size());
    R_length = 0;

  }
  void init(){
    arma::uvec L0 = {0}, R0 = {1};
    init(L0,R0);
  }
  void init(arma::uvec L0, arma::uvec R0){
    L_length = L0.n_elem;
    L.head(L0.n_elem) = L0;
    R_length = R0.n_elem;
    R.head(R0.n_elem) = R0;
    initialized = true;
  }
  arma::uvec getL(){ return(L.head(L_length)); }
  arma::uvec getR(){ return(R.head(R_length)); }

  bool hasNext(){
    int n = nodes.size();
    arma::uvec uL = getL();
    arma::uvec uR = getR();
    if( (uL.n_elem == 1) & (uR.n_elem + 1 == n) & (uL[0] == n-1) ){
      return(false);
    }
    return(true);
  }
  void nextBalance(){
    int n = nodes.size();
    arma::uvec uL = getL();
    arma::uvec uR = getR();
    arma::uvec O = arma::zeros<arma::uvec>(n);
    O(uL).fill(1);
    O(uR).fill(2);
    do{
      int pos = 0;
      while(O[pos] == 2){
        O[pos] = 0;
        pos++;
      }
      O[pos]++;
      uL = find(O == 1);
      uR = find(O == 2);
    } while ( (uL.n_elem == 0) | (uR.n_elem == 0) );
    init(uL,uR);
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

  void rnd_init(){
    int n = nodes.size();

    int first = (int)floor(n * arma::randu(1)[0]);
    int delta = 1 + (int)floor((n-1) * arma::randu(1)[0]);
    int second = (first + delta) % n;
    arma::vec v_initial(2);
    v_initial(0) = first;
    v_initial(1) = second;
    arma::vec v_random = arma::floor(3*arma::vec(n).randu());


    v_random[v_initial[0]] = 1;
    v_random[v_initial[1]] = 2;
    arma::uvec L0 = find(v_random == 1);
    arma::uvec R0 = find(v_random == 2);
    init(L0,R0);
  }

  void print(){
    if(!initialized){
      Rcpp::Rcout << "Balance not initilised";
      return;
    }
    Rcpp::Rcout << "L:";
    for(unsigned int i=0; i < L_length; i++){
      Rcpp::Rcout << " " << L[i];
    }
    Rcpp::Rcout << ", R:";
    for(unsigned int i=0; i < R_length; i++){
      Rcpp::Rcout << " " << R[i];
    }
    Rcpp::Rcout << std::endl;
  }
};

class EvaluateBalance {
protected:
  Balance *bal;

public:
  EvaluateBalance(Balance *bal_){
    bal = bal_;
  }
  virtual void print_state(){
    Rcpp::Rcout << "Print method for this heuristic not defined";
  }
};
