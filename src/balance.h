#ifndef balance_H
#define balance_H

#include <RcppArmadillo.h>

class Balance {
  arma::uvec L;
  arma::uvec R;

public:
  unsigned int L_length;
  unsigned int R_length;
  int N;
  std::map<int,arma::uvec> nodes;

  Balance(){}
  Balance(std::map<int,arma::uvec> nodes0){
    nodes = nodes0;
    N = nodes.size();
    L = arma::uvec(N);
    L_length = 0;
    R = arma::uvec(N);
    R_length = 0;
  }
  // Initialiasation ready to iterate all elements sequentially
  void init(){
    arma::uvec L0 = {0}, R0 = {1};
    init(L0,R0);
  }
  void init(arma::uvec L0, arma::uvec R0){
    L_length = L0.n_elem;
    L.head(L0.n_elem) = L0;
    R_length = R0.n_elem;
    R.head(R0.n_elem) = R0;
  }
  std::map<int,arma::uvec> top(){
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
    return(node0);
  }
  std::map<int,arma::uvec> left(){
    arma::uvec uL = getL();
    arma::uvec uR = getL();
    std::map<int,arma::uvec> node0;
    for(unsigned int i = 0; i<uL.n_elem; i++){
      node0[i] = arma::uvec(nodes[uL[i]]);
    }
    return(node0);
  }
  std::map<int,arma::uvec> right(){
    arma::uvec uL = getR();
    arma::uvec uR = getR();
    std::map<int,arma::uvec> node0;
    for(unsigned int i = 0; i<uR.n_elem; i++){
      node0[i] = arma::uvec(nodes[uR[i]]);
    }
    return(node0);
  }
  arma::uvec getL(){ return(L.head(L_length)); }
  arma::uvec getR(){ return(R.head(R_length)); }
  void addL(unsigned int I){
    L(L_length) = I;
    L_length++;
  }
  void addR(int I){
    R(R_length) = I;
    R_length++;
  }
  void removeL(int I){
    arma::uvec uLnew = find(L.head(L_length) != I);
    L_length--;
    L.head(L_length) = uLnew;
  }
  void removeR(int I){
    arma::uvec uRnew = find(R.head(R_length) != I);
    R_length--;
    R.head(R_length) = uRnew;
  }
  bool hasNext(){
    arma::uvec uL = L.head(L_length);
    arma::uvec uR = R.head(R_length);
    return(!( (uL.n_elem == 1) & (uR.n_elem + 1 == N) & (uL[0] == N-1) ));
  }
  void nextBalance(){
    arma::uvec uL = L.head(L_length);
    arma::uvec uR = R.head(R_length);
    arma::uvec O = arma::zeros<arma::uvec>(N);
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
  int get_n(){
    double n = 0;
    for(unsigned int i=0; i<N;n+=nodes[i++].n_elem);
    return(n);
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

    int first = (int)floor(N * arma::randu(1)[0]);
    int delta = 1 + (int)floor((N-1) * arma::randu(1)[0]);
    int second = (first + delta) % N;
    arma::vec v_initial(2);
    v_initial(0) = first;
    v_initial(1) = second;
    arma::vec v_random = arma::floor(3*arma::vec(N).randu());


    v_random[v_initial[0]] = 1;
    v_random[v_initial[1]] = 2;
    arma::uvec L0 = find(v_random == 1);
    arma::uvec R0 = find(v_random == 2);
    init(L0,R0);
  }

  void print(){
    Rcpp::Rcout << "Elements: ";
    for(unsigned int i=0; i<nodes.size();i++){
      Rcpp::Rcout << "{";
      for(int j = 0; j<nodes[i].n_elem;j++) Rcpp::Rcout << " " << nodes[i][j];
      Rcpp::Rcout << " } ";
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
  int D;
public:
  Balance bal;
  EvaluateBalance(){};
  EvaluateBalance(int D_){
    std::map<int,arma::uvec> nodes0;
    for(int i=0;i<D_;i++){
      nodes0[i] = arma::uvec(1);
      nodes0[i][0] = i;
    }
    bal = Balance(nodes0);
    D = D_;
  }
  EvaluateBalance(std::map<int,arma::uvec> nodes0, int D_){
    bal = Balance(nodes0);
    D = D_;
  }
  double eval(){
    return eval(&bal);
  }
  virtual double eval(Balance *bal){
    Rcpp::Rcout << "Calling parent, return 0" << std::endl;
    return 0;
  }
  arma::vec getBalance(Balance *bal_){
    arma::uvec L = bal_->getL();
    arma::uvec R = bal_->getR();

    arma::uvec partition = arma::zeros<arma::uvec>(D);
    for(unsigned int i = 0; i< bal_->L_length; i++){
      partition(bal_->nodes[L[i]]).fill(1);
    }
    for(unsigned int i = 0; i< bal_->R_length; i++){
      partition(bal_->nodes[R[i]]).fill(2);
    }
    arma::uvec iL = find(partition == 1);
    arma::uvec iR = find(partition == 2);

    double nL = (double)iL.n_elem;
    double nR = (double)iR.n_elem;

    arma::vec b = arma::zeros(D);
    b(iL).fill(1/nL * sqrt(nL*nR/(nL+nR)));
    b(iR).fill(-1/nR * sqrt(nL*nR/(nL+nR)));
    return(b);
  }
  arma::vec getBalance(){
    arma::uvec L = bal.getL();
    arma::uvec R = bal.getR();

    arma::uvec partition = arma::zeros<arma::uvec>(D);
    for(unsigned int i = 0; i< bal.L_length; i++){
      partition(bal.nodes[L[i]]).fill(1);
    }
    for(unsigned int i = 0; i< bal.R_length; i++){
      partition(bal.nodes[R[i]]).fill(2);
    }
    arma::uvec iL = find(partition == 1);
    arma::uvec iR = find(partition == 2);

    double nL = (double)iL.n_elem;
    double nR = (double)iR.n_elem;

    arma::vec b = arma::zeros(D);
    b(iL).fill(1/nL * sqrt(nL*nR/(nL+nR)));
    b(iR).fill(-1/nR * sqrt(nL*nR/(nL+nR)));
    return(b);
  }
  void print_state(){
    bal.print();
    Rcpp::Rcout << eval() << std::endl;
  }
  double setLocalSearch(){ //arma::vec v
    int iter = 0;
    while(best_improve()){
      Rcpp::checkUserInterrupt();
    }
    return eval();
  }
  double setOptimal(){
    bal.init();

    double best_score = eval();
    arma::uvec best_L = arma::uvec(bal.getL());
    arma::uvec best_R = arma::uvec(bal.getR());


    unsigned int iter = 1;
    while(bal.hasNext()){
      if(iter % 10000 == 0){
        R_CheckUserInterrupt();
      }
      iter++;
      bal.nextBalance();

      double score = eval();
      if(score > best_score){
        best_score = score;
        best_L = arma::uvec(bal.getL());
        best_R = arma::uvec(bal.getR());
      }
    }

    bal.init(best_L, best_R);
    return best_score;
  }
  double v_addL(int I){

    arma::uvec uLnew = arma::uvec(1+bal.L_length);
    uLnew.head(bal.L_length) = bal.getL();
    uLnew(bal.L_length) = I;
    arma::uvec uRnew = arma::uvec(bal.getR());

    Balance bnew = Balance(bal.nodes);
    bnew.init(uLnew, uRnew);
    return eval(&bnew);

  }
  double v_addR(int I){

    arma::uvec uLnew = arma::uvec(bal.getL());
    arma::uvec uRnew = arma::uvec(1+bal.R_length);
    uRnew.head(bal.R_length) = bal.getR();
    uRnew(bal.R_length) = I;


    Balance bnew = Balance(bal.nodes);
    bnew.init(uLnew, uRnew);
    return eval(&bnew);

  }
  double v_removeL(int I){

    arma::uvec uLnew = find(bal.getL() != I);
    arma::uvec uRnew = arma::uvec(bal.getR());

    Balance bnew = Balance(bal.nodes);
    bnew.init(uLnew, uRnew);
    return eval(&bnew);

  }
  double v_removeR(int I){

    arma::uvec uLnew = arma::uvec(bal.getL());
    arma::uvec uRnew = find(bal.getR() != I);

    Balance bnew = Balance(bal.nodes);
    bnew.init(uLnew, uRnew);
    double score = eval(&bnew);
    return score;

  }
  bool linear_transformation(arma::mat M){
    return true;
  }
  bool best_improve(){

    arma::uvec best_L = arma::uvec(bal.getL());
    arma::uvec best_R = arma::uvec(bal.getR());

    arma::uvec uL = arma::uvec(bal.getL());
    arma::uvec uR = arma::uvec(bal.getR());


    arma::uvec O = arma::zeros<arma::uvec>(bal.N);
    O(uL).fill(1);
    O(uR).fill(1);


    arma::uvec uO = find(O == 0);

    double max_score = eval();
    //Rcpp::Rcout << "Initial score: " << max_score << std::endl;

    bool add = false;
    bool left = false;
    int which = -1;

    for(unsigned int i=0; i < uO.n_elem; i++){
      double new_score = v_addL(uO[i]);
      if(new_score > max_score){
        max_score = new_score;
        add = true;
        left = true;
        which = uO[i];
      }
      new_score = v_addR(uO[i]);
      if(new_score > max_score){
        max_score = new_score;
        add = true;
        left = false;
        which = uO[i];
      }
    }
    if(uL.n_elem > 1){
      for(unsigned int i=0; i<uL.n_elem; i++){
        double new_score = v_removeL(uL[i]);
        if(new_score > max_score){
          max_score = new_score;
          add = false;
          left = true;
          which = uL[i];
        }
      }
    }
    if(uR.n_elem > 1){
      for(unsigned int i=0; i<uR.n_elem; i++){
        double new_score = v_removeR(uR[i]);
        if(new_score > max_score){
          max_score = new_score;
          add = false;
          left = false;
          which = uR[i];
        }
      }
    }
    //Rcpp::Rcout << "Add: " << add << " Left: " << left << " I: " << which << std::endl;
    //Rcpp::Rcout << "Best score: " << max_score << std::endl;
    if(which > -1){
      if(add){
        if(left){
          bal.addL((unsigned int)which);
        }else{
          bal.addR((unsigned int)which);
        }
      }else{ // remove
        if(left){
          bal.removeL((unsigned int)which);
        }else{
          bal.removeR((unsigned int)which);
        }
      }
      return true;
    }
    return false;
  }

};

#endif
