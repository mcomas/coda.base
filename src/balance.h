#include <RcppArmadillo.h>

class Balance {

  arma::uvec L;
  arma::uvec R;

public:
  std::map<int,arma::uvec> nodes;
  unsigned n;
  unsigned int L_length;
  unsigned int R_length;
  Balance (std::map<int,arma::uvec> nodes0){
    nodes = nodes0;
    n = nodes.size();

    L = arma::uvec(n);
    L_length = 0;
    R = arma::uvec(n);
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
    return(!( (uL.n_elem == 1) & (uR.n_elem + 1 == n) & (uL[0] == n-1) ));
  }
  void nextBalance(){
    arma::uvec uL = L.head(L_length);
    arma::uvec uR = R.head(R_length);
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
  double eval(){
    return eval(bal);
  }
  virtual double eval(Balance *bal){
    return 0;
  }
  void print_state(){
    bal->print();
    // Rcpp::Rcout << logX.cols(bal->getL()) << std::endl;
    // Rcpp::Rcout << logX.cols(bal->getR());
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
    bal->init();
    //print_state();

    double best_score = eval();
    arma::uvec best_L = arma::uvec(bal->getL());
    arma::uvec best_R = arma::uvec(bal->getR());


    unsigned int iter = 1;
    while(bal->hasNext()){
      if(iter % 10000 == 0){
        R_CheckUserInterrupt();
      }
      iter++;
      bal->nextBalance();
      //print_state();

      double score = eval();
      if(score > best_score){
        best_score = score;
        best_L = arma::uvec(bal->getL());
        best_R = arma::uvec(bal->getR());
      }
    }

    bal->init(best_L, best_R);
    return best_score;
    //print_state();
  }
  double v_addL(int I){

    arma::uvec uLnew = arma::uvec(1+bal->L_length);
    uLnew.head(bal->L_length) = bal->getL();
    uLnew(bal->L_length) = I;
    arma::uvec uRnew = arma::uvec(bal->getR());

    Balance bnew = Balance(bal->nodes);
    bnew.init(uLnew, uRnew);
    return eval(&bnew);

  }
  double v_addR(int I){

    arma::uvec uLnew = arma::uvec(bal->getL());
    arma::uvec uRnew = arma::uvec(1+bal->R_length);
    uRnew.head(bal->R_length) = bal->getR();
    uRnew(bal->R_length) = I;


    Balance bnew = Balance(bal->nodes);
    bnew.init(uLnew, uRnew);
    return eval(&bnew);

  }
  double v_removeL(int I){

    arma::uvec uLnew = find(bal->getL() != I);
    arma::uvec uRnew = arma::uvec(bal->getR());

    Balance bnew = Balance(bal->nodes);
    bnew.init(uLnew, uRnew);
    return eval(&bnew);

  }
  double v_removeR(int I){

    arma::uvec uLnew = arma::uvec(bal->getL());
    arma::uvec uRnew = find(bal->getR() != I);

    Balance bnew = Balance(bal->nodes);
    bnew.init(uLnew, uRnew);
    double score = eval(&bnew);
    return score;

  }
  bool linear_transformation(arma::mat M){

  }
  bool best_improve(){
    // Rcpp::Rcout << eval() << std::endl;
    // int n = get_n();

    arma::uvec best_L = arma::uvec(bal->getL());
    arma::uvec best_R = arma::uvec(bal->getR());

    arma::uvec uL = arma::uvec(bal->getL());
    arma::uvec uR = arma::uvec(bal->getR());


    arma::uvec O = arma::zeros<arma::uvec>(bal->n);
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
          bal->addL((unsigned int)which);
        }else{
          bal->addR((unsigned int)which);
        }
      }else{ // remove
        if(left){
          bal->removeL((unsigned int)which);
        }else{
          bal->removeR((unsigned int)which);
        }
      }
      return true;
    }
    return false;
  }

};
